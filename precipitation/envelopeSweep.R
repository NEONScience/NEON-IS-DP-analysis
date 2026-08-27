source('~/github/NEON-IS-data-processing/flow/flow.precip.aepg.smooth/def.precip.depth.smooth.R')

# =============================================================================
# envelopeSweep.R
# Runs the dynamic-envelope derivation (same logic as envelopeChecks.R, minus
# the smoothing/surrogate steps) across every site x every sampled date under
# DataBaseDir, to see which sites' configured Envelope still matches what the
# data actually supports, and which ones swing around a lot (unstable).
# =============================================================================

DataBaseDir <- '~/pfs/precipWeighing_thresh_select_ts_pad_smoother'
ThreshFile  <- file.path(DataBaseDir,'thresholds.csv')
OutCsv      <- '~/pfs/envelopeSweep_results.csv'

## ---- Helpers (same logic as envelopeChecks.R) -----------------------------
def.dir.in.site.date <- function(baseDir,date,site){
  siteDir <- Sys.glob(file.path(path.expand(baseDir),date,paste0('precip-weighing_',site,'*')))
  if(length(siteDir) != 1) return(NULL)
  sensDir <- list.dirs(siteDir,recursive=FALSE)
  if(length(sensDir) != 1) return(NULL)
  cfgDir <- list.dirs(sensDir,recursive=FALSE)
  if(length(cfgDir) != 1) return(NULL)
  return(cfgDir)
}

def.thresh.site <- function(threshFile,site){
  thresh <- utils::read.csv(path.expand(threshFile),stringsAsFactors=FALSE)
  row <- thresh[thresh$Named.Location.Name == site,]
  if(nrow(row) == 0) return(NULL)
  if(nrow(row) > 1) row <- row[1,]
  list(WndwAgr=row$WndwAgr,RangeSizeHour=row$RangeSizeHour,Envelope=row$Envelope,Recharge=row$Recharge)
}

# Runs the dynamic-envelope derivation for one site/date. Returns one summary row.
def.envelope.check <- function(site,date){
  out <- list(Site=site,Date=date,EnvelopeConfig=NA,EnvelopeDerived=NA,RechargeConfig=NA,
              RechargeDerived=NA,nDayNoRain=NA,nDayTotal=NA,Error=NA)

  thresh <- def.thresh.site(ThreshFile,site)
  if(is.null(thresh)){ out$Error <- 'no threshold row'; return(out) }
  out$EnvelopeConfig <- thresh$Envelope
  out$RechargeConfig <- thresh$Recharge

  DirIn <- def.dir.in.site.date(DataBaseDir,date,site)
  if(is.null(DirIn)){ out$Error <- 'no data dir'; return(out) }

  res <- base::tryCatch({
    Envelope <- thresh$Envelope
    Recharge <- thresh$Recharge
    WndwAgr <- thresh$WndwAgr
    WndwAgrNumc <- as.numeric(stringr::str_extract(string=WndwAgr,pattern='[0-9]+'))
    rangeSize <- if(stringr::str_detect(WndwAgr,'min')) thresh$RangeSizeHour*(60/WndwAgrNumc) else thresh$RangeSizeHour/WndwAgrNumc

    log <- NEONprocIS.base::def.log.init(Lvl='error')
    dirInData <- fs::path(DirIn,'data')
    fileData <- base::list.files(dirInData,pattern='.parquet',full.names=FALSE)
    fileManifest <- base::list.files(dirInData,pattern='manifest',full.names=TRUE)
    if(length(fileManifest) == 0) stop('no manifest')

    data <- NEONprocIS.base::def.read.parq.ds(fileIn=fs::path(dirInData,fileData),
                                              VarTime='readout_time',RmvDupl=TRUE,Df=TRUE,log=log)

    data$strain_gauge1_stability[is.na(data$strainGauge1Depth)] <- NA
    data$strain_gauge2_stability[is.na(data$strainGauge2Depth)] <- NA
    data$strain_gauge3_stability[is.na(data$strainGauge3Depth)] <- NA
    data <- data %>% dplyr::mutate(strainGaugeDepth = base::rowMeans(x=base::cbind(strainGauge1Depth,strainGauge2Depth,strainGauge3Depth),na.rm=F),
                                   strainGaugeStability = base::rowSums(x=base::cbind(strain_gauge1_stability,strain_gauge2_stability,strain_gauge3_stability),na.rm=F)==3)
    data$strainGaugeDepth[data$strainGaugeStability == FALSE] <- NA

    strainGaugeDepthAgr <- data %>%
      dplyr::mutate(startDateTime = lubridate::floor_date(as.POSIXct(readout_time,tz='UTC'),unit=WndwAgr)) %>%
      dplyr::group_by(startDateTime) %>%
      dplyr::summarise(strainGaugeDepth = mean(strainGaugeDepth,na.rm=T),.groups='drop')

    dataHourly <- strainGaugeDepthAgr %>%
      dplyr::mutate(startDateTime = lubridate::floor_date(startDateTime,unit='hour')) %>%
      dplyr::group_by(startDateTime) %>%
      dplyr::summarise(strainGaugeDepthMin = min(strainGaugeDepth,na.rm=T),
                       strainGaugeDepthMax = max(strainGaugeDepth,na.rm=T),.groups='drop')
    dataDaily <- dataHourly %>%
      dplyr::mutate(startDateTime = lubridate::floor_date(startDateTime,unit='day')) %>%
      dplyr::group_by(startDateTime) %>%
      dplyr::summarise(strainGaugeDepthChg = tail(strainGaugeDepthMax,1)-head(strainGaugeDepthMin,1),.groups='drop')

    setNoRain <- (dataDaily$strainGaugeDepthChg < 0.25*Envelope) & (dataDaily$strainGaugeDepthChg > -1*Envelope)
    setNoRain[is.na(setNoRain)] <- FALSE
    out$nDayTotal <- nrow(dataDaily)
    out$nDayNoRain <- sum(setNoRain)

    if(any(setNoRain)){
      dayNoRain <- dataDaily$startDateTime[setNoRain]
      timeDay <- lubridate::floor_date(strainGaugeDepthAgr$startDateTime,unit='day')
      envelopeComp <- data.frame(day=unique(timeDay),envelope=as.numeric(NA))
      for(idxDay in unique(timeDay)){
        if(!(idxDay %in% dayNoRain)) next
        setDay <- timeDay == idxDay
        envelopeComp$envelope[envelopeComp$day == idxDay] <-
          max(strainGaugeDepthAgr$strainGaugeDepth[setDay],na.rm=TRUE) - min(strainGaugeDepthAgr$strainGaugeDepth[setDay],na.rm=TRUE)
      }
      envelopeMax <- max(envelopeComp$envelope,na.rm=TRUE)
      if(!is.na(envelopeMax)){
        out$EnvelopeDerived <- envelopeMax
        out$RechargeDerived <- if(envelopeMax > Recharge/3) 3*envelopeMax else Recharge
      }
    }
    out
  },error=function(e){ out$Error <<- conditionMessage(e); out })

  return(res)
}

## ---- Sweep every site x every sampled date --------------------------------
sites <- unique(utils::read.csv(path.expand(ThreshFile),stringsAsFactors=FALSE)$Named.Location.Name)
dates <- c('2025/07/01','2025/08/15','2025/09/15','2025/10/15','2025/11/15','2025/12/15',
          '2026/01/15','2026/02/15','2026/03/15','2026/04/15','2026/05/15')

results <- do.call(rbind,lapply(sites,function(site){
  do.call(rbind,lapply(dates,function(date){
    message(site,' ',date)
    as.data.frame(def.envelope.check(site,date),stringsAsFactors=FALSE)
  }))
}))

utils::write.csv(results,path.expand(OutCsv),row.names=FALSE)
message('Wrote ',nrow(results),' rows to ',OutCsv)
