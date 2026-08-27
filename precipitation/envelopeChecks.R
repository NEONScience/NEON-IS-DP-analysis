rm(list=setdiff(ls(),c('arg','log')))
source('~/github/NEON-IS-data-processing/flow/flow.precip.aepg.smooth/def.ucrt.agr.precip.bench.R')
source('~/github/NEON-IS-data-processing/flow/flow.precip.aepg.smooth/def.precip.depth.smooth.R')

# =============================================================================
# envelopeChecks.R
# Plug-and-play driver to sanity-check precip smoothing thresholds against a
# sampled site-day. Defaults for all tunable thresholds are pulled from
# thresholds.csv; override any of them in the block below to test a change.
# =============================================================================

## ---- USER INPUT: pick site + date (must exist under DataBaseDir) ---------
Site <- 'WREF'                     # 4-letter site code
Date <- '2024/10/15'               # yyyy/mm/dd
ShowSurrogates <- FALSE            # TRUE also plots the 30 uncertainty surrogate depth/bench traces

DataBaseDir <- '~/pfs/precipWeighing_thresh_select_ts_pad_smoother'
ThreshFile  <- file.path(DataBaseDir,'thresholds.csv')
DirOutBase  <- '~/pfs/out_envelopeChecks'

## ---- Locate the input CFGLOC directory for the chosen site/date ----------
def.dir.in.site.date <- function(baseDir,date,site){
  siteDir <- Sys.glob(file.path(path.expand(baseDir),date,paste0('precip-weighing_',site,'*')))
  if(length(siteDir) != 1){
    stop(paste0('Expected exactly 1 site dir for ',site,' on ',date,', found ',length(siteDir)))
  }
  sensDir <- list.dirs(siteDir,recursive=FALSE) # aepg600m or aepg600m_heated
  if(length(sensDir) != 1){
    stop(paste0('Expected exactly 1 sensor-type dir under ',siteDir,', found ',length(sensDir)))
  }
  cfgDir <- list.dirs(sensDir,recursive=FALSE) # CFGLOC#####
  if(length(cfgDir) != 1){
    stop(paste0('Expected exactly 1 CFGLOC dir under ',sensDir,', found ',length(cfgDir)))
  }
  return(cfgDir)
}

DirIn <- def.dir.in.site.date(DataBaseDir,Date,Site)
message('DirIn: ',DirIn)

## ---- Pull default thresholds for this site from thresholds.csv -----------
def.thresh.site <- function(threshFile,site){
  thresh <- utils::read.csv(path.expand(threshFile),stringsAsFactors=FALSE)
  row <- thresh[thresh$Named.Location.Name == site,]
  if(nrow(row) == 0){
    stop(paste0('No threshold row found for site ',site,' in ',threshFile))
  }
  if(nrow(row) > 1){
    warning(paste0(nrow(row),' threshold rows found for site ',site,' - using the first'))
    row <- row[1,]
  }
  list(
    ThshCountHour    = row$ThshCountHour,
    WndwAgr          = row$WndwAgr,
    ChangeFactor     = row$ChangeFactor,
    ExtremePrecipMax = row$ExtremePrecipMax,
    Quant            = row$Quant,
    RangeSizeHour    = row$RangeSizeHour,
    Envelope         = row$Envelope,
    ChangeFactorEvap = row$ChangeFactorEvap,
    Recharge         = row$Recharge,
    ThshChange       = row$ThshChange
  )
}

thresh <- def.thresh.site(ThreshFile,Site)

## ---- Override any threshold here to test a change (defaults come from the CSV above) ----
thresh$WndwAgr          <- '5 min'
#thresh$RangeSizeHour    <- 60
# thresh$ThshCountHour    <- 15
# thresh$Quant            <- 0.5
# thresh$ThshChange       <- 0.2
# thresh$ChangeFactor     <- 1
# thresh$ChangeFactorEvap <- 0.5
# thresh$Recharge         <- 20
# thresh$ExtremePrecipMax <- 60
# thresh$Envelope         <- 2

# unpack into the variable names the processing code below expects
ThshCountHour    <- thresh$ThshCountHour
WndwAgr          <- thresh$WndwAgr
ChangeFactor     <- thresh$ChangeFactor
ExtremePrecipMax <- thresh$ExtremePrecipMax
Quant            <- thresh$Quant
RangeSizeHour    <- thresh$RangeSizeHour
Envelope         <- thresh$Envelope
ChangeFactorEvap <- thresh$ChangeFactorEvap
Recharge         <- thresh$Recharge
ThshChange       <- thresh$ThshChange

DirSubCopy <- NULL

# Start logging
log <- NEONprocIS.base::def.log.init(Lvl='debug')

# Gather info about the input directory and create the output directory.
InfoDirIn <- NEONprocIS.base::def.dir.splt.pach.time(DirIn,log=log)
dirInData <- fs::path(DirIn,'data')
dirInFlags <- fs::path(DirIn,'flags')
dirOut <- fs::path(DirOutBase,InfoDirIn$dirRepo)
dirOutData <- fs::path(dirOut,'data')
dirOutFlags <- fs::path(dirOut,'flags')
NEONprocIS.base::def.dir.crea(DirBgn = dirOut,
                              DirSub = c('data','flags'),
                              log = log)

# Copy with a symbolic link the desired subfolders
DirSubCopy <- base::unique(base::setdiff(DirSubCopy,c('data','flags')))
if(base::length(DirSubCopy) > 0){

  NEONprocIS.base::def.dir.copy.symb(DirSrc=fs::path(DirIn,DirSubCopy),
                                     DirDest=dirOut,
                                     LnkSubObj=FALSE,
                                     log=log)
}

#adjust thresholds based on WndwAgr unit
WndwAgrNumc <- as.numeric(stringr::str_extract(string = WndwAgr, pattern = '[0-9]+'))
if(stringr::str_detect(WndwAgr, 'min')) {
  ThshCount <- ThshCountHour * (60/WndwAgrNumc)
  rangeSize <- RangeSizeHour*(60/WndwAgrNumc)   #!!! POTENTIAL FOR MAKING AN INPUT VARIABLE !!!
  if(WndwAgrNumc > 60 | WndwAgrNumc < 5){
    log$fatal('averaging unit must be between 5 minutes and one hour')
    stop()
  }
} else if ((stringr::str_detect(WndwAgr, 'hour')) ){
  ThshCount <- ThshCountHour/WndwAgrNumc
  rangeSize <- RangeSizeHour/WndwAgrNumc #account for evap in last 24 hours
  if(WndwAgrNumc > 1 | WndwAgrNumc < (5/60)){
    log$fatal('averaging unit must be between 5 minutes and one hour')
    stop()
  }
} else {
  log$fatal('averaging unit needs to be in minutes (min) or hours (hour)')
  stop()
}

# Take stock of our data files.
# !! Try to make more generic, while excluding the manifest.txt file
fileData <- base::list.files(dirInData,pattern='.parquet',full.names=FALSE)
fileFlagsPlau <- base::list.files(dirInFlags,pattern='Plausibility.parquet',full.names=FALSE)
fileFlagsCal <- base::list.files(dirInFlags,pattern='Cal.parquet',full.names=FALSE)

# Check for a manifest file. Ensures full pad.
fileManifest <- base::list.files(dirInData,pattern='manifest',full.names=TRUE)
if(length(fileManifest) == 0){
  stop('No manifest file. Cannot continue.')
} else {
  dayExpc <- base::readLines(fileManifest)
  dayHave <- unlist(lapply(fileData,FUN=function(fileIdx){
    mtch <- regexec(pattern='[0-9]{4}-[0-1]{1}[0-9]{1}-[0-3]{1}[0-9]{1}',fileIdx)[[1]]
    if(mtch != -1){
      dayHaveIdx <- substr(fileIdx,mtch,mtch+attr(mtch,"match.length")-1)
      return(dayHaveIdx)
    } else {
      return(NULL)
    }
  }))
  dayChk <- dayExpc %in% dayHave
  if(!all(dayChk)){
    stop(paste0('Timeseries pad incomplete. Missing the following days: ',
                    paste0(dayExpc[!dayChk],collapse=', ')))

  }
}



# Read the datasets
data <- NEONprocIS.base::def.read.parq.ds(fileIn=fs::path(dirInData,fileData),
                                          VarTime='readout_time',
                                          RmvDupl=TRUE,
                                          Df=TRUE,
                                          log=log)

flagsPlau <- NEONprocIS.base::def.read.parq.ds(fileIn=fs::path(dirInFlags,fileFlagsPlau),
                                          VarTime='readout_time',
                                          RmvDupl=TRUE,
                                          Df=TRUE,
                                          log=log)

flagsCal <- NEONprocIS.base::def.read.parq.ds(fileIn=fs::path(dirInFlags,fileFlagsCal),
                                               VarTime='readout_time',
                                               RmvDupl=TRUE,
                                               Df=TRUE,
                                               log=log)

flags <- dplyr::full_join(flagsPlau, flagsCal, by =  'readout_time')

#combine three gauges into one flagging variable. If any are 1 all flagged, any -1 all flagged, else not flagged
flagNames <- names(flags)[grepl(unique(names(flags)), pattern = 'strainGauge')]
flagNames <- unique(sub(pattern='strainGauge[1-3]Depth',replacement='',x=flagNames))
qfs<- flags[, 'readout_time', drop = F]
for (name in flagNames){
  flags_sub <- flags[,grepl(names(flags), pattern = name)]
  flagVar <- paste0('strainGaugeDepth', name)
  qfs[[flagVar]] <- NA
  flag_0 <- rowSums(flags_sub == 0, na.rm = T)
  qfs[[flagVar]][flag_0 == ncol(flags_sub)] <- 0
  flag_1 <- rowSums(flags_sub == 1, na.rm = T)
  qfs[[flagVar]][flag_1 >=1] <- 1
  flags_neg1 <- rowSums(flags_sub == -1, na.rm = T)
  qfs[[flagVar]][is.na(qfs[[flagVar]]) & flags_neg1 >=1] <- -1
  qfs[[flagVar]][is.na(qfs[[flagVar]])] <- -1
}

# Aggregate depth streams into a single depth.
data$strain_gauge1_stability[is.na(data$strainGauge1Depth)] <- NA
data$strain_gauge2_stability[is.na(data$strainGauge2Depth)] <- NA
data$strain_gauge3_stability[is.na(data$strainGauge3Depth)] <- NA
data <- data %>% dplyr::mutate(strainGaugeDepth = base::rowMeans(x=base::cbind(strainGauge1Depth,
                                                                               strainGauge2Depth,
                                                                               strainGauge3Depth), na.rm = F),
                               strainGaugeStability = base::rowSums(x=base::cbind(strain_gauge1_stability,
                                                                                  strain_gauge2_stability,
                                                                                  strain_gauge3_stability), na.rm = F)==3)
data$strainGaugeDepth[data$strainGaugeStability == FALSE] <- NA

#if there are no heater streams add them in as NA
if(!('internal_temperature' %in% names(data))){data$internal_temperature <- as.numeric(NA)}
if(!('inlet_temperature' %in% names(data))){data$inlet_temperature <- as.numeric(NA)}
if(!('orifice_heater_flag' %in% names(data))){data$orifice_heater_flag <- as.numeric(NA)}

# Add the suspectCal flag to the data so that it can be time-averaged and fed into the final QF
data$strainGaugeDepthSuspectCalQF <- qfs$strainGaugeDepthSuspectCalQF

# Do time averaging
strainGaugeDepthAgr <- data %>%
  dplyr::mutate(startDateTime = lubridate::floor_date(as.POSIXct(readout_time, tz = 'UTC'), unit = WndwAgr)) %>%
  dplyr::mutate(endDateTime = lubridate::ceiling_date(as.POSIXct(readout_time, tz = 'UTC'), unit = WndwAgr,change_on_boundary=TRUE)) %>%
  dplyr::group_by(startDateTime,endDateTime) %>%
  dplyr::summarise(strainGaugeDepth = mean(strainGaugeDepth, na.rm = T),
                   strainGaugeStability = dplyr::if_else(all(is.na(strainGaugeStability)),NA,all(strainGaugeStability==TRUE, na.rm = T)),
                   inletTemperature = mean(inlet_temperature, na.rm = T),
                   internalTemperature = mean(internal_temperature, na.rm = T), 
                   orificeHeaterFlag = max(orifice_heater_flag, na.rm = T), #used to see if heater was on when temps were above heating threshold (heaterErrorQF)
                   inletHeater1QM = round((length(which(orifice_heater_flag == 100))/dplyr::n())*100,1),
                   inletHeater2QM = round((length(which(orifice_heater_flag == 110))/dplyr::n())*100,1), 
                   inletHeater3QM = round((length(which(orifice_heater_flag == 111))/dplyr::n())*100,1),
                   inletHeaterNAQM = round((length(which(is.na(orifice_heater_flag)))/dplyr::n())*100,1),
                   suspectCalQF = max(strainGaugeDepthSuspectCalQF,na.rm = T),
                   .groups = 'drop') 

# Aggregate flags
flagsAgr <- strainGaugeDepthAgr %>% dplyr::select(startDateTime, endDateTime)
flagsAgr$insuffDataQF <- 0
flagsAgr$extremePrecipQF <- 0
flagsAgr$dielNoiseQF <- 0
flagsAgr$strainGaugeStabilityQF <- 0
flagsAgr$strainGaugeStabilityQF[is.na(strainGaugeDepthAgr$strainGaugeStability)] <- -1 
flagsAgr$strainGaugeStabilityQF[strainGaugeDepthAgr$strainGaugeStability == FALSE] <- 1 # Probably make informational flag b/c we removed unstable values
flagsAgr$heaterErrorQF <- 0
flagsAgr$evapDetectedQF <- 0
flagsAgr$heaterErrorQF[strainGaugeDepthAgr$internalTemperature > -6 & 
                         strainGaugeDepthAgr$internalTemperature < 2 & 
                         strainGaugeDepthAgr$inletTemperature < strainGaugeDepthAgr$internalTemperature] <- 1
flagsAgr$heaterErrorQF[strainGaugeDepthAgr$internalTemperature > 6 & strainGaugeDepthAgr$orificeHeaterFlag > 0] <- 1
flagsAgr$heaterErrorQF[is.na(strainGaugeDepthAgr$internalTemperature) |
                         is.na(strainGaugeDepthAgr$inletTemperature) |
                         is.na(strainGaugeDepthAgr$orificeHeaterFlag)] <- -1


# Dynamic Envelope
# Do computation of no-rain days in order to apply a dynamic envelope calculation
# Require that there be no change in depth between the start and end of the day that is
# greater than the pre-defined envelope
EnvelopeConfig <- Envelope # keep the CSV/override value for the summary comparison below
RechargeConfig <- Recharge

dataHourly <- strainGaugeDepthAgr %>%
  dplyr::mutate(startDateTime = lubridate::floor_date(startDateTime, unit = 'hour')) %>%
  dplyr::group_by(startDateTime) %>%
  dplyr::summarise(strainGaugeDepthMed = median(strainGaugeDepth, na.rm = T),
                   strainGaugeDepthMin = min(strainGaugeDepth, na.rm = T), 
                   strainGaugeDepthMax = max(strainGaugeDepth, na.rm = T))
dataDaily <- dataHourly %>%
  dplyr::mutate(startDateTime = lubridate::floor_date(startDateTime, unit = 'day')) %>%
  dplyr::group_by(startDateTime) %>%
  dplyr::summarise(strainGaugeDepthChg = tail(strainGaugeDepthMax,1)-head(strainGaugeDepthMin,1))
setNoRain <- (dataDaily$strainGaugeDepthChg < 0.25*Envelope) & (dataDaily$strainGaugeDepthChg > -1*Envelope)
setNoRain[is.na(setNoRain)] <- FALSE

# Get the envelope if we have determined days without rain
if(any(setNoRain)){
  dayNoRain <- dataDaily$startDateTime[setNoRain]
  timeDay <- lubridate::floor_date(strainGaugeDepthAgr$startDateTime,unit='day')

  envelopeComp <- data.frame(day=unique(timeDay),envelope=as.numeric(NA))

  for (idxDay in unique(timeDay)){
    if(!(idxDay %in% dayNoRain)){
      next
    }
    setDay <- timeDay == idxDay # row indices for this day
    envelopeIdx <- max(strainGaugeDepthAgr$strainGaugeDepth[setDay],na.rm=TRUE)-min(strainGaugeDepthAgr$strainGaugeDepth[setDay],na.rm=TRUE)
    envelopeComp$envelope[envelopeComp$day == idxDay] <- envelopeIdx
  }

  # Take the max envelope
  envelopeMax<- max(envelopeComp$envelope,na.rm=TRUE)
  if(!is.na(envelopeMax)){
    Envelope <- envelopeMax
  }
  
  # if envelope is larger than Recharge threshold adjust recharge.
  if(Envelope > Recharge/3){
    Recharge <- 3*Envelope
  }
}


#initialize fields
strainGaugeDepthAgr$bench <- as.numeric(NA)
strainGaugeDepthAgr$precip <- FALSE #add TRUE when rain detected
strainGaugeDepthAgr$precipType <- as.numeric(NA)

# Initialize & pre-allocate surrogate data for uncertainty analysis
numRow <- nrow(strainGaugeDepthAgr)
nSurr <- 30
surr <- matrix(as.numeric(NA),nrow=numRow,ncol=nSurr)
nameVarDepthS <- paste0('strainGaugeDepthS',seq_len(nSurr))
nameVarBenchS <- paste0('benchS',seq_len(nSurr))
nameVarPrecipS <- paste0('precipS',seq_len(nSurr))
nameVarPrecipTypeS <- paste0('precipTypeS',seq_len(nSurr))
nameVarPrecipBulkS <- paste0('precipBulkS',seq_len(nSurr))
strainGaugeDepthAgr[,nameVarDepthS] <- as.numeric(NA)
strainGaugeDepthAgr[,nameVarBenchS] <- as.numeric(NA)
strainGaugeDepthAgr[,c(nameVarPrecipS)] <- FALSE
strainGaugeDepthAgr[,nameVarPrecipTypeS] <- as.character(NA)
strainGaugeDepthAgr[,nameVarPrecipBulkS] <- as.numeric(NA)

# fixed seed (matches production wrap.precip.aepg.smooth.R) so IAAFT surrogate generation is reproducible and doesn't intermittently fail
set.seed(6366)

for(idxSurr in c(0,seq_len(nSurr))){

  if (idxSurr == 0){
    message(paste0('Running original timeseries'))
    nameVarDepth <- 'strainGaugeDepth'
    nameVarBench <- 'bench'
    nameVarPrecip <- 'precip'
    nameVarPrecipType <- 'precipType'
    nameVarPrecipBulk <- 'precipBulk'
    
    strainGaugeDepthS <- strainGaugeDepthAgr$strainGaugeDepth
    
  } else {
    
    message(paste0('Running Surrogate ',idxSurr))
    nameVarDepth <- paste0('strainGaugeDepthS',idxSurr)
    nameVarBench <- paste0('benchS',idxSurr)
    nameVarPrecip <- paste0('precipS',idxSurr)
    nameVarPrecipType <- paste0('precipTypeS',idxSurr)
    nameVarPrecipBulk <- paste0('precipBulkS',idxSurr)
    
    # If this is the first surrogate, create them
    if(idxSurr == 1){
      
      depthMinusBench <- strainGaugeDepthAgr$strainGaugeDepth - strainGaugeDepthAgr$bench # remove the computed benchmark
      setNotNa <- !is.na(depthMinusBench) # Remove all NA
      
      # Remove rangeSize amount of data at beginning of the timeseries, which can have untracked changes in benchmark that corrupt the surrogates
      setNotNa[1:rangeSize] <- FALSE 
      setNotNa[floor(numRow-rangeSize/2+1):numRow] <- FALSE
      
      # Create surrogates
      surrFill <- base::try(multifractal::iaaft(x=depthMinusBench[setNotNa],N=nSurr),silent=F)
      if("try-error" %in% base::class(surrFill)){
        log$warn('Surrogate generation failed (likely not enough data). Uncertainty estimates cannot be generated. A check will be done later to ensure precip is also NA. If not, an error will occur.')
        break
      }
      strainGaugeDepthAgr[setNotNa,nameVarDepthS] <- strainGaugeDepthAgr$bench[setNotNa] + surrFill    # Add the surrogates to the benchmark

      # Backfill rangeSize amount of data at beginning and end of timeseries to the original strain gauge depth to form a complete timeseries
      strainGaugeDepthAgr[1:rangeSize,nameVarDepthS] <- strainGaugeDepthAgr$strainGaugeDepth[1:rangeSize] 
      strainGaugeDepthAgr[floor(numRow-rangeSize/2+1):numRow,nameVarDepthS] <- strainGaugeDepthAgr$strainGaugeDepth[floor(numRow-rangeSize/2+1):numRow] 
    }
    
    strainGaugeDepthS <- strainGaugeDepthAgr[[nameVarDepth]]
    
  }  
  
  # Run the smoothing algorithm
  precipSmooth <- def.precip.depth.smooth(dateTime=strainGaugeDepthAgr$startDateTime,
                                          gaugeDepth=strainGaugeDepthS,
                                          RangeSize=rangeSize,
                                          Quant=Quant,
                                          ThshCount=ThshCount,
                                          Envelope=Envelope,
                                          ThshChange=ThshChange,
                                          ChangeFactor=ChangeFactor,
                                          ChangeFactorEvap=ChangeFactorEvap,
                                          Recharge=Recharge,
                                          log=log)
  
  
  # Reassign outputs
  strainGaugeDepthAgr[[nameVarBench]] <- precipSmooth$bench
  strainGaugeDepthAgr[[nameVarPrecip]] <- precipSmooth$precip 
  strainGaugeDepthAgr[[nameVarPrecipType]] <- precipSmooth$precipType
  if(idxSurr == 0){
    # Only save evapDetectedQF if we processed the original data
    flagsAgr$evapDetectedQF <- precipSmooth$evapDetectedQF
  }
  
  # Compute precip
  strainGaugeDepthAgr[[nameVarPrecipBulk]] <- c(base::diff(strainGaugeDepthAgr[[nameVarBench]]),as.numeric(NA))
  strainGaugeDepthAgr[[nameVarPrecipBulk]][strainGaugeDepthAgr[[nameVarPrecipBulk]] < 0] <- 0
  strainGaugeDepthAgr[[nameVarPrecipType]] <- c(strainGaugeDepthAgr[[nameVarPrecipType]][2:numRow],as.numeric(NA)) # Shift precip type to align with precip

} # End loop around surrogates


# Compute the uncertainty in precip based on the variability in computed benchmark of the surrogates
# The uncertainty of a sum or difference is equal to their individual uncertainties added in quadrature.
nameVar <- names(strainGaugeDepthAgr)
nameVarBenchS <- c('bench',nameVar[grepl('benchS[0-9]',nameVar)])
strainGaugeDepthAgr$benchS_std <- matrixStats::rowSds(as.matrix(strainGaugeDepthAgr[,nameVarBenchS]))
strainGaugeDepthAgr$precipS_std <- sqrt(strainGaugeDepthAgr$benchS_std^2 + lag(strainGaugeDepthAgr$benchS_std, 1)^2)
strainGaugeDepthAgr$precipS_std <- c(strainGaugeDepthAgr$precipS_std[2:numRow],NA)
strainGaugeDepthAgr$precipS_u95 <- strainGaugeDepthAgr$precipS_std*2


# Compute post-precip quality flags (same frequency as original settings)
flagsAgr$insuffDataQF[is.na(strainGaugeDepthAgr$precipBulk)] <- 1

# Envelope == Massive --> Flag all the data
if(all(flagsAgr$insuffDataQF == 1)){
  flagsAgr$dielNoiseQF <- -1
} else if(Envelope > 10){
  flagsAgr$dielNoiseQF <- 1
}

# Clean up flag logic for NA data
flagsAgr$evapDetectedQF[flagsAgr$insuffDataQF == 1] <- -1
flagsAgr$extremePrecipQF[flagsAgr$insuffDataQF == 1] <- -1

# Join flagsAgr into strainGaugeDepthAgr
strainGaugeDepthAgr <- dplyr::full_join(strainGaugeDepthAgr, flagsAgr, by = c('startDateTime', 'endDateTime'))

# Aggregate to the hour
statsAgrHour <- strainGaugeDepthAgr %>%
    mutate(startDateTime = lubridate::floor_date(startDateTime, '1 hour')) %>%
    mutate(endDateTime = lubridate::ceiling_date(endDateTime, '1 hour')) %>%
    group_by(startDateTime,endDateTime) %>%
    summarise(
              precipBulk = sum(precipBulk),
              insuffDataQF = max(insuffDataQF, na.rm = T),
              extremePrecipQF = max(extremePrecipQF, na.rm = T),
              heaterErrorQF = ifelse(all(is.na(heaterErrorQF)),
                                     -1,
                                     ifelse(sum(heaterErrorQF==1, na.rm=T) >= 0.5*dplyr::n(),
                                            1,
                                            ifelse(all(heaterErrorQF==-1),
                                                   -1,
                                                   0))),
              dielNoiseQF = max(dielNoiseQF, na.rm = T), # Just a placeholder. Computed below.
              strainGaugeStabilityQF = max(strainGaugeStabilityQF, na.rm = T),
              evapDetectedQF = max(evapDetectedQF, na.rm = T),
              inletHeater1QM = mean(inletHeater1QM, na.rm = T),
              inletHeater2QM = mean(inletHeater2QM, na.rm = T),
              inletHeater3QM = mean(inletHeater3QM, na.rm = T),
              inletHeaterNAQM = mean(inletHeaterNAQM, na.rm = T),
              suspectCalQF = max(suspectCalQF,na.rm = T),
              .groups = 'drop'
    )
statsAgrHour$suspectCalQF[!(statsAgrHour$suspectCalQF %in% c(-1,0,1))] <- -1

# Flag for max precip over 60-min - based on hourly totals
statsAgrHour$extremePrecipQF[statsAgrHour$precipBulk > ExtremePrecipMax] <- 1

# Compute hourly final quality flag
statsAgrHour$finalQF <- 0
flags_sub <- statsAgrHour[,c('insuffDataQF','extremePrecipQF', 'heaterErrorQF','suspectCalQF')]
flag_1 <- rowSums(flags_sub == 1, na.rm = T) 
statsAgrHour$finalQF[flag_1 >=1] <- 1 


# Aggregate to the day
statsAgrDay <- statsAgrHour %>%
  mutate(startDate = lubridate::floor_date(startDateTime, '1 day')) %>%
  mutate(endDate = lubridate::ceiling_date(endDateTime, '1 day')) %>%
  group_by(startDate,endDate) %>%
  summarise(precipBulk = sum(precipBulk),
            insuffDataQF = max(insuffDataQF, na.rm = T),
            extremePrecipQF = max(extremePrecipQF, na.rm = T),
            heaterErrorQF = ifelse(all(is.na(heaterErrorQF)),
                                   -1,
                                   ifelse(sum(heaterErrorQF==1, na.rm=T) >= 0.5*dplyr::n(),
                                          1,
                                          ifelse(all(heaterErrorQF==-1),
                                                 -1,
                                                 0))),
            dielNoiseQF = max(dielNoiseQF, na.rm = T),
            strainGaugeStabilityQF = max(strainGaugeStabilityQF, na.rm = T),
            evapDetectedQF = max(evapDetectedQF, na.rm = T),
            inletHeater1QM = mean(inletHeater1QM, na.rm = T),
            inletHeater2QM = mean(inletHeater2QM, na.rm = T),
            inletHeater3QM = mean(inletHeater3QM, na.rm = T),
            inletHeaterNAQM = mean(inletHeaterNAQM, na.rm = T),
            suspectCalQF = max(suspectCalQF,na.rm = T),
            .groups = 'drop'
  )
statsAgrDay$suspectCalQF[!(statsAgrDay$suspectCalQF %in% c(-1,0,1))] <- -1

# Compute daily final quality flag
statsAgrDay$finalQF <- 0
flags_sub <- statsAgrDay[,c('insuffDataQF','extremePrecipQF', 'heaterErrorQF','suspectCalQF')]
flag_1 <- rowSums(flags_sub == 1, na.rm = T) 
statsAgrDay$finalQF[flag_1 >= 1] <- 1 


# Aggregate uncertainty to the hour and day
# Report daily precip, and uncertainty for the central day
# We can use the same equation here, adding the uncertainties for the start and
# end of the day in quadrature, with the caveat that the benchmark does not drop 
# over the course of the hour/day. If this occurs we need to compute for each leg of 
# a flat or increasing benchmark, summing the legs in quadrature

# Hourly
hours <- seq.POSIXt(from=strainGaugeDepthAgr$startDateTime[1],to=strainGaugeDepthAgr$startDateTime[numRow],by='hour')
ucrtAgrHour <- lapply(hours,FUN=function(hourIdx){
  setUcrt <- which(strainGaugeDepthAgr$startDateTime >= hourIdx & 
                     strainGaugeDepthAgr$startDateTime <= (hourIdx + as.difftime(1,units='hours'))) # include first point of next day, because that is the point from which the difference is taken
  ucrtAgr <- def.ucrt.agr.precip.bench(strainGaugeDepthAgr$bench[setUcrt],strainGaugeDepthAgr$benchS_std[setUcrt])
  return(ucrtAgr)
})
statsAgrHour$ucrtExp <- 2*unlist(ucrtAgrHour)
statsAgrHour$ucrtExp[is.na(statsAgrHour$precipBulk)] <- as.numeric(NA) # last value will always be NA because we don't have the start of the next day



# Daily
days <- seq.POSIXt(from=strainGaugeDepthAgr$startDateTime[1],to=strainGaugeDepthAgr$startDateTime[numRow],by='day')
ucrtAgrDay <- lapply(days,FUN=function(dayIdx){
  setUcrt <- which(strainGaugeDepthAgr$startDateTime >= dayIdx & 
                   strainGaugeDepthAgr$startDateTime <= (dayIdx + as.difftime(1,units='days'))) # include first point of next day, because that is the point from which the difference is taken
  ucrtAgr <- def.ucrt.agr.precip.bench(strainGaugeDepthAgr$bench[setUcrt],strainGaugeDepthAgr$benchS_std[setUcrt])
  return(ucrtAgr)
})
statsAgrDay$ucrtExp <- 2*unlist(ucrtAgrDay)
statsAgrDay$ucrtExp[is.na(statsAgrDay$precipBulk)] <- as.numeric(NA) # last value will always be NA because we don't have the start of the next day




# Final check whether uncertainty estimates are provided for all non-NA daily precip output for central 3 days
dayOut <- InfoDirIn$time+as.difftime(c(-1,0,1),units='days')
setOutHour <- statsAgrHour$startDateTime >= min(dayOut) & 
  statsAgrHour$startDateTime < (max(dayOut) + as.difftime(1,units='days'))
setOutDay <- statsAgrDay$startDate >= min(dayOut) & 
  statsAgrDay$startDate < (max(dayOut) + as.difftime(1,units='days'))
# Filter the data for this output day
statsAgrHourIdx <- statsAgrHour[setOutHour,]
statsAgrDayIdx <- statsAgrDay[setOutDay,]

if(!(all.equal(is.na(statsAgrHourIdx$ucrtExp),is.na(statsAgrHourIdx$precipBulk)) == TRUE)){
  stop('There is a mismatch between NA values in hourly precipBulk vs. ucrtExp')
}
if(!(all.equal(is.na(statsAgrDayIdx$ucrtExp),is.na(statsAgrDayIdx$precipBulk)) == TRUE)){
  stop('There is a mismatch between NA values in daily precipBulk vs. ucrtExp')
}

## ---- THRESHOLD CHECK SUMMARY ----------------------------------------------
pctFlag <- function(x) round(100*mean(x==1,na.rm=TRUE),1)
cat('\n==================== THRESHOLD CHECK SUMMARY ====================\n')
cat('Site:',Site,' Date:',Date,' WndwAgr:',WndwAgr,'\n\n')
cat(sprintf('%-20s %12s %12s\n','Parameter','Config','Used/Derived'))
cat(sprintf('%-20s %12s %12s\n','Envelope',round(EnvelopeConfig,3),round(Envelope,3)))
cat(sprintf('%-20s %12s %12s\n','Recharge',round(RechargeConfig,3),round(Recharge,3)))
cat(sprintf('%-20s %12s\n','ThshCountHour',ThshCountHour))
cat(sprintf('%-20s %12s\n','ChangeFactor',ChangeFactor))
cat(sprintf('%-20s %12s\n','ChangeFactorEvap',ChangeFactorEvap))
cat(sprintf('%-20s %12s\n','ThshChange',ThshChange))
cat(sprintf('%-20s %12s\n','ExtremePrecipMax',ExtremePrecipMax))
cat(sprintf('%-20s %12s\n','Quant',Quant))
cat(sprintf('%-20s %12s\n','RangeSizeHour',RangeSizeHour))
cat('\n---- % of hourly intervals flagged ----\n')
cat('insuffDataQF:   ',pctFlag(statsAgrHour$insuffDataQF),'%\n')
cat('extremePrecipQF:',pctFlag(statsAgrHour$extremePrecipQF),'%\n')
cat('heaterErrorQF:  ',pctFlag(statsAgrHour$heaterErrorQF),'%\n')
cat('suspectCalQF:   ',pctFlag(statsAgrHour$suspectCalQF),'%\n')
cat('dielNoiseQF:    ',pctFlag(statsAgrHour$dielNoiseQF),'%\n')
cat('evapDetectedQF: ',pctFlag(statsAgrHour$evapDetectedQF),'%\n')
cat('finalQF:        ',pctFlag(statsAgrHour$finalQF),'%\n')
cat('\nTotal precipBulk over period (mm):',round(sum(statsAgrDay$precipBulk,na.rm=TRUE),2),'\n')
cat('===================================================================\n\n')

## ---- PLOTS: eyeball whether the current thresholds look reasonable -------
titleTag <- paste0(Site,' ',Date,' (',WndwAgr,')')

# Clean view: aggregated gauge depth vs. the computed benchmark, with detected precip points marked -
# the main plot for checking whether Envelope/ThshCountHour/etc. are still triggering precip sensibly.
precipPts <- strainGaugeDepthAgr[!is.na(strainGaugeDepthAgr$precip) & strainGaugeDepthAgr$precip == TRUE,]
dfClean <- strainGaugeDepthAgr[,c('startDateTime','strainGaugeDepth','bench')] %>%
  tidyr::pivot_longer(cols=-startDateTime,names_to='variable',values_to='value')
print(plotly::plot_ly() %>%
  plotly::add_lines(data=dfClean,x=~startDateTime,y=~value,color=~variable) %>%
  plotly::add_markers(data=precipPts,x=~startDateTime,y=~bench,name='precip detected',marker=list(color='red',size=6)) %>%
  plotly::layout(title=paste0(titleTag,' - Depth vs. Benchmark')))

# Full diagnostic view including the uncertainty surrogates (set ShowSurrogates <- FALSE above to skip)
if(ShowSurrogates){
df <- strainGaugeDepthAgr[,c('startDateTime','strainGaugeDepth',paste0('strainGaugeDepthS',seq_len(nSurr)),'bench',paste0('benchS',seq_len(nSurr)))] %>%
  tidyr::pivot_longer(cols=-startDateTime,names_to='variable',values_to='value')
plotly::plot_ly(data=df,x=~startDateTime,y=~value,color=~variable,mode='lines') %>%
  plotly::layout(title=paste0(titleTag,' - With Surrogates'))
}

df <- strainGaugeDepthAgr[,c('startDateTime','precipBulk','precipS_u95')] %>%
  tidyr::pivot_longer(cols=-startDateTime,names_to='variable',values_to='value')
plotly::plot_ly(data=df,x=~startDateTime,y=~value,color=~variable,mode='lines') %>%
  plotly::layout(title=titleTag)

df <- statsAgrHour[,c('startDateTime','precipBulk','ucrtExp')] %>%
  tidyr::pivot_longer(cols=-startDateTime,names_to='variable',values_to='value')
plotly::plot_ly(data=df,x=~startDateTime,y=~value,color=~variable,mode='lines') %>%
  plotly::layout(title=paste0(titleTag,' - Hourly'))

df <- statsAgrDay[,c('startDate','precipBulk','ucrtExp')] %>%
  tidyr::pivot_longer(cols=-startDate,names_to='variable',values_to='value')
plotly::plot_ly(data=df,x=~startDate,y=~value,color=~variable,mode='lines')  %>%
  plotly::layout(title=paste0(titleTag,' - Daily'))
