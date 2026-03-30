# Compute and compare Mahalanobis distance for intercomparison sites vs. 
# all sites to evaluate representativeness of intercomparison.

df=data.frame(site_id=c('BLUE','BONA','CPER','GUAN','HARV','OSBS','PRIN','REDB','SRER','TOOL','WOOD','WREF','YELL','ARIK','CLBJ','KONZ','NIWO','ONAQ','ORNL','PUUM','SCBI','SJER','TALL','UNDE'),
              eligible=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
              MAT=c(16,-1,8,25,8,20,18,10,20,-4,5,8,0,11,18,12,0,9,15,13,13,17,17,3),
              MAP=c(980,399,370,1168,967,1290,840,713,290,331,490,2530,509,449,840,860,758,388,1222,2685,1054,270,1350,854),
              envelope=c(3.15,3.1,5.6,0.9,2.2,2.1,1.3,6,1.1,3,1.3,2.1,11.8,1.4,2.3,2.7,4.3,1.3,2.4,1.1,5.9,2.1,2.1,1.5)
            )
vars <- c("MAT", "MAP","envelope")


# Compute distances for all sites (center/cov from all sites by default)
md <- compute_mahalanobis(df, vars, scale = TRUE, ridge = 0.05)

# Quick check: compare subset vs. all distribution
cmp <- compare_subset_md(df, vars, subset_flag = "eligible",
                         scale = TRUE, ridge = 0.05)

cmp$summary_subset
cmp$summary_all
cmp$ks_test

# Optional: simple visualization
op <- par(mfrow = c(1, 2))
hist(cmp$md_all, breaks = 20, col = "grey85", border = "white",
     main = "Mahalanobis Distance (All Sites)", xlab = "MD")
hist(cmp$md_subset, breaks = 15, col = "#3182bd", border = "white",
     main = "Mahalanobis Distance (Eligible Subset)", xlab = "MD")
par(op)
