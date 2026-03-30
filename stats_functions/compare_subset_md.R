# Compare Mahalanobis distance distributions between a subset and all sites.
# - df: data with numeric meteorology columns and a logical subset flag column
# - vars: meteorology variables to use
# - subset_flag: column name (logical) indicating sites eligible for intercomparison
# Returns a list with distances, summaries, and KS test.
compare_subset_md <- function(df, vars, subset_flag,
                              scale = TRUE, na_action = "impute_median",
                              ridge = 0.05, squared = FALSE) {
  stopifnot(is.logical(df[[subset_flag]]))
  md <- compute_mahalanobis(df, vars,
                            scale = scale,
                            na_action = na_action,
                            ridge = ridge,
                            squared = squared)
  sub_idx <- which(df[[subset_flag]] & !is.na(md))
  all_idx <- which(!is.na(md))
  
  md_subset <- md[sub_idx]
  md_all    <- md[all_idx]
  
  ks <- stats::ks.test(md_subset, md_all)
  
  list(
    md = md,
    md_subset = md_subset,
    md_all = md_all,
    summary_subset = summary(md_subset),
    summary_all = summary(md_all),
    ks_test = ks
  )
}
