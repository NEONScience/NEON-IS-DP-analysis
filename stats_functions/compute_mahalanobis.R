
# Compute Mahalanobis distances with options for scaling, NA handling, and covariance regularization.
# - df: data.frame or tibble with numeric variables
# - vars: character vector of column names to use
# - center: optional numeric vector of length p (defaults to column means of df[vars])
# - covmat: optional covariance matrix (defaults to cov of df[vars])
# - scale: if TRUE, standardize columns to mean 0, sd 1 before computing covmat/center
# - na_action: "omit" (drop rows with any NA in vars) or "impute_median" (median per column)
# - ridge: non-negative scalar added to the diagonal of covariance (stabilizes inversion if p ~ n or collinearity)
# - squared: if TRUE, return squared MD; else Euclidean MD in the covariance metric
compute_mahalanobis <- function(df, vars,
                                center = NULL,
                                covmat = NULL,
                                scale = TRUE,
                                na_action = c("omit", "impute_median"),
                                ridge = 0,
                                squared = FALSE) {
  stopifnot(is.data.frame(df))
  na_action <- match.arg(na_action)
  
  X <- df[, vars, drop = FALSE]
  
  # ensure numeric
  if (!all(vapply(X, is.numeric, logical(1)))) {
    stop("All 'vars' must be numeric.")
  }
  
  # NA handling
  if (na_action == "omit") {
    keep <- stats::complete.cases(X)
    X <- X[keep, , drop = FALSE]
    dropped <- which(!keep)
  } else {
    # Median imputation (robust for skewed met vars)
    for (j in seq_along(vars)) {
      v <- X[[j]]
      if (anyNA(v)) {
        med <- stats::median(v, na.rm = TRUE)
        v[is.na(v)] <- med
        X[[j]] <- v
      }
    }
    dropped <- integer(0)
  }
  
  # Optionally standardize columns (recommended when variables are on different scales)
  if (scale) {
    X_scaled <- scale(X)                    # attributes retain center and scale
    scale_center <- attr(X_scaled, "scaled:center")
    scale_scale  <- attr(X_scaled, "scaled:scale")
  } else {
    X_scaled <- as.matrix(X)
    scale_center <- rep(0, ncol(X_scaled))
    scale_scale  <- rep(1, ncol(X_scaled))
  }
  
  # Determine center in the working (possibly scaled) space
  if (is.null(center)) {
    center_work <- colMeans(X_scaled)
  } else {
    center <- as.numeric(center)
    if (length(center) != ncol(X_scaled)) {
      stop("'center' must have length equal to number of variables.")
    }
    # If user supplied center in the original scale AND we scaled features,
    # map it to the scaled space.
    if (scale) {
      center_work <- (center - scale_center) / scale_scale
    } else {
      center_work <- center
    }
  }
  
  # Determine covariance matrix in the working space
  if (is.null(covmat)) {
    # Use unbiased covariance (n-1 denominator)
    S <- stats::cov(X_scaled)
  } else {
    S <- as.matrix(covmat)
    # If user supplied cov in original scale and we scaled features,
    # map it to the scaled space: S_scaled = D^{-1} S D^{-1}, where D = diag(sd)
    if (scale) {
      Dinv <- diag(1 / scale_scale, nrow = length(scale_scale))
      S <- Dinv %*% S %*% Dinv
    }
  }
  
  # Ridge regularization for numerical stability
  if (ridge < 0) stop("'ridge' must be non-negative.")
  if (ridge > 0) {
    S <- S + diag(ridge, nrow = ncol(S))
  }
  
  # Invert covariance; use solve with safeguard
  S_inv <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(S_inv)) {
    # Fall back to a pseudo-inverse if needed (no extra packages)
    S_svd <- svd(S)
    tol <- max(dim(S)) * max(S_svd$d) * .Machine$double.eps
    d_inv <- ifelse(S_svd$d > tol, 1 / S_svd$d, 0)
    S_inv <- S_svd$u %*% (diag(d_inv, length(d_inv)) %*% t(S_svd$v))
  }
  
  # Compute squared MD: (x - mu)^T S^{-1} (x - mu)
  diffs <- sweep(X_scaled, 2, center_work, FUN = "-")
  md2 <- rowSums((diffs %*% S_inv) * diffs)
  
  # Prepare output aligned to input rows
  out <- rep(NA_real_, nrow(df))
  if (na_action == "omit" && length(dropped)) {
    idx <- which(stats::complete.cases(df[, vars, drop = FALSE]))
    out[idx] <- if (squared) md2 else sqrt(md2)
  } else {
    out[] <- if (squared) md2 else sqrt(md2)
  }
  
  # Attach attributes for transparency/repro
  attr(out, "center_used_scaled_space") <- center_work
  attr(out, "cov_used_scaled_space") <- S
  attr(out, "scaled") <- scale
  attr(out, "ridge") <- ridge
  attr(out, "dropped_rows_due_to_NA") <- if (length(dropped)) dropped else integer(0)
  out
}

