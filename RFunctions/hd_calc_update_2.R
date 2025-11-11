hd_calc_update_2 <- function (data, reference, biomarkers) 
{
  ref1 = as.matrix(reference[, c(biomarkers, "weights")])
  dat = as.matrix(data[, biomarkers])
  dat = na.omit(dat)
  ref1 = na.omit(ref1)
  ref = ref1[,biomarkers]
  
  for (j in 1:ncol(dat)) {
    dat[, j] <- (dat[, j] - weighted.mean(ref[, j], ref1[,"weights"], na.rm = TRUE))/sqrt(diag(cov.wt(ref, wt = ref1[,"weights"])$cov)[j])
  }
  for (j in 1:ncol(ref)) {
    ref[, j] <- (ref[, j] - weighted.mean(ref[, j], ref1[,"weights"], na.rm = TRUE))/sqrt(diag(cov.wt(ref, wt = ref1[,"weights"])$cov)[j])
  }
  if (nrow(ref) == 1) {
    warning("The reference matrix must have more than one row")
  }
  else {
    means = apply(ref, 2, weighted.mean, w = ref1[,"weights"], na.rm = TRUE)
    cv_mat = cov.wt(ref, wt = ref1[,"weights"])$cov
  }
  if (nrow(dat) == 1) {
    warning("The function does not work with single-row data")
  }
  else {
    dat = as.matrix(dat)
    hd = rep(NA, nrow(dat))
    for (x in 1:nrow(dat)) {
      hd[x] <- sqrt((dat[x, ] - means) %*% solve(cv_mat) %*% 
                      (dat[x, ] - means))
    }
  }
  dat = data %>% select(sampleID, all_of(biomarkers)) %>% 
    na.omit()
  dat$hd = hd/sd(hd)
  dat$hd_log = log(hd)/sd(log(hd))
  nobs = sum(!is.na(dat$hd))
  dat = left_join(data, dat[, c("sampleID", "hd", "hd_log")], 
                  by = "sampleID")
  fit = list(mcov = means, cov_mat = cv_mat, nobs = nobs)
  hd = list(data = as.data.frame(dat), fit = fit)
  class(hd) = append(class(hd), "hd")
  return(hd)
}
