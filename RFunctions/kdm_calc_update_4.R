# Without including Sba2 in calculation
# Calculate weighted variance of chronological age

#build a formula
form <- function(y,x){
  return(as.formula(paste0(y,'~',x)))
}

#get effects from linear models on biomarkers
get_effs <- function(mod){
  m = summary(mod)
  res = data.frame(q=coef(mod)['(Intercept)'],k=coef(mod)[2],s=NA,r=NA)
  if('glm' %in% attr(mod,'class')){
    #sd of residuals is the RMSE
    res$s = sd(mod$residuals)
    #calculate s=pseudo r-squared
    res$r = 1-(mod$deviance/mod$null.deviance)
  }else{
    res$s = m$sigma
    res$r = m$r.squared
  }
  
  return(res)
}

kdm_calc_update_4 <- function (data, biomarkers, fit = NULL, s_ba2 = NULL) 
{
  dat = data
  rm(data)
  bm = biomarkers
  rm(biomarkers)
  design = survey::svydesign(id = ~1, weights = ~weights, data = dat)
  if (is.null(fit)) {
    lm_age = lapply(bm, function(marker) {
      survey::svyglm(form(marker, "age"), design = design, 
                     family = gaussian())
    })
    agev = do.call(rbind, lapply(lm_age, get_effs)) %>% 
      mutate(bm = bm, r1 = abs((k/s) * sqrt(r)), r2 = abs(k/s), 
             n2 = (k/s)^2)
    rm(lm_age)
    age_range = range(dat$age, na.rm = TRUE)
    age_var <- Hmisc::wtd.var(dat$age, weights = dat$weights)
    rchar = sum(agev$r1)/sum(agev$r2)
    s_r = ((1 - (rchar^2))/(rchar^2)) * (age_var/nrow(agev))
  }
  else {
    agev = fit$lm_age
    s_r = fit$s_r
  }
  n1 = dat[, bm]
  for (m in colnames(n1)) {
    row = which(agev$bm == m)
    obs = (dat[, m] - agev[row, "q"]) * (agev[row, "k"]/(agev[row, 
                                                              "s"]^2))
    n1[, m] = (obs)
  }
  BA_nmiss = apply(n1, 1, function(x) sum(is.na(x)))
  BA_obs = length(bm) - BA_nmiss
  BAe_n = rowSums(n1, na.rm = TRUE)
  BAe_d = sum(agev$n2, na.rm = TRUE)
  dat = dat %>% mutate(BA_eo = BAe_n/BAe_d, BA_e = (BA_eo/(BA_obs)) * 
                         length(bm))
  dat$BA_CA = unlist(dat$BA_e - dat$age)
  t1 = (dat$BA_CA - mean(dat$BA_CA, na.rm = TRUE))^2
  s2 = mean(t1, na.rm = TRUE)
  nobs = sum(!is.na(dat$BA_CA))
  if (is.null(s_ba2)) {
    s_ba2 = s2 - s_r
  }
  else {
    s_ba2 = s_ba2
  }
  dat$kdm = unlist((BAe_n)/(BAe_d))
  dat$kdm = ifelse(BA_nmiss > 2, NA, dat$kdm)
  dat$kdm_advance = dat$kdm - dat$age
  dat$BA_eo = NULL
  dat$BA_e = NULL
  dat$BA_CA = NULL
  fit = list(lm_age = agev, s_r = s_r, s_ba2 = s_ba2, s2 = s2, 
             nobs = nobs)
  kdm = list(data = as.data.frame(dat), fit = fit)
  class(kdm) = append(class(kdm), "kdm")
  return(kdm)
}
