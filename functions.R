library(dplyr)
library(purrr)
library(tibble)
library(lme4)
library(clubSandwich)
library(splines)

# Tables and Primary Inference 
calc_cif_at_time <- function(data,time,status,eval_time = 12,failcode,cencode = 0,cohort_level,cohort_var,digits = 1,ci_method = c("loglog", "wald")) {
  ci_method <- match.arg(ci_method)
  
  library(dplyr)
  library(cmprsk)
  library(gtsummary)
  
  d <- data
  
  if (!is.null(cohort_level) && !is.null(cohort_var)) {
    d <- d |> filter(.data[[cohort_var]] == cohort_level)
  }
  
  ci <- cmprsk::cuminc(
    ftime = d[[time]],
    fstatus = d[[status]],
    cencode = cencode
  )
  
  # Locate the relapse curve
  nm <- names(ci)[names(ci) == as.character(failcode)]
  
  if (length(nm) == 0) {
    nm <- names(ci)[grepl(paste0(" ", failcode, "$|^", failcode, "$"), names(ci))]
  }
  
  x <- ci[[nm[1]]]
  
  # CIF is a step function: use estimate at the latest time <= eval_time
  idx <- max(which(x$time <= eval_time))
  
  if (!is.finite(idx)) {
    return(NA_character_)
  }
  
  est <- x$est[idx]
  se <- sqrt(x$var[idx])
  z <- qnorm(0.975)
  
  if (ci_method == "wald") {
    lcl <- max(0, est - z * se)
    ucl <- min(1, est + z * se)
  }
  
  if (ci_method == "loglog") {
    # log(-log(1 - CIF)) transformation
    if (est <= 0 || est >= 1 || is.na(se)) {
      lcl <- max(0, est - z * se)
      ucl <- min(1, est + z * se)
    } else {
      theta <- log(-log(1 - est))
      se_theta <- se / ((1 - est) * abs(log(1 - est)))
      
      lcl <- 1 - exp(-exp(theta - z * se_theta))
      ucl <- 1 - exp(-exp(theta + z * se_theta))
    }
  }
  
  paste0(
    gtsummary::style_percent(est, digits = digits),
    " (95% CI ",
    gtsummary::style_percent(lcl, digits = digits),
    "–",
    gtsummary::style_percent(ucl, digits = digits),
    ")"
  )
}
fmt_binary <- function(data, outcome, cohort_level, cohort_var) {
  d <- data |> filter(.data[[cohort_var]] == cohort_level)
  n_event <- sum(d[[outcome]] == 1, na.rm = TRUE)
  n_total <- sum(!is.na(d[[outcome]]))
  pct <- 100 * n_event / n_total
  
  paste0(n_event, "/", n_total, " (", style_number(pct, digits = 1), "%)")
}
fmt_continuous <- function(data, outcome, cohort_level, cohort_var) {
  d <- data |> filter(.data[[cohort_var]] == cohort_level)
  sts = quantile(d[[outcome]], probs=c(0.25, 0.5, 0.75), na.rm = T)
  paste0(sts[2], " (", sts[1], '-', sts[3], ")")
}
fmt_median_surv <- function(data, time, event, cohort_level, cohort_var) {
  d <- data |> filter(.data[[cohort_var]] == cohort_level)
  
  fit <- survfit(
    as.formula(paste0("Surv(", time, ", ", event, ") ~ 1")),
    data = d
  )
  
  s <- summary(fit)$table
  
  med <- unname(s["median"])
  lcl <- unname(s["0.95LCL"])
  ucl <- unname(s["0.95UCL"])
  
  fmt <- function(x) {
    ifelse(is.na(x), "NR", style_number(x, digits = 1))
  }
  
  paste0(fmt(med), " (", fmt(lcl), "–", fmt(ucl), ")")
}
fmt_wilcox <- function(data, outcome, cohort_levels, cohort_var) {
  twilc <- wilcox.test(
    x= (data %>% filter(cohort_2lvl==cohort_levels[2]) %>% arrange(subclass))[[outcome]],
    y= (data %>% filter(cohort_2lvl!=cohort_levels[2])%>% arrange(subclass))[[outcome]],
    paired=T,
    alternative=c("two.sided")
  )
  paste0('p=', style_number(twilc$p.value, digits=2))
}
fmt_or_unmatch <- function(data, outcome, cohort_var) {
  f <- as.formula(paste0('`', outcome, '`', " ~ ", cohort_var ))
  
  fit <- glm(f, data = data, family='binomial')
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term != "(Intercept)") |>
    transmute(
      effect = paste0(
        "OR ",
        style_number(estimate, digits = 2),
        " (",
        style_number(conf.low, digits = 2),
        "–",
        style_number(conf.high, digits = 2),
        ")"
      )
    ) |>
    pull(effect)
}
fmt_or_pval_unmatch <- function(data, outcome, cohort_var) {
  f <- as.formula(paste0('`', outcome, '`', " ~ ", cohort_var ))
  
  fit <- glm(f, data = data, family='binomial')
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term != "(Intercept)") |>
    transmute(
      pval = paste0(
        "p=",
        style_number(p.value, digits = 2)
      )
    ) |>
    pull(pval)
}
fmt_or <- function(data, outcome, cohort_var) {
  f <- as.formula(paste0('`', outcome, '`', " ~ ", cohort_var, "+ strata(subclass)" ))
  
  fit <- clogit(f, data = data)
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term != "(Intercept)") |>
    transmute(
      effect = paste0(
        "OR ",
        style_number(estimate, digits = 2),
        " (",
        style_number(conf.low, digits = 2),
        "–",
        style_number(conf.high, digits = 2),
        ")"
      )
    ) |>
    pull(effect)
}
fmt_or_pval <- function(data, outcome, cohort_var) {
  f <- as.formula(paste0('`', outcome, '`', " ~ ", cohort_var, "+ strata(subclass)" ))
  
  fit <- clogit(f, data = data)
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term != "(Intercept)") |>
    transmute(
      pval = paste0(
        "p=",
        style_number(p.value, digits = 2)
      )
    ) |>
    pull(pval)
}
fmt_hr <- function(data, time, event, cohort_var) {
  f <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ", cohort_var, "+ strata(subclass)"))
  
  fit <- coxph(f, data = data, cluster=subclass)
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    transmute(
      effect = paste0(
        "HR ",
        style_number(estimate, digits = 2),
        " (",
        style_number(conf.low, digits = 2),
        "–",
        style_number(conf.high, digits = 2),
        ")"
      )
    ) |>
    pull(effect)
}
fmt_hr_pval <- function(data, time, event, cohort_var) {
  f <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ", cohort_var, "+ strata(subclass)"))
  
  fit <- coxph(f, data = data, cluster=subclass)
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    transmute(
      pval = paste0(
        "p=",
        style_number(p.value, digits = 2)
      )
    ) |>
    pull(pval)
}
fmt_cihr <- function(data, time, event, cohort_var, etype) {
  f <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ."))
  fgp <- finegray(f, etype=etype, data=data %>% select(tt_pfs_m, ev_rel, subclass, cohort_2lvl))
  
  f1 <- as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", cohort_var, "+ strata(subclass)"))
  fit <- coxph(f1, weight=fgwt, data = fgp, cluster=subclass)
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    transmute(
      effect = paste0(
        "HR ",
        style_number(estimate, digits = 2),
        " (",
        style_number(conf.low, digits = 2),
        "–",
        style_number(conf.high, digits = 2),
        ")"
      )
    ) |>
    pull(effect)
}
fmt_cihr_pval <- function(data, time, event, cohort_var, etype) {
  f <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ."))
  fgp <- finegray(f, etype=etype, data=data %>% select(tt_pfs_m, ev_rel, subclass, cohort_2lvl))
  
  f1 <- as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", cohort_var, "+ strata(subclass)"))
  fit <- coxph(f1, weight=fgwt, data = fgp, cluster=subclass)
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    transmute(
      pval = paste0(
        "p=",
        style_number(p.value, digits = 2)
      )
    ) |>
    pull(pval)
}


# Wald test across strata for objects with complex interactions and robust estimation of paired/clustered variance. 
wald_contrast <- function(fit, L, label) {
  b <- coef(fit)
  V <- vcov(fit)
  
  L_full <- setNames(rep(0, length(b)), names(b))
  L_full[names(L)] <- L
  
  log_hr <- sum(L_full * b)
  se <- sqrt(as.numeric(t(L_full) %*% V %*% L_full))
  z <- log_hr / se
  p <- 2 * pnorm(-abs(z))
  
  data.frame(
    contrast = label,
    log_HR = log_hr,
    HR = exp(log_hr),
    lower_95 = exp(log_hr - 1.96 * se),
    upper_95 = exp(log_hr + 1.96 * se),
    z = z,
    p = p
  )
}
hr_curve_from_spline <- function(fit, logcrp_values) {
  b <- coef(fit)
  V <- vcov(fit)
  
  new0 <- data.frame(
    trt = 0,
    logcrp = logcrp_values
  )
  
  new1 <- data.frame(
    trt = 1,
    logcrp = logcrp_values
  )
  
  X0 <- model.matrix(delete.response(terms(fit)), new0)
  X1 <- model.matrix(delete.response(terms(fit)), new1)
  
  X0 <- X0[, names(b), drop = FALSE]
  X1 <- X1[, names(b), drop = FALSE]
  
  Xdiff <- X1 - X0
  
  log_hr <- as.numeric(Xdiff %*% b)
  se <- sqrt(rowSums((Xdiff %*% V) * Xdiff))
  
  data.frame(
    logcrp = logcrp_values,
    HR = exp(log_hr),
    lower_95 = exp(log_hr - 1.96 * se),
    upper_95 = exp(log_hr + 1.96 * se)
  )
}

library(ggsci)
colors <- c(
  nejm_red="#BC3C29FF",
  nejm_blue="#0072B5FF",
  nejm_orange="#E18727FF",
  nejm_green="#20854EFF",
  nejm_purple="#7876B1FF",
  
  blue="#1F5AA6",
  green="#087F8C"
)

## Linear mixed effects modeling functions
# Helper functions
.get_vcov <- function(fit, dat, cluster_var = NULL) {
  
  beta <- lme4::fixef(fit)
  
  V <- if (!is.null(cluster_var)) {
    tryCatch(
      as.matrix(
        clubSandwich::vcovCR(
          fit,
          cluster = dat[[cluster_var]],
          type = "CR2"
        )
      ),
      error = function(e) as.matrix(vcov(fit))
    )
  } else {
    as.matrix(vcov(fit))
  }
  
  V[names(beta), names(beta), drop = FALSE]
}
.add_reference_covars <- function(newdat, fit, covariates) {
  
  mf <- model.frame(fit)
  
  for (v in covariates) {
    
    if (!v %in% names(mf)) next
    
    if (is.factor(mf[[v]])) {
      newdat[[v]] <- levels(mf[[v]])[1]
    } else if (is.numeric(mf[[v]])) {
      newdat[[v]] <- ifelse(
        grepl("_c$", v),
        0,
        median(mf[[v]], na.rm = TRUE)
      )
    } else {
      newdat[[v]] <- mf[[v]][which(!is.na(mf[[v]]))[1]]
    }
  }
  
  for (v in intersect(names(newdat), names(mf))) {
    if (is.factor(mf[[v]])) {
      newdat[[v]] <- factor(
        as.character(newdat[[v]]),
        levels = levels(mf[[v]])
      )
    }
  }
  
  newdat
}
.get_fixed_X <- function(fit, newdat) {
  
  mf <- model.frame(fit)
  beta <- lme4::fixef(fit)
  
  fixed_terms <- delete.response(terms(fit, fixed.only = TRUE))
  
  X <- model.matrix(
    fixed_terms,
    data = newdat,
    contrasts.arg = attr(lme4::getME(fit, "X"), "contrasts"),
    xlev = stats::.getXlevels(fixed_terms, mf)
  )
  
  if (!all(names(beta) %in% colnames(X))) {
    stop("Prediction matrix is missing model coefficient columns.")
  }
  
  X[, names(beta), drop = FALSE]
}
.cell_w <- function(newdat, group_var, group_value, time_var, time_value) {
  
  if (is.numeric(newdat[[time_var]])) {
    idx <- as.character(newdat[[group_var]]) == group_value &
      abs(newdat[[time_var]] - time_value) < 1e-8
  } else {
    idx <- as.character(newdat[[group_var]]) == group_value &
      as.character(newdat[[time_var]]) == as.character(time_value)
  }
  
  if (sum(idx) != 1) {
    stop(
      "Could not identify exactly one cell for ",
      group_var, " = ", group_value,
      ", ", time_var, " = ", time_value
    )
  }
  
  as.numeric(idx)
}
.calc_contrast <- function(w,X,beta,V,contrast,contrast_class,visit_f = NA_character_,day_rel_car = NA_real_,baseline_visit = NA_character_,baseline_day = NA_real_) {
  
  cvec <- as.numeric(w %*% X)
  est <- sum(cvec * beta)
  var_est <- as.numeric(cvec %*% V %*% cvec)
  se <- sqrt(pmax(var_est, 0))
  z <- est / se
  p <- 2 * pnorm(-abs(z))
  
  tibble(
    visit_f = visit_f,
    day_rel_car = day_rel_car,
    baseline_visit = baseline_visit,
    baseline_day = baseline_day,
    contrast_class = contrast_class,
    contrast = contrast,
    estimate_log2 = est,
    se = se,
    z = z,
    p = p,
    fold_change = 2^est,
    lower_95_log2 = est - 1.96 * se,
    upper_95_log2 = est + 1.96 * se,
    lower_95_fc = 2^(est - 1.96 * se),
    upper_95_fc = 2^(est + 1.96 * se)
  )
}
.add_fit_meta <- function(out, dat, fit, cluster_var = NULL) {
  
  out %>%
    mutate(
      n_obs = nrow(dat),
      n_patients = n_distinct(dat$record_id),
      n_clusters = if (!is.null(cluster_var)) n_distinct(dat[[cluster_var]]) else NA_integer_,
      singular_fit = lme4::isSingular(fit, tol = 1e-4)
    )
}

# LMM Calls
fit_cat_trt_lmm <- function(dat_lab,covars,cluster_var='subclass',fml) {

  fit <- lme4::lmer(
    fml,
    data = dat_lab,
    REML = FALSE
  )
  
  beta <- lme4::fixef(fit)
  V <- .get_vcov(fit, dat_lab, cluster_var)
  
  mf <- model.frame(fit)
  
  newdat <- expand.grid(
    trt_f = levels(mf$trt_f),
    visit_f = levels(mf$visit_f),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  newdat <- .add_reference_covars(newdat, fit, covars)
  X <- .get_fixed_X(fit, newdat)
  
  w_cell <- function(trt_value, visit_value) {
    .cell_w(newdat, "trt_f", trt_value, "visit_f", visit_value)
  }
  
  visits <- levels(mf$visit_f)
  baseline_visit <- if ("Baseline" %in% visits) "Baseline" else NA_character_
  
  out <- list()
  
  for (v in visits) {
    
    w_trt <- w_cell("anakinra", v) - w_cell("control", v)
    
    out[[length(out) + 1]] <- .calc_contrast(
      w_trt, X, beta, V,
      contrast = paste0("Anakinra vs control at ", v),
      contrast_class = "treatment_effect",
      visit_f = v,
      baseline_visit = baseline_visit
    )
    
    if (!is.na(baseline_visit) && v != baseline_visit) {
      
      w_base <- w_cell("anakinra", baseline_visit) -
        w_cell("control", baseline_visit)
      
      out[[length(out) + 1]] <- .calc_contrast(
        w_trt - w_base, X, beta, V,
        contrast = paste0("DID: ", v, " vs ", baseline_visit),
        contrast_class = "did_vs_baseline",
        visit_f = v,
        baseline_visit = baseline_visit
      )
    }
  }
  
  bind_rows(out) %>%
    .add_fit_meta(dat_lab, fit, cluster_var)
}
fit_cat_schedule_lmm <- function(dat_lab,covars,cluster_var = "subclass", fml) {
  
  fit <- lme4::lmer(
    fml,
    data = dat_lab,
    REML = FALSE
  )
  
  beta <- lme4::fixef(fit)
  V <- .get_vcov(fit, dat_lab, cluster_var)
  
  mf <- model.frame(fit)
  
  newdat <- expand.grid(
    akra_cohort_f = levels(mf$akra_cohort_f),
    visit_f = levels(mf$visit_f),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  newdat <- .add_reference_covars(newdat, fit, covars)
  X <- .get_fixed_X(fit, newdat)
  
  w_cell <- function(group_value, visit_value) {
    .cell_w(newdat, "akra_cohort_f", group_value, "visit_f", visit_value)
  }
  
  visits <- levels(mf$visit_f)
  baseline_visit <- if ("Baseline" %in% visits) "Baseline" else NA_character_
  
  te_c1 <- function(v) {
    w_cell("cohort1", v) - w_cell("cohort1 control", v)
  }
  
  te_c2 <- function(v) {
    w_cell("cohort2", v) - w_cell("cohort2 control", v)
  }
  
  out <- list()
  
  for (v in visits) {
    
    w_c1 <- te_c1(v)
    w_c2 <- te_c2(v)
    
    out[[length(out) + 1]] <- .calc_contrast(
      w_c1, X, beta, V,
      contrast = paste0("Cohort1 anakinra vs cohort1 control at ", v),
      contrast_class = "cohort1_treatment_effect",
      visit_f = v,
      baseline_visit = baseline_visit
    )
    
    out[[length(out) + 1]] <- .calc_contrast(
      w_c2, X, beta, V,
      contrast = paste0("Cohort2 anakinra vs cohort2 control at ", v),
      contrast_class = "cohort2_treatment_effect",
      visit_f = v,
      baseline_visit = baseline_visit
    )
    
    out[[length(out) + 1]] <- .calc_contrast(
      w_c2 - w_c1, X, beta, V,
      contrast = paste0("Cohort2 treatment effect minus cohort1 treatment effect at ", v),
      contrast_class = "cohort2_minus_cohort1_treatment_effect",
      visit_f = v,
      baseline_visit = baseline_visit
    )
    
    if (!is.na(baseline_visit) && v != baseline_visit) {
      
      w_c1_base <- te_c1(baseline_visit)
      w_c2_base <- te_c2(baseline_visit)
      
      w_c1_did <- w_c1 - w_c1_base
      w_c2_did <- w_c2 - w_c2_base
      
      out[[length(out) + 1]] <- .calc_contrast(
        w_c1_did, X, beta, V,
        contrast = paste0("Cohort1 DID: ", v, " vs ", baseline_visit),
        contrast_class = "cohort1_did",
        visit_f = v,
        baseline_visit = baseline_visit
      )
      
      out[[length(out) + 1]] <- .calc_contrast(
        w_c2_did, X, beta, V,
        contrast = paste0("Cohort2 DID: ", v, " vs ", baseline_visit),
        contrast_class = "cohort2_did",
        visit_f = v,
        baseline_visit = baseline_visit
      )
      
      out[[length(out) + 1]] <- .calc_contrast(
        w_c2_did - w_c1_did, X, beta, V,
        contrast = paste0("Cohort2 DID minus cohort1 DID: ", v, " vs ", baseline_visit),
        contrast_class = "cohort2_minus_cohort1_did",
        visit_f = v,
        baseline_visit = baseline_visit
      )
    }
  }
  
  bind_rows(out) %>%
    .add_fit_meta(dat_lab, fit, cluster_var)
}
fit_cont_trt_lmm <- function(dat_lab,cluster_var,df_time,eval_days,baseline_day,covars,fml) {
  
  fit <- lme4::lmer(
    fml,
    data = dat_lab,
    REML = FALSE
  )
  
  beta <- lme4::fixef(fit)
  V <- .get_vcov(fit, dat_lab, cluster_var)
  
  mf <- model.frame(fit)
  pred_days <- sort(unique(c(eval_days, baseline_day)))
  
  newdat <- expand.grid(
    trt_f = levels(mf$trt_f),
    day_rel_car = pred_days,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  newdat <- .add_reference_covars(newdat, fit, covars)
  X <- .get_fixed_X(fit, newdat)
  
  w_cell <- function(trt_value, day_value) {
    .cell_w(newdat, "trt_f", trt_value, "day_rel_car", day_value)
  }
  
  out <- list()
  
  for (d in eval_days) {
    
    w_trt <- w_cell("anakinra", d) - w_cell("control", d)
    
    out[[length(out) + 1]] <- .calc_contrast(
      w_trt, X, beta, V,
      contrast = paste0("Anakinra vs control at day ", d),
      contrast_class = "treatment_effect",
      day_rel_car = d,
      baseline_day = baseline_day
    )
    
    if (abs(d - baseline_day) > 1e-8) {
      
      w_base <- w_cell("anakinra", baseline_day) -
        w_cell("control", baseline_day)
      
      out[[length(out) + 1]] <- .calc_contrast(
        w_trt - w_base, X, beta, V,
        contrast = paste0("DID: day ", d, " vs day ", baseline_day),
        contrast_class = "did_vs_baseline",
        day_rel_car = d,
        baseline_day = baseline_day
      )
    }
  }
  
  bind_rows(out) %>%
    .add_fit_meta(dat_lab, fit, cluster_var)

}
fit_cont_schedule_lmm <- function(dat_lab,cluster_var,df_time,eval_days,baseline_day,fml,covars,cohort_var) {
  fit <- lme4::lmer(
    fml,
    data = dat_lab,
    REML = FALSE
  )
  
  beta <- lme4::fixef(fit)
  V <- .get_vcov(fit, dat_lab, cluster_var)
  
  mf <- model.frame(fit)
  pred_days <- sort(unique(c(eval_days, baseline_day)))
  
  newdat <- expand.grid(
    group = levels(mf[[cohort_var]]),
    day_rel_car = pred_days,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  names(newdat)[names(newdat) == "group"] <- cohort_var
  
  newdat <- .add_reference_covars(newdat, fit, covars)
  X <- .get_fixed_X(fit, newdat)
  
  w_cell <- function(group_value, day_value) {
    .cell_w(newdat, cohort_var, group_value, "day_rel_car", day_value)
  }
  
  te_c1 <- function(d) {
    w_cell("cohort1", d) - w_cell("cohort1 control", d)
  }
  
  te_c2 <- function(d) {
    w_cell("cohort2", d) - w_cell("cohort2 control", d)
  }
  
  out <- list()
  
  for (d in eval_days) {
    
    w_c1 <- te_c1(d)
    w_c2 <- te_c2(d)
    
    out[[length(out) + 1]] <- .calc_contrast(
      w_c1, X, beta, V,
      contrast = paste0("Cohort1 anakinra vs cohort1 control at day ", d),
      contrast_class = "cohort1_treatment_effect",
      day_rel_car = d,
      baseline_day = baseline_day
    )
    
    out[[length(out) + 1]] <- .calc_contrast(
      w_c2, X, beta, V,
      contrast = paste0("Cohort2 anakinra vs cohort2 control at day ", d),
      contrast_class = "cohort2_treatment_effect",
      day_rel_car = d,
      baseline_day = baseline_day
    )
    
    out[[length(out) + 1]] <- .calc_contrast(
      w_c2 - w_c1, X, beta, V,
      contrast = paste0("Cohort2 treatment effect minus cohort1 treatment effect at day ", d),
      contrast_class = "cohort2_minus_cohort1_treatment_effect",
      day_rel_car = d,
      baseline_day = baseline_day
    )
    
    if (abs(d - baseline_day) > 1e-8) {
      
      w_c1_base <- te_c1(baseline_day)
      w_c2_base <- te_c2(baseline_day)
      
      w_c1_did <- w_c1 - w_c1_base
      w_c2_did <- w_c2 - w_c2_base
      
      out[[length(out) + 1]] <- .calc_contrast(
        w_c1_did, X, beta, V,
        contrast = paste0("Cohort1 DID: day ", d, " vs day ", baseline_day),
        contrast_class = "cohort1_did",
        day_rel_car = d,
        baseline_day = baseline_day
      )
      
      out[[length(out) + 1]] <- .calc_contrast(
        w_c2_did, X, beta, V,
        contrast = paste0("Cohort2 DID: day ", d, " vs day ", baseline_day),
        contrast_class = "cohort2_did",
        day_rel_car = d,
        baseline_day = baseline_day
      )
      
      out[[length(out) + 1]] <- .calc_contrast(
        w_c2_did - w_c1_did, X, beta, V,
        contrast = paste0("Cohort2 DID minus cohort1 DID: day ", d, " vs day ", baseline_day),
        contrast_class = "cohort2_minus_cohort1_did",
        day_rel_car = d,
        baseline_day = baseline_day
      )
    }
  }
  
  bind_rows(out) %>%
    .add_fit_meta(dat_lab, fit, cluster_var)
}
fit_cont_cohort_lmm <- function(dat_lab,df_time,eval_days,baseline_day,fml,covars,cohort_var) {
  fit <- lme4::lmer(
    fml,
    data = dat_lab,
    REML = FALSE
  )
  
  beta <- lme4::fixef(fit)
  V <- .get_vcov(fit, dat_lab, cluster_var = NULL)
  
  mf <- model.frame(fit)
  pred_days <- sort(unique(c(eval_days, baseline_day)))
  
  newdat <- expand.grid(
    group = levels(mf[[cohort_var]]),
    day_rel_car = pred_days,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  names(newdat)[names(newdat) == "group"] <- cohort_var
  
  newdat <- .add_reference_covars(newdat, fit, covars)
  X <- .get_fixed_X(fit, newdat)
  
  w_cell <- function(group_value, day_value) {
    .cell_w(newdat, cohort_var, group_value, "day_rel_car", day_value)
  }
  
  out <- list()
  
  for (d in eval_days) {
    
    w_c2_vs_c1 <- w_cell("cohort2", d) - w_cell("cohort1", d)
    
    out[[length(out) + 1]] <- .calc_contrast(
      w_c2_vs_c1, X, beta, V,
      contrast = paste0("Cohort2 vs cohort1 at day ", d),
      contrast_class = "cohort2_vs_cohort1",
      day_rel_car = d,
      baseline_day = baseline_day
    )
    
    if (abs(d - baseline_day) > 1e-8) {
      
      w_base <- w_cell("cohort2", baseline_day) -
        w_cell("cohort1", baseline_day)
      
      out[[length(out) + 1]] <- .calc_contrast(
        w_c2_vs_c1 - w_base, X, beta, V,
        contrast = paste0("Cohort2 vs cohort1 DID: day ", d, " vs day ", baseline_day),
        contrast_class = "cohort2_vs_cohort1_did",
        day_rel_car = d,
        baseline_day = baseline_day
      )
    }
  }
  
  bind_rows(out) %>%
    .add_fit_meta(dat_lab, fit, cluster_var = NULL) %>%
    mutate(
      n_cohort1 = n_distinct(dat_lab$record_id[dat_lab[[cohort_var]] == "cohort1"]),
      n_cohort2 = n_distinct(dat_lab$record_id[dat_lab[[cohort_var]] == "cohort2"])
    )
}

# Wrapper function for LMM
run_labs <- function(data, fit_fun, ...) {
  
  labs <- sort(unique(as.character(data$lab)))
  
   purrr::map_dfr(
    labs,
    function(x) {
      
      dat_lab <- data %>%
        filter(lab == x)
      
      res <- fit_fun(dat_lab, ...)
      
      if (is.null(res)) {
        return(NULL)
      }
      
      res %>%
        mutate(lab = as.character(x), .before = 1)
    }
  ) %>%
    group_by(contrast) %>%
    mutate(q = p.adjust(p, method = "BH")) %>%
    ungroup() %>%
    arrange(contrast, q)
  
}

# Plot cytokine/lab trajectories. 
plot_lab_trajectory_cont_trt <- function(lab_name) {
  
  plot_dat <- dat1 %>%
    filter(lab == lab_name)
  
  ggplot(plot_dat, aes(x = day_rel_car, y = npx, color=trt_f)) +
    #geom_line(alpha = 0.25) +
    #geom_point(alpha = 0.45) +
    stat_summary(
      aes(group = trt_f),
      fun = mean,
      geom = "line",
      linewidth = 1.2
    ) +
    stat_summary(
      aes(group = trt_f),
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.15
    ) +
    #facet_wrap(~ trt_f) +
    labs(
      x = NULL,
      y = "NPX",
      title = paste("Longitudinal NPX:", lab_name)
    ) +
    theme_classic()
}
plot_lab_trajectory_cont_cohort <- function(data, lab_name) {
  
  dat1 <- data 
  plot_dat <- dat1 %>%
    filter(lab == lab_name) %>%
    filter(!grepl('control', akra_cohort_f))
  
  ggplot(plot_dat, aes(x = day_rel_car, y = npx, color=akra_cohort_f)) +
    #geom_line(alpha = 0.25) +
    #geom_point(alpha = 0.45) +
    stat_summary(
      aes(group = akra_cohort_f),
      fun = mean,
      geom = "line",
      linewidth = 1.2
    ) +
    stat_summary(
      aes(group = akra_cohort_f),
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.15
    ) +
    labs(
      x = NULL,
      y = "NPX",
      title = paste("Longitudinal NPX:", lab_name)
    ) +
    theme_classic()
}
plot_lab_trajectory_cont_cohort_grid <- function(lab_name) {
  
  plot_dat <- dat1 %>%
    filter(lab == lab_name)
  
  ggplot(plot_dat, aes(x = day_rel_car, y = npx, color=trt_f)) +
    #geom_line(alpha = 0.25) +
    #geom_point(alpha = 0.45) +
    stat_summary(
      aes(group = trt_f),
      fun = mean,
      geom = "line",
      linewidth = 1.2
    ) +
    stat_summary(
      aes(group = trt_f),
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.15
    ) +
    facet_wrap(~ match_groups) +
    labs(
      x = NULL,
      y = "NPX",
      title = paste("Longitudinal NPX:", lab_name)
    ) +
    theme_classic()
}

