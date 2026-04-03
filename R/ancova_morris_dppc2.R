ancova_morris_dppc2 <- function(data,
                                pre,
                                post,
                                condition,
                                treat_level   = "treatment",
                                control_level = "control",
                                ci_level      = 0.95,
                                correct       = TRUE) {
  
  require(dplyr)
  require(MDMA)
  require(janitor)
  
  # capture column names as strings
  pre_col       <- rlang::as_name(rlang::ensym(pre))
  post_col      <- rlang::as_name(rlang::ensym(post))
  condition_col <- rlang::as_name(rlang::ensym(condition))
  
  # ensure condition has control_level as reference
  dat <- data %>%
    mutate(
      {{ condition }} := stats::relevel(
        factor({{ condition }}),
        ref = control_level
      )
    )
  
  # ANCOVA: post ~ condition + pre
  form <- stats::as.formula(
    paste(post_col, "~", condition_col, "+", pre_col)
  )
  fit  <- lm(formula = form, data = dat)
  
  coefs <- summary(fit)$coef
  # coefficient for treatment vs control
  coef_name <- paste0(condition_col, treat_level)
  if (!coef_name %in% rownames(coefs)) {
    stop("Could not find coefficient for treatment level '", treat_level, "'.")
  }
  p_val <- coefs[coef_name, "Pr(>|t|)"]
  
  # split data by group
  dat_treat <- dat %>% filter({{ condition }} == treat_level)
  dat_ctrl  <- dat %>% filter({{ condition }} == control_level)
  
  # Morris' d_ppc2 with Hedges' correction
  es <- dPPC2(
    preT    = dat_treat %>% pull({{ pre }}),
    preC    = dat_ctrl  %>% pull({{ pre }}),
    posT    = dat_treat %>% pull({{ post }}),
    posC    = dat_ctrl  %>% pull({{ post }}),
    correct = correct,
    CIlevel = ci_level
  ) %>%
    select(
      d_ppc2   = d,
      ci_lower = lower.bound,
      ci_upper = upper.bound
    ) %>%
    round_half_up(2) %>%
    mutate(p = p_val)
  
  es
}