#!/usr/bin/env Rscript
fit_script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
source(file.path(dirname(normalizePath(fit_script)),'fit_checks.R'))
# Refit the retained 2005-2006 exploratory analyses from the new raw-data cohort.
# Historical, data-informed nonlinear selection is reproduced transparently.
# These are not prespecified confirmatory analyses. BH families remain separate:
# four mortality modifiers, three nonlinear screens, eight participant attributes.
options(stringsAsFactors=FALSE,survey.lonely.psu='adjust')
suppressPackageStartupMessages({library(survey);library(survival)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop('Supply exactly one local work directory')
WORK <- normalizePath(args[1],mustWork=TRUE)
COHORT <- file.path(WORK,'2005-2006','cohort','analytic_cohort.csv.gz')
OUT <- file.path(WORK,'exploratory')
if(dir.exists(OUT)) stop('Refusing to overwrite exploratory results')
if(!dir.create(OUT)) stop('Cannot create exploratory result directory')
scalar <- function(b,se,df,exponentiate=FALSE) {
  stopifnot(is.finite(b),is.finite(se),se>0,df>0)
  critical <- if(is.infinite(df)) qnorm(.975) else qt(.975,df)
  p <- if(is.infinite(df)) 2*pnorm(-abs(b/se)) else 2*pt(-abs(b/se),df)
  tr <- if(exponentiate) exp else identity
  data.frame(estimate=tr(b),ci_low=tr(b-critical*se),ci_high=tr(b+critical*se),p=p)
}
DESIGN_DF <- 15
as_flag <- function(x) {
  y <- tolower(as.character(x))
  !is.na(y) & y %in% c("true", "1", "yes")
}

weighted_mean <- function(x, w) sum(x * w) / sum(w)

fmt_p <- function(x) {
  ifelse(x < 0.001, "<.001", sub("^0", "", formatC(x, digits = 3, format = "f")))
}

d <- read.csv(gzfile(COHORT), check.names = FALSE)
for (nm in c("primary_complete_case", "reference_domain")) d[[nm]] <- as_flag(d[[nm]])
d$sex <- factor(d$sex, levels = c("male", "female"))
d$race4 <- factor(
  d$race_ethnicity4,
  levels = c("Non-Hispanic White", "Hispanic", "Non-Hispanic Black", "Other")
)
d$education3 <- factor(
  d$education3,
  levels = c("more than high school", "high school/GED", "under high school")
)
d$smoking3 <- factor(d$smoking3, levels = c("never", "former", "current"))
d$death <- as.integer(d$MORTSTAT)

primary <- d[d$primary_complete_case, ]
stopifnot(nrow(primary) == 2737, sum(primary$death) == 505)

age_center <- weighted_mean(primary$age, primary$WTMEC2YR)
pir_center <- weighted_mean(primary$PIR, primary$WTMEC2YR)
bmi_center <- weighted_mean(
  primary$BMI[is.finite(primary$BMI)],
  primary$WTMEC2YR[is.finite(primary$BMI)]
)

d$age10_c <- (d$age - age_center) / 10
d$PIR_c <- d$PIR - pir_center
d$BMI5_c <- (d$BMI - bmi_center) / 5

health_vars <- c("diabetes_history", "CVD_history", "cancer_history")
health <- d[, health_vars, drop = FALSE]
health[health == ""] <- NA_character_
d$any_baseline_disease <- NA_character_
d$any_baseline_disease[apply(health == "yes", 1, any, na.rm = TRUE)] <- "yes"
all_no <- apply(health == "no", 1, all)
d$any_baseline_disease[all_no] <- "no"
d$any_baseline_disease <- factor(d$any_baseline_disease, levels = c("no", "yes"))

base_formula <- paste(
  "Surv(follow_up_years, death) ~ Delta * sex + Level + age +",
  "race4 + education3 + PIR + smoking3"
)

specs <- list(
  age = list(
    label = "Age",
    unit = "per 10-year higher age",
    term = "Delta:age10_c",
    extra = "Delta:age10_c",
    domain = d$primary_complete_case,
    values = c("40 years" = (40 - age_center) / 10, "60 years" = (60 - age_center) / 10)
  ),
  bmi = list(
    label = "Body mass index",
    unit = "per 5-kg/m² higher BMI",
    term = "Delta:BMI5_c",
    extra = "BMI5_c + Delta:BMI5_c",
    domain = d$primary_complete_case & is.finite(d$BMI5_c),
    values = c("BMI 25 kg/m²" = (25 - bmi_center) / 5, "BMI 30 kg/m²" = (30 - bmi_center) / 5)
  ),
  pir = list(
    label = "Poverty-income ratio",
    unit = "per 1-unit higher ratio",
    term = "Delta:PIR_c",
    extra = "Delta:PIR_c",
    domain = d$primary_complete_case,
    values = c("Ratio 1.3" = 1.3 - pir_center, "Ratio 3.5" = 3.5 - pir_center)
  ),
  health = list(
    label = "Any baseline disease",
    unit = "yes versus no",
    term = "Delta:any_baseline_diseaseyes",
    extra = "any_baseline_disease + Delta:any_baseline_disease",
    domain = d$primary_complete_case & !is.na(d$any_baseline_disease),
    values = c("No baseline disease" = 0, "Any baseline disease" = 1)
  )
)

modifier_rows <- list();modifier_support <- list()
for(id in names(specs)) {
  spec <- specs[[id]]
  dd <- d;dd$analysis_domain <- spec$domain
  des <- subset(svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC2YR,nest=TRUE,data=dd),analysis_domain)
  fit <- checked_model(svycoxph(as.formula(paste(base_formula,'+',spec$extra)),design=des,method='efron'))
  stopifnot(is.null(fit$fail),all(is.finite(coef(fit))),all(is.finite(vcov(fit))))
  b <- unname(coef(fit)[spec$term]);se <- unname(sqrt(vcov(fit)[spec$term,spec$term]))
  rdf <- fit$degf.resid;ddf <- degf(des)
  modifier_support[[id]] <- data.frame(modifier=id,formula=paste(deparse(formula(fit)),collapse=' '),
     n=nrow(des$variables),deaths=sum(des$variables$death),design_df=ddf,residual_df=rdf)
  for(conv in c('full_design_t','residual_t','normal')) {
    df <- switch(conv,full_design_t=ddf,residual_t=rdf,normal=Inf)
    modifier_rows[[length(modifier_rows)+1]] <- cbind(data.frame(modifier=spec$label,scale=spec$unit,
      n=nrow(des$variables),deaths=sum(des$variables$death),convention=conv,denominator_df=df,
      design_df=ddf,residual_df=rdf,coefficient=b,design_based_se=se),scalar(b,se,df,TRUE))
  }
}
modifier_rows <- do.call(rbind,modifier_rows)
modifier_rows$q_bh <- ave(modifier_rows$p,modifier_rows$convention,FUN=function(p)p.adjust(p,'BH'))
write.csv(modifier_rows,file.path(OUT,'mortality_modifiers.csv'),row.names=FALSE)
write.csv(do.call(rbind,modifier_support),file.path(OUT,'modifier_support.csv'),row.names=FALSE)
write.csv(data.frame(age=age_center,PIR=pir_center,BMI=bmi_center),file.path(OUT,'modifier_centers.csv'),row.names=FALSE)

# Participant-level cross-sectional exploratory correlates.
DESIGN_DF <- 15
CHARACTERISTIC_TESTS <- c(
  "age", "sex", "race_ethnicity", "education", "PIR",
  "BMI", "smoking", "baseline_disease"
)

as_flag <- function(x) {
  y <- tolower(as.character(x))
  !is.na(y) & y %in% c("true", "1", "yes")
}

weighted_mean <- function(x, w) sum(x * w) / sum(w)

weighted_sd <- function(x, w) {
  mu <- weighted_mean(x, w)
  sqrt(sum(w * (x - mu)^2) / sum(w))
}

fmt_p <- function(x) {
  ifelse(x < 0.001, "<.001", sub("^0", "", formatC(x, digits = 3, format = "f")))
}

make_disease <- function(d) {
  h <- d[, c("diabetes_history", "CVD_history", "cancer_history"), drop = FALSE]
  h[h == ""] <- NA_character_
  out <- rep(NA_character_, nrow(d))
  out[apply(h == "yes", 1, any, na.rm = TRUE)] <- "yes"
  out[apply(h == "no", 1, all)] <- "no"
  factor(out, levels = c("no", "yes"))
}

set_variables <- function(d) {
  d$primary_complete_case <- as_flag(d$primary_complete_case)
  d$sex <- factor(d$sex, levels = c("male", "female"))
  d$race4 <- factor(
    d$race_ethnicity4,
    levels = c("Non-Hispanic White", "Hispanic", "Non-Hispanic Black", "Other")
  )
  d$education3 <- factor(
    d$education3,
    levels = c("more than high school", "high school/GED", "under high school")
  )
  d$smoking3 <- factor(d$smoking3, levels = c("never", "former", "current"))
  d$baseline_disease <- make_disease(d)
  d$age10_c50 <- (d$age - 50) / 10
  d$BMI5_c25 <- (d$BMI - 25) / 5
  d$analysis_domain <- d$primary_complete_case &
    is.finite(d$BMI5_c25) & !is.na(d$baseline_disease)
  d
}

make_design <- function(d) {
  svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC2YR,
    nest = TRUE,
    data = d
  )
}

tidy_coefficients <- function(fit, model_id, term_labels) {
  b <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  rdf <- fit$df.residual
  critical <- qt(0.975, rdf)
  p <- 2 * pt(abs(b / se), df = rdf, lower.tail = FALSE)
  keep <- intersect(names(term_labels), names(b))
  data.frame(
    model_id = model_id,
    term = keep,
    label = unname(term_labels[keep]),
    estimate_delta_sd = unname(b[keep]),
    design_based_se = unname(se[keep]),
    ci_low = unname(b[keep] - critical * se[keep]),
    ci_high = unname(b[keep] + critical * se[keep]),
    p_value = unname(p[keep]),
    residual_df = rdf,
    stringsAsFactors = FALSE
  )
}

term_test <- function(fit, model_id, characteristic, formula) {
  test <- regTermTest(fit, formula)
  data.frame(
    model_id = model_id,
    characteristic = characteristic,
    numerator_df = unname(test$df),
    denominator_df = unname(test$ddf),
    F_statistic = as.numeric(test$Ftest),
    p_value = as.numeric(test$p),
    stringsAsFactors = FALSE
  )
}

weighted_r2 <- function(fit, domain) {
  y <- domain$variables$Delta
  w <- weights(domain)
  mu <- weighted_mean(y, w)
  1 - sum(w * (y - fitted(fit))^2) / sum(w * (y - mu)^2)
}

contrast_row <- function(fit, label, L, reference = "Age 50 years") {
  b <- coef(fit)
  v <- vcov(fit)
  full <- setNames(rep(0, length(b)), names(b))
  full[names(L)] <- L
  estimate <- sum(full * b)
  se <- sqrt(as.numeric(t(full) %*% v %*% full))
  rdf <- fit$df.residual
  critical <- qt(0.975, rdf)
  data.frame(
    label = label,
    reference = reference,
    estimate_delta_sd = estimate,
    design_based_se = se,
    ci_low = estimate - critical * se,
    ci_high = estimate + critical * se,
    p_value = 2 * pt(abs(estimate / se), df = rdf, lower.tail = FALSE),
    residual_df = rdf,
    stringsAsFactors = FALSE
  )
}

d <- set_variables(read.csv(gzfile(COHORT), check.names = FALSE))
design_full <- make_design(d)
domain <- subset(design_full, analysis_domain)
stopifnot(
  nrow(domain$variables) == 2721,
  degf(domain) == DESIGN_DF,
  length(unique(domain$variables$SDMVSTRA)) == 15,
  nrow(unique(domain$variables[c("SDMVSTRA", "SDMVPSU")])) == 30
)

core_formula_linear <- Delta ~
  Level + age10_c50 + sex + race4 + education3 + PIR

nonlinear_specs <- list(
  age = list(
    formula = update(core_formula_linear, . ~ . + I(age10_c50^2)),
    term = ~I(age10_c50^2)
  ),
  PIR = list(
    formula = update(core_formula_linear, . ~ . + I(PIR^2)),
    term = ~I(PIR^2)
  ),
  BMI = list(
    formula = update(core_formula_linear, . ~ . + BMI5_c25 + I(BMI5_c25^2)),
    term = ~I(BMI5_c25^2)
  )
)

nonlinear_rows <- lapply(names(nonlinear_specs), function(id) {
  spec <- nonlinear_specs[[id]]
  fit <- checked_model(svyglm(spec$formula, design = domain))
  test <- regTermTest(fit, spec$term)
  data.frame(
    variable = id,
    quadratic_estimate = unname(tail(coef(fit), 1)),
    quadratic_se = unname(tail(sqrt(diag(vcov(fit))), 1)),
    denominator_df = unname(test$ddf),
    F_statistic = as.numeric(test$Ftest),
    p_value = as.numeric(test$p),
    stringsAsFactors = FALSE
  )
})
nonlinear <- do.call(rbind, nonlinear_rows)
nonlinear$q_value_bh <- p.adjust(nonlinear$p_value, method = "BH")
nonlinear$retained <- nonlinear$variable == "age" &
  nonlinear$q_value_bh < 0.05

stopifnot(
  nonlinear$retained[nonlinear$variable == "age"],
  !any(nonlinear$retained[nonlinear$variable != "age"])
)

core_formula <- update(core_formula_linear, . ~ . + I(age10_c50^2))
formulas <- list(
  core = core_formula,
  bmi_extension = update(core_formula, . ~ . + BMI5_c25),
  smoking_extension = update(core_formula, . ~ . + smoking3),
  disease_extension = update(core_formula, . ~ . + baseline_disease),
  full_sensitivity = update(
    core_formula,
    . ~ . + BMI5_c25 + smoking3 + baseline_disease
  )
)
fits <- lapply(formulas, function(f) checked_model(svyglm(f, design = domain)))

stopifnot(
  fits$core$df.residual == 5,
  fits$bmi_extension$df.residual == 4,
  fits$smoking_extension$df.residual == 3,
  fits$disease_extension$df.residual == 4,
  fits$full_sensitivity$df.residual == 1
)

term_labels <- c(
  "Level" = "Overall activity level, per 1 SD",
  "age10_c50" = "Age linear component, per 10 years",
  "I(age10_c50^2)" = "Age quadratic component",
  "sexfemale" = "Women versus men",
  "race4Hispanic" = "Hispanic versus non-Hispanic White",
  "race4Non-Hispanic Black" = "Non-Hispanic Black versus non-Hispanic White",
  "race4Other" = "Other race and ethnicity versus non-Hispanic White",
  "education3high school/GED" = "High school or GED versus more than high school",
  "education3under high school" = "Under high school versus more than high school",
  "PIR" = "Poverty-income ratio, per 1 unit",
  "BMI5_c25" = "Body mass index, per 5 kg/m²",
  "smoking3former" = "Former versus never smoking",
  "smoking3current" = "Current versus never smoking",
  "baseline_diseaseyes" = "Any baseline disease versus none"
)

coefficients <- do.call(
  rbind,
  lapply(names(fits), function(id) tidy_coefficients(fits[[id]], id, term_labels))
)

omnibus <- rbind(
  term_test(fits$core, "core", "age", ~age10_c50 + I(age10_c50^2)),
  term_test(fits$core, "core", "sex", ~sex),
  term_test(fits$core, "core", "race_ethnicity", ~race4),
  term_test(fits$core, "core", "education", ~education3),
  term_test(fits$core, "core", "PIR", ~PIR),
  term_test(fits$bmi_extension, "bmi_extension", "BMI", ~BMI5_c25),
  term_test(fits$smoking_extension, "smoking_extension", "smoking", ~smoking3),
  term_test(
    fits$disease_extension, "disease_extension",
    "baseline_disease", ~baseline_disease
  )
)
stopifnot(setequal(omnibus$characteristic, CHARACTERISTIC_TESTS))
omnibus$q_value_bh <- p.adjust(omnibus$p_value, method = "BH")

cr <- list();orr <- list()
for(id in names(fits)) {
  fit <- fits[[id]]
  stopifnot(isTRUE(fit$converged),all(is.finite(coef(fit))),all(is.finite(vcov(fit))))
  rows <- coefficients[coefficients$model_id==id,]
  for(i in seq_len(nrow(rows))) {
    r <- rows[i,]
    for(conv in c('full_design_t','residual_t','normal')) {
      df <- switch(conv,full_design_t=degf(domain),residual_t=fit$df.residual,normal=Inf)
      cr[[length(cr)+1]] <- cbind(data.frame(model=id,term=r$term,label=r$label,n=nrow(domain$variables),
        convention=conv,denominator_df=df,design_df=degf(domain),residual_df=fit$df.residual,
        coefficient=r$estimate_delta_sd,design_based_se=r$design_based_se),
        scalar(r$estimate_delta_sd,r$design_based_se,df))
    }
  }
}
for(i in seq_len(nrow(omnibus))) {
  r <- omnibus[i,]
  for(conv in c('full_design_F','residual_F','chi_square')) {
    df <- switch(conv,full_design_F=degf(domain),residual_F=r$denominator_df,chi_square=Inf)
    p <- if(is.infinite(df))pchisq(r$F_statistic*r$numerator_df,r$numerator_df,lower.tail=FALSE) else pf(r$F_statistic,r$numerator_df,df,lower.tail=FALSE)
    orr[[length(orr)+1]] <- data.frame(model=r$model_id,characteristic=r$characteristic,n=nrow(domain$variables),
      F=r$F_statistic,numerator_df=r$numerator_df,denominator_df=df,convention=conv,p=p)
  }
}
orr <- do.call(rbind,orr)
orr$q_bh <- ave(orr$p,orr$convention,FUN=function(p)p.adjust(p,'BH'))
write.csv(do.call(rbind,cr),file.path(OUT,'participant_correlates.csv'),row.names=FALSE)
write.csv(orr,file.path(OUT,'participant_omnibus.csv'),row.names=FALSE)
write.csv(nonlinear,file.path(OUT,'historical_nonlinear_screen.csv'),row.names=FALSE)
support <- do.call(rbind,lapply(names(fits),function(id)data.frame(model=id,n=nrow(domain$variables),
  design_df=degf(domain),residual_df=fits[[id]]$df.residual,formula=paste(deparse(formula(fits[[id]])),collapse=' '))))
write.csv(support,file.path(OUT,'participant_support.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'session_info.txt'))
writeLines('All retained exploratory models refitted successfully.',file.path(OUT,'complete.txt'))
