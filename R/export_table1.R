#!/usr/bin/env Rscript
# Current Table 1 rebuilt from raw-derived cohorts, not existing table cells.
arg <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
HERE <- dirname(normalizePath(arg))
source(file.path(HERE,'model_common.R'))
OUT <- file.path(WORK,'displays')
dir.create(OUT,showWarnings=FALSE)
if(file.exists(file.path(OUT,'table_1_two_panel_characteristics.csv'))) stop('Refusing to overwrite Table 1')
design0304 <- subset(design(groups[['2003-2004']]),primary_complete_case)
design0506 <- subset(design(groups[['2005-2006']]),primary_complete_case)
design_pooled <- subset(design(groups[['hip']]),primary_complete_case)
designs <- list('2011-2012'=subset(design(groups[['2011-2012']]),primary_complete_case),
  '2013-2014'=subset(design(groups[['2013-2014']]),primary_complete_case),
  'Pooled 2011-2014'=subset(design(groups[['wrist']]),primary_complete_case))
fmt_mean_sd <- function(design, variable, digits = 1) {
  f <- as.formula(paste0("~", variable))
  mean_value <- as.numeric(coef(svymean(f, design, na.rm = TRUE)))
  sd_value <- sqrt(as.numeric(coef(svyvar(f, design, na.rm = TRUE))))
  sprintf(paste0("%.", digits, "f ± %.", digits, "f"), mean_value, sd_value)
}

fmt_category <- function(design, variable, level) {
  d <- design$variables
  n <- sum(d[[variable]] == level, na.rm = TRUE)
  indicator <- as.numeric(d[[variable]] == level)
  proportion <- as.numeric(coef(svymean(~indicator, design, na.rm = TRUE)))
  sprintf("%s (%.1f%%)", format(n, big.mark = ","), 100 * proportion)
}

table1 <- data.frame(
  Characteristic = c(
    "Participants, n",
    "Deaths, n",
    "Age, years",
    "Female sex, n (weighted %)",
    "Race and ethnicity, n (weighted %)",
    "  Hispanic",
    "  Non-Hispanic White",
    "  Non-Hispanic Black",
    "  Other race and ethnicity",
    "Education, n (weighted %)",
    "  More than high school",
    "  High school or GED",
    "  Less than high school",
    "Poverty-income ratio",
    "Smoking, n (weighted %)",
    "  Never",
    "  Former",
    "  Current",
    "Self-reported activity index, units/day",
    "Device activity, counts per wear minute",
    "Mortality follow-up, years"
  ),
  stringsAsFactors = FALSE
)

fill_column <- function(design) {
  d <- design$variables
  c(
    format(nrow(d), big.mark = ","),
    format(sum(d$death), big.mark = ","),
    fmt_mean_sd(design, "age"),
    fmt_category(design, "sex", "female"),
    "",
    fmt_category(design, "race4", "Hispanic"),
    fmt_category(design, "race4", "Non-Hispanic White"),
    fmt_category(design, "race4", "Non-Hispanic Black"),
    fmt_category(design, "race4", "Other"),
    "",
    fmt_category(design, "education3", "more than high school"),
    fmt_category(design, "education3", "high school/GED"),
    fmt_category(design, "education3", "under high school"),
    fmt_mean_sd(design, "PIR"),
    "",
    fmt_category(design, "smoking3", "never"),
    fmt_category(design, "smoking3", "former"),
    fmt_category(design, "smoking3", "current"),
    fmt_mean_sd(design, "SR_day"),
    fmt_mean_sd(design, "PAX_counts_per_wear_minute__primary"),
    fmt_mean_sd(design, "follow_up_years")
  )
}

table1[["2003-2004"]] <- fill_column(design0304)
table1[["2005-2006"]] <- fill_column(design0506)
table1[["Pooled 2003-2006"]] <- fill_column(design_pooled)

mean_sd <- fmt_mean_sd
factor_value_percent <- function(design, variable, level) {
  indicator <- as.numeric(design$variables[[variable]] == level)
  estimate <- coef(svymean(~indicator, design, na.rm = TRUE))[1] * 100
  n <- sum(design$variables[[variable]] == level, na.rm = TRUE)
  sprintf("%s (%.1f%%)", format(n, big.mark = ","), estimate)
}

wrist_panel <- data.frame(
  Characteristic = c(
    "Participants, n",
    "Deaths, n",
    "Age, years",
    "Female sex, n (weighted %)",
    "Race and ethnicity, n (weighted %)",
    "  Hispanic",
    "  Non-Hispanic White",
    "  Non-Hispanic Black",
    "  Other race and ethnicity",
    "Education, n (weighted %)",
    "  More than high school",
    "  High school or GED",
    "  Less than high school",
    "Poverty-income ratio",
    "Smoking, n (weighted %)",
    "  Never",
    "  Former",
    "  Current",
    "Self-reported activity-equivalent time, min/day",
    "Device activity, mean daily triaxial MIMS",
    "Mortality follow-up, years"
  ),
  stringsAsFactors = FALSE
)

fill_wrist_panel <- function(design) {
  x <- design$variables
  c(
    format(nrow(x), big.mark = ","),
    format(sum(x$death), big.mark = ","),
    mean_sd(design, "age"),
    factor_value_percent(design, "sex", "female"),
    "",
    factor_value_percent(design, "race4", "Hispanic"),
    factor_value_percent(design, "race4", "Non-Hispanic White"),
    factor_value_percent(design, "race4", "Non-Hispanic Black"),
    factor_value_percent(design, "race4", "Other"),
    "",
    factor_value_percent(design, "education3", "more than high school"),
    factor_value_percent(design, "education3", "high school/GED"),
    factor_value_percent(design, "education3", "under high school"),
    mean_sd(design, "PIR"),
    "",
    factor_value_percent(design, "smoking3", "never"),
    factor_value_percent(design, "smoking3", "former"),
    factor_value_percent(design, "smoking3", "current"),
    mean_sd(design, "SR_day"),
    mean_sd(design, "PAX_mean_daily_mims", 0),
    mean_sd(design, "follow_up_years")
  )
}

wrist_panel[["2011-2012"]] <- fill_wrist_panel(designs[["2011-2012"]])
wrist_panel[["2013-2014"]] <- fill_wrist_panel(designs[["2013-2014"]])
wrist_panel[["Pooled 2011-2014"]] <-
  fill_wrist_panel(designs[["Pooled 2011-2014"]])

hip_panel <- table1
stopifnot(
  identical(which(hip_panel$Characteristic != wrist_panel$Characteristic), c(19L,20L)),
  hip_panel$Characteristic[19] == 'Self-reported activity index, units/day',
  wrist_panel$Characteristic[19] == 'Self-reported activity-equivalent time, min/day',
  hip_panel$Characteristic[20] == 'Device activity, counts per wear minute',
  wrist_panel$Characteristic[20] == 'Device activity, mean daily triaxial MIMS'
)

panel_a_label <- data.frame(
  Characteristic =
    "Panel A. Hip-worn uniaxial accelerometry (NHANES 2003-2006)",
  `2003-2004` = "", `2005-2006` = "", `Pooled 2003-2006` = "",
  check.names = FALSE
)
panel_b_label <- data.frame(
  Characteristic =
    "Panel B. Wrist-worn triaxial accelerometry (NHANES 2011-2014)",
  `2011-2012` = "", `2013-2014` = "", `Pooled 2011-2014` = "",
  check.names = FALSE
)
panel_a_header <- as.data.frame(
  as.list(names(hip_panel)), check.names = FALSE
)
names(panel_a_header) <- names(hip_panel)
panel_b_header <- as.data.frame(
  as.list(names(wrist_panel)), check.names = FALSE
)
names(panel_b_header) <- names(wrist_panel)

table_1_two_panel <- rbind(
  panel_a_label,
  panel_a_header,
  hip_panel,
  setNames(panel_b_label, names(hip_panel)),
  setNames(panel_b_header, names(hip_panel)),
  setNames(wrist_panel, names(hip_panel))
)
write.csv(
  table_1_two_panel,
  file.path(OUT, "table_1_two_panel_characteristics.csv"),
  row.names = FALSE, quote = TRUE
)
