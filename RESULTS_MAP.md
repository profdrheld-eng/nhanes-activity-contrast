# Result-to-code map

The original 43 display components (37 tables and 6 figures) are listed first. Subsequent sections map S32–S36 and the diagnostic source values for S37. Reporting checklists S38–S39 are in the journal Supplementary File; they are not computational outputs.
The existing supplement numbering is retained; there is no additional S16
or unsuffixed S24 to generate. Each named numerical CSV below is under
`work/models/` unless another folder is shown. Display files go to `work/displays/`.
Main tables are formatted by `export_main_tables.py`, supplement tables by
`export_supplement_tables.py`, and post-hoc Tables S30–S31 by
`export_no_transport.py`. These exporters read freshly computed results,
not the Word manuscript or historical coefficients.

Source chain for numerical results:
`config/sources.csv` → verified local NHANES/NCHS files → hip processing where
applicable → `nhanes_activity/cohorts.py` → four cycle cohorts →
`prepare_models.py` → model/diagnostic stages below → table/figure exporters.
The manifest and DATA_SOURCES.md identify each cycle's exact input modules.
METHODS.md specifies domains, transformations, covariates and inference.

| Document target | Computational producer | Intermediate / figure | Selection or interpretation |
|---|---|---|---|
| Table 1 | `R/export_table1.R` | `cohorts → models/revision_cohort.csv.gz` | mortality complete-case characteristics |
| Table 2 | `R/run_models.R; R/check_regime_nuisance.R` | `mortality_results; regime_nuisance_results` | primary common slopes; all-nuisance regime comparison |
| Table 3 | `R/run_models.R; R/check_regime_nuisance.R` | `mortality_results; regime_nuisance_results` | selected sensitivities; five-year regime comparison |
| Table 4 | `R/run_models.R` | `age_results; joint_tests` | age 80 versus 50, primary definition |
| Table S1 | `config/method_tables.json` | `config/sources.csv; cohort/questionnaire source variables` | source specification |
| Table S2 | `prepare_models.py` | `participant_flow` | sequential cohort domains |
| Table S3 | `prepare_models.py` | `inclusion_comparison` | included/excluded examined adults |
| Table S4 | `config/method_tables.json` | `nhanes_activity/questionnaires.py` | questionnaire specification |
| Table S5 | `config/method_tables.json` | `nhanes_activity/pax.py; questionnaires.aggregate_wrist` | device specification |
| Table S6 | `prepare_models.py` | `contrast_distribution` | paired reference domain, raw normal-score difference |
| Table S7 | `export_score_diagnostics.py` | `rank_ties_and_bounds` | raw values; log-input diagnostic exported separately |
| Table S8 | `export_score_diagnostics.py` | `contrast_zero_reference` | equal-rank location after centering |
| Table S9 | `prepare_models.py` | `correlations` | reference and complete-case domains |
| Table S10 | `R/run_models.R` | `model_support` | primary common/sex-interaction sample and df |
| Table S11 | `R/run_models.R` | `mortality_results` | non-primary sensitivities, three references |
| Table S12 | `R/run_additional_sensitivities.R; R/response_weights.R` | `additional_sensitivity_results` | device/response/time-structure variants |
| Table S13 | `R/check_regime_nuisance.R` | `regime_nuisance_results` | two nuisance specifications, horizons, three references |
| Table S14 | `R/run_models.R; R/check_regime_nuisance.R` | `joint_tests; regime_nuisance_tests` | joint Wald tests |
| Table S15 | `R/run_retained_sensitivities.R; export_retained.py` | `legacy_sensitivities_unified` | retained 2005–2006 sensitivity definitions |
| Table S17 | `R/run_diagnostics.R` | `ph_diagnostics` | conventional weighted Schoenfeld |
| Table S18 | `R/run_diagnostics.R` | `time_varying` | survey step-time and cluster-robust log-time |
| Table S19 | `R/run_diagnostics.R` | `diagnostic_support` | conditioning and approximate influence |
| Table S19a | `R/run_diagnostics.R` | `influence` | coefficient-specific approximate influence |
| Table S20 | `R/run_models.R` | `spline_support` | weighted knots/ranges and tail n/deaths |
| Table S21 | `R/run_models.R` | `age_results` | age endpoint contrasts, all definitions |
| Table S22 | `R/check_age_support.R` | `age_support_sensitivity` | age endpoints/domains |
| Table S23 | `R/age_sex_overlay.R` | `age_sex_overlay` | pooled overlay sample/model specification |
| Table S24a | `R/run_exploratory.R` | `exploratory/mortality_modifiers` | all inference references |
| Table S24b | `R/run_exploratory.R` | `exploratory/participant_correlates` | all inference references |
| Table S24c | `R/run_exploratory.R` | `exploratory/participant_omnibus` | all inference references |
| Table S25 | `R/run_exploratory.R` | `exploratory/mortality_modifiers` | full-design t and modifier BH family |
| Table S26 | `R/run_exploratory.R` | `exploratory/participant_omnibus` | full-design F and characteristic BH family |
| Table S27 | `R/run_exploratory.R` | `exploratory/participant_correlates` | selected core and extension contrasts |
| Table S28 | `R/run_models.R` | `mortality_results` | primary common/interaction, three references |
| Table S29 | `R/run_retained_sensitivities.R; export_retained.py` | `wrist_mortality_sensitivity_results` | retained wrist sensitivities |
| Figure 1 | `R/plot_density.R` | `models/revision_cohort.csv.gz → Figure_1_two_regime_density.png` | paired reference sample by regime/sex |
| Figure 2 | `R/plot_age_main.R` | `age_sex_overlay → Figure_2_age_overlay.png` | pooled common and exploratory sex-specific age curves |
| Figure S1 | `R/plot_flow.R` | `participant_flow → participant_flow_sankey.png` | four cycle flow counts |
| Figure S2 | `R/plot_supplement.R` | `spline_curves; spline_support; joint_tests → mortality_spline_common.png` | common slope curves |
| Figure S3 | `R/plot_supplement.R` | `spline_curves; spline_support; joint_tests → mortality_spline_interaction.png` | sex-specific curves |
| Figure S4 | `R/plot_supplement.R` | `age_curves → age_curves_alternative_definitions.png` | cycle-specific alternative definitions |

The main mortality claims use the `__common__primary` full-design-t rows,
and the direct regime comparison uses `all_nuisance_era__Inf`. The five-year
regime sensitivity uses `all_nuisance_era__5`. Age endpoint claims use
`__age__primary` and their global age tests. Their displayed estimates,
intervals, p-values, sample/event counts and degrees of freedom are covered by
Tables 2–4; unrounded numerical results remain in the local model exports.

The main-text wrist response-balance statement is regenerated by
`R/response_balance.R` from the freshly fitted response-model domain and
weights. It exports `response/wrist/balance.csv` and `balance_summary.csv`.
The target comprises 9,844 eligible complete-predictor records; 7,814 responders
are compared with that target before and after weighting. Standardizers are
the target weighted population SD for continuous variables and target
Bernoulli SD for categorical indicators. Maximum absolute SMDs are 0.060310
before and 0.001599 after weighting (reported as 0.060 and 0.002). This describes
observed covariate balance, not guaranteed absence of selection bias.

## Additional post-hoc transport sensitivity

`R/run_no_transport.R` reads the freshly prepared model cohort, calls
`R/no_transport_scores.R` within the original hip-cycle reference domains,
and refits the original and no-transport common/sex mortality models and
80-versus-50 age contrasts. `export_no_transport.py` writes Tables S30/S31
to `displays/no_transport_tables.json` and individual table CSVs.
`models/no_transport_score_summary.csv` records domain sizes and weighted
score correlations. `models/no_transport_scores.csv.gz` is a local
participant-level verification file and must never be published.

## Reporting additions

- S32–S33: `prepare_reporting.py`, reporting/missingness.csv and criterion_failures.csv. Four explicit domains; nonsequential counts.
- S34: reporting/followup.csv, unweighted principal-case person-years and follow-up summaries.
- S35: reporting/quintile_cutpoints.csv and quintile_counts.csv, original reference-domain weighted boundaries and membership verification.
- S36: `R/run_reporting.R`, reporting/mortality_adjustment.csv and mortality_support.csv. Six samples, three adjustment levels; 18 fits with principal-result parity checks.
- All five formatted tables and notes: `export_reporting.py`, displays/reporting_tables.json and reporting_table_notes.json.

| Table S37 | `R/run_reporting_diagnostics.R` | `reporting_diagnostics/{schoenfeld,step_time,support}.csv` | all 18 S36 models; conventional and full-survey tests kept distinct |
