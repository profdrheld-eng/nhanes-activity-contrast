# Analysis definitions and interpretation

This package reproduces the current NHANES activity-contrast manuscript. Its
analysis history includes choices made after earlier results were available.
It is not presented as a preregistered confirmatory analysis. Statistical
reproduction and defensible interpretation are separate requirements.

## Sources and domains

The source manifest identifies 34 official files, their cycles, sizes and
SHA-256 hashes. The hip period covers 2003–2004 and 2005–2006; the wrist period
covers 2011–2012 and 2013–2014. Public linked mortality is the 2019 release.
No person-level data or fitted objects are distributed with this code.

The examined-adult domain requires age at least 20, examination participation,
positive examination weight and complete masked survey identifiers. Eligible
mortality records have linkage eligibility, known vital status and positive
examination-based follow-up. Follow-up years equal released examination
follow-up months divided by 12. The score-reference domain additionally
requires complete self-report and an eligible device score. The primary
mortality domain then requires age, sex, race/ethnicity, education, poverty-income
ratio and smoking. Source joins enforce their intended key cardinalities.

Complete-case selection does not eliminate selection bias. Reference samples,
mortality model samples and age-model samples differ; exports report each
denominator rather than treating them as interchangeable. Age models preserve
historical additional BMI and disease-completeness restrictions, with separate
domain-sensitivity fits.

## Activity measurements

The hip source contains minute-level uniaxial counts; the requested protocol
used an AM-7164 on the right hip during waking hours. The wrist protocol used
GT3X+ triaxial recording across day and night, preferably on the non-dominant
wrist. These are different measurement regimes, not interchangeable instruments.
See [hip documentation](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAXRAW_D.htm)
and [wrist documentation](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/PAXDAY_G.htm).

The study's hip algorithm retains valid calibrated minutes with counts in
[0,32767). Nonwear is a candidate run of at least 60 minutes: zero counts may
contain one- or two-minute runs strictly between 0 and 100 counts, immediately
bounded by zeros. Missing/invalid positions and day boundaries split runs.
This is not a requirement for 60 uninterrupted zero-count minutes. A valid
day has at least 600 wear minutes, and the primary participant definition
requires four valid days. Counts per wear minute use the ratio of summed
counts to summed wear time over valid days. Steps are ancillary. The code
also reconstructs three-day, five-day, calibration-code 1-or-2 and 90-minute
nonwear alternatives, with new reference samples and scores.

The wrist analysis uses released daily summaries, not a new reconstruction
from 80-Hz recordings. A valid day has PAXWWMD≥600, PAXVMD>0 and PAXMTSD≥0.
The primary score averages PAXMTSD across at least four valid days. PAXMTSD
summarizes quality-valid minutes, including valid minutes outside wake wear;
it is **not a wake-only total**. Alternatives use one valid day or summed MIMS
divided by summed valid minutes. The official
[daily-variable definitions](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/PAXDAY_G.htm)
distinguish wake-wear minutes, valid minutes and total daily MIMS.

Questionnaire formulas and bounds are in `questionnaires.py` and Tables S1/S4
in `config/method_tables.json`. Hip scoring combines transport, home/yard,
moderate leisure and twice vigorous leisure, expressed on a 30-day basis
and divided by 30. Wrist scoring combines the five work/transport/recreation
domains on a weekly basis, doubles vigorous contributions and divides by seven.
Negative gates yield zero; unable, refused, unknown or incomplete positive
responses remain missing. Contradictory hip leisure gates and activity rows
fail explicitly. Multiple leisure rows are expected and summed, not deduplicated.
The [PAQIAF codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAQIAF_D.htm)
documents intensity coding and repeated activity records.

These are operational activity-equivalent indices, not validated absolute
duration estimates. For example, the hip transport formula multiplies the
converted frequency by PAD080, while the questionnaire describes PAD080 as
total minutes on active days and records frequency as times. Repeated daily
events can therefore receive additional weight. The script preserves the
reported study definition, not an assertion that this multiplication measures
true daily duration. The transport minimum of ten minutes is also a study
restriction: the [official PAQ codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAQ_D.htm)
contains shorter released responses. Changing these definitions would change
the estimand and requires a separately reported sensitivity/revision.

## Rank construction and numerical precision

Within each cycle's paired reference domain, self-report and device values
are log1p-transformed, converted to examination-weighted tied midranks and then
inverse-normal scores. A tied group's percentile is its preceding weight plus
half its own weight, divided by total reference weight. Percentiles are bounded
to [0.0001,0.9999]. The difference and average of the two normal scores are each
centered and divided by their weighted population standard deviation. These
are `Delta` (signed contrast) and `Level` (common activity component).

Higher contrast values indicate a shift toward a higher self-report rank relative
to device rank within the reference cycle. It does not demonstrate misreporting, absolute
measurement error or that either instrument is a criterion standard. After
centering, equal ranks need not correspond to contrast zero. The mean/difference
construction does not guarantee independence after selection or adjustment.

Full floating-point precision is retained when exporting and rereading current
cohorts. Raw-value ties, log-input ties and ties introduced by decimal serialization
must not be conflated. The four earlier processing sensitivities retained in S15
explicitly reproduce the historical twelve-significant-digit serialization;
current S12 analyses use full precision. No previous coefficients are loaded
as inputs. The old S7 descriptive CSV reader introduced two different raw-value
tie counts; manuscript correction is tracked in the release audit.

## Survey design and models

Examination weights are used because the analysis combines examination and
questionnaire information. Each compatible two-cycle period uses WTMEC2YR/2.
Cycle-qualified stratum and PSU identifiers keep cycles distinct. Domains are
formed through survey-design subsetting. Multiplying all weights in one
single-cycle fit by a common positive constant does not change the fitted
contrasts. The [NHANES weighting tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx)
provides the source guidance; the executable model specification is authoritative
for this study's precise domains and predictors.

Survey-weighted Cox models use Efron ties, examination as time origin, the
contrast, common activity component and core covariates. Pooled fits stratify
baseline hazards by cycle. Sex-specific slopes are obtained from the interaction
model using named linear contrasts and their full covariance, not by adding
standard errors. The preferred direct wrist-to-hip ratio of hazard ratios permits nuisance
slopes to differ by measurement period. It is a ratio of contrast hazard ratios,
not a participant mortality hazard ratio or a causal instrument effect.

Early-death exclusions remove specified deaths while retaining examination as
time origin. Landmark analyses instead require survival past the landmark and
subtract the landmark time. Administrative censoring changes event definitions
at the stated horizon while retaining the original sample. These are distinct
analyses. Matched BMI/disease comparisons hold sample support fixed.

Three-knot restricted cubic splines use weighted 10th/50th/90th percentile knots.
Plots compare with contrast zero and show pointwise intervals, not simultaneous
bands. Joint nonlinearity tests do not establish nonlinearity separately in both
sexes. Tail support and conditioning/influence diagnostics are exported.

Age regressions are survey-weighted linear models with linear and quadratic age,
the common activity component, sex and socioeconomic covariates, plus cycle when
pooled. The illustrative age-80 versus age-50 contrast is post hoc. Wrist public
age 80 denotes the top-coded 80+ group. These cross-sectional curves do not measure
individual aging. Separate sex-by-age curves are exploratory; interval overlap
is not a formal interaction test.

## Inference and sensitivity limits

The main scalar intervals and probabilities use full design degrees of freedom
(PSUs minus strata). Residual-model t and normal references are also exported
without changing coefficients or covariance. Joint tests use Wald/q with the
corresponding F reference, or the unscaled Wald statistic with a chi-square
reference. Low residual degrees of freedom can produce very wide intervals;
an unestimable reference is not silently replaced by a favorable alternative.

Response sensitivities refit survey logistic models for the 2005–2006 hip cycle
and pooled wrist period. Eligible complete questionnaire/predictor records
include responders and nonresponders. Predictors are age, sex, race/ethnicity,
education, poverty-income ratio, BMI and smoking, plus wrist cycle. Responders'
examination weights are multiplied by weighted response prevalence divided by
predicted response probability. Weights are not trimmed. Mortality outcome and
follow-up duration are not predictors, although linkage/follow-up eligibility
is part of the domain. Model-based response weighting cannot guarantee removal
of selection bias or recovery of unmeasured predictors.

Reported standard errors condition on generated scores and fitted response
weights. They do not propagate uncertainty through the entire scoring/response
pipeline. No full-pipeline bootstrap is claimed. Conventional weighted Schoenfeld
and cluster-robust exact log-time checks are distinguished from survey step-time
tests; none proves all proportional-hazards assumptions. Dfbeta is approximate
coefficient-scale influence, not exact survey deletion or standardized DFBETAS.

Retained exploratory analyses disclose the historical nonlinear screen and
separate BH families (four mortality modifiers, three nonlinear screens, eight
participant characteristics). BH adjustment does not cover every analysis in
the paper or undo data-informed model selection. Numerical parity, synthetic
tests and visual checks support reproducibility, not causal validity or immunity
to all alternative analytical choices.

## Post-hoc hip-index sensitivity without transport

A questionnaire/code audit identified that PAQ050Q counts transport occasions,
whereas PAD080 reports total minutes on active days. Their product is retained
in the original operational index but is not a validated duration estimate.
The index is therefore labelled in index units per day, not transport minutes.

The additional sensitivity removes transport, using
(home/yard + moderate leisure + 2 × vigorous leisure)/30 for 2003–2006.
Original reference domains, including transport completeness, are retained.
Self-report ranks and both Delta and Level are recomputed within each original
cycle reference domain with examination weights; device ranks are unchanged.
No participants are added. Mortality and age models retain their original
samples, weights, covariates, follow-up and inference conventions. Both original
and alternative models are refitted. This was specified after the audit finding
and before inspecting its results, not prospectively before study development.

Tables S30 and S31 report same-sample mortality and age contrasts. Each estimate
uses the SD and common activity component of its own score definition. Changes
are descriptive and must not be tested by treating paired model estimates as
independent. The hip-cycle age variable is top-coded at 85; age 80 is below that threshold. Similar estimates do
not validate the original transport-duration arithmetic.

## Post-hoc reporting additions

Tables S32–S33 provide nonsequential missing counts and overlapping criterion failures across the examined adult survey, paired-activity reference, principal mortality and historical age-analysis domains. Counts and percentages are unweighted; unconstructed scores outside eligible reference domains must not be interpreted as item nonresponse. The prepared model variables include principal covariates, health composites, score variants and applicable device alternatives. Raw questionnaire routing is documented separately in the source/scoring tables. Domain definitions and original scores are unchanged.

Table S34 sums PERMTH_EXM/12 in the principal cases without survey weighting; it reports arithmetic mean and sample SD, distinct from weighted Table 1. Table S35 uses the same cycle reference sample, weights, arithmetic and empirical inverse-CDF rule as prepare_models.py; ties go to the lower quintile. Full-precision cutpoints and exact membership checks are exported.

Table S36 compares a contrast-only model, a contrast-plus-average model and the original principal model on identical principal complete cases. All retain the survey design and pooled cycle-specific baseline hazards. Efron ties and the existing design-df t convention are unchanged. The unadjusted coefficient is not conditional on average activity; it targets a different association from the principal coefficient. Changes across adjustments do not establish confounding or mediation. No additional multiplicity correction or hypothesis selection was introduced.

## Post-hoc time-dependence diagnostics for S36

`R/run_reporting_diagnostics.R` evaluates all 18 S36 fits, without selecting models by diagnostic p values. Weighted Cox refits reproduce each coefficient within 1e-8. Kaplan–Meier-transformed Schoenfeld tests cover every term and the global model; they are conventional diagnostics, not full-survey-design tests. Contrast plots are exported locally.

A separate counting-process model splits follow-up at five years and adds only a contrast-by-late-time interaction to each original adjustment set. Deaths at exactly five years belong to the early interval. The horizon follows the existing diagnostic horizon, not the new results. The full survey base is retained before domain subsetting; excluded records receive arbitrary positive times solely to retain survey information and never enter model risk sets. Without the interaction, split fits must reproduce both coefficients and survey covariance within 1e-8, with identical events, person-time and design degrees of freedom. Early HRs, late HRs and late/early ratios use full design-df t inference; the late HR variance includes the coefficient covariance.

These focused step tests retain constant nuisance slopes and cannot exclude other forms of time variation. Conventional diagnostics flag contrast nonproportionality in the 2005–2006 and pooled hip unadjusted/average-component models. No five-year step test has p<0.05 (minimum 0.0616). This does not overturn the conventional signals. S36 remains a descriptive comparison of distinct associations, not a mediation analysis. All diagnostics are exploratory and unadjusted for multiplicity.

## Terminology and source documentation

Publication tables label exponentiated coefficient differences as ratios of hazard ratios (HR ratios). Legacy strings such as `Female/Male slope ratio` and `Wrist/Hip slope ratio` remain internal CSV lookup keys for reproducibility; they never denote a quotient of log-hazard coefficients. Exporters translate these keys to the publication labels without changing values. The NCHS PAXDAY_G and PAXDAY_H codebooks were first published in November 2020; 2011–2012 and 2013–2014 identify the survey cycles.

Rank transformations define relative scales; they do not by themselves establish measurement equivalence, improve every aspect of robustness, or add information beyond the paired activity measures. The measurements were paired at baseline but did not cover identical observation periods. Between-period comparisons combine changes in instruments, processing, population and follow-up, and cannot identify a single cause of the observed difference.

## Additional exploratory age-by-period comparison (S40)

The primary additional domain retains the historical period-specific age-analysis
eligibility and restricts ages to 20–79, avoiding different public-use upper-age
coding. Sensitivities (i) apply the same BMI and combined diabetes/CVD/cancer
completeness rule to both periods at ages 20–79 and (ii) apply this harmonized rule
with age capped at 80 in both periods. A positive disease component is sufficient
for positive composite status; a negative composite requires all components negative.
Original cycle-specific contrast and common-component scores are retained.

Gaussian survey regression uses (age−50)/10, its square, common activity component,
sex, race/ethnicity, education, PIR and cycle intercepts. Both age terms and all
nuisance coefficients may differ by period; cycle intercepts absorb the period
intercept difference. Two interaction coefficients are jointly tested. The design
is constructed before subsetting, using cycle-specific strata/PSUs and half the
two-year examination weights. The primary Wald F/t convention uses 62 survey
design degrees of freedom. Residual-df (39) and asymptotic alternatives are also
exported, not selected according to significance. Age-75 versus age-50 contrasts
are calculated within period and directly differenced (wrist minus hip), with
pointwise 95% intervals. The asymptotic joint statistic is chi-square=2×F, df=2.

The fully interacted model is independently checked against separate period fits
for both age coefficients and covariance (absolute tolerance 1e−10). The original
age-80 contrasts and sample sizes must reproduce before new tests are run. Tests
are post hoc and unadjusted for multiplicity. Cross-sectional period differences
cannot isolate device placement or estimate within-person aging. Original Figure 2
and Table 4 retain their full historical age domains and age-80 interpretation.
