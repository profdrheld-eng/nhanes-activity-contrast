# Validation scope and evidence

## Version 1.0.0 release check (20 September 2026)

The complete `all` pipeline was rerun from the 34 manifest-verified public raw
files in a new work directory, using the documented clean Python 3.9.6 and
isolated R 4.6.0 environments. All seven stages completed: hip-minute processing,
cohort construction, model preparation, model fitting, displays, reporting and
the additional age comparison. This supersedes the earlier staged-replay limits
for the current code version; the dated records below remain as audit history.

The newly generated outputs matched all 45 computational table blocks
in the working manuscript and supplement, including both Table S40 panels.
S38/S39 are editorial reporting checklists and were checked separately rather
than presented as generated statistical results. Forty-four selected original
main-text numerical assertions and the added age-comparison result paragraph
also matched the new exports. Regressions against the preceding validated stage
outputs passed for 59 CSV files (567 column comparisons;
continuous values: absolute and relative tolerance 1e-8; integer-valued and
categorical values: exact). All six figure PNGs were byte-identical to the
previously reviewed corrected figures. The analysis source hashes were unchanged
throughout this final run. All 55 Python synthetic tests and seven R test scripts
passed. These are internal computational checks, not external validation.

The release package uses an explicit allowlist of 76 code, configuration,
documentation and synthetic-test files. It excludes participant-level files,
fitted objects, manuscripts, review comments and local audit paths. A review of
the existing six-commit history found no participant-data files, private paths
or matches to the tested credential patterns; this is not a formal security
audit. Historical documentation remains visible in Git history.

Original code and repository documentation are released under MIT with software
citation metadata in CITATION.cff. Version `v1.0.0` identifies the code used for
the paper. Later changes require a new version. The earlier private-hosting and
undecided-license notes below describe the development stage, not this release.

Historical release preparation (18 September 2026). The second critical review corrected a density-coordinate error and added the author-approved post-hoc transport sensitivity. Historical clean-Python raw replay and model/export integration passed, but reproduction alone did not detect the Figure 1 coordinate error described below. All reconstruction and cross-checks described below were performed within the study project; they do not constitute external validation. This is an internal reproducibility audit,
not external peer review or a guarantee of statistical correctness.

## Input and processing checks

- The manifest covers 34 public NHANES/NCHS inputs, totaling 1,005,324,383 bytes.
  Local files match every recorded size and SHA-256. Official links were
  checked for HTTP 200 and matching size; remote HEAD checks do not substitute
  for validating the bytes a reproducer downloads.
- Raw hip reconstruction processed 72,250,027 minute rows for 7,176 participants
  in 2003–2004 and 74,874,095 rows for 7,455 participants in 2005–2006. No participant
  required structural quarantine in these input snapshots.
- Four reconstructed cohorts were compared against retained study cohorts,
  including 396 complete-column checks. Membership, categorical codes and
  count-valued fields are checked exactly. Historical floating-point cohort
  serialization requires explicit numerical tolerances.
- Primary mortality samples by cycle are 2,732/680,2,737/505,3,860/400 and
  4,018/314 participants/deaths. The pooled totals are 5,469/1,185 and 7,878/714.
  Paired-reference samples are 2,871,2,844,4,185 and 4,329. These are different
  domains, not interchangeable denominators.

## Environment and numerical reproduction

The tested Python environment is Python 3.9.6 with all seven versions fixed in
requirements.txt. The tested R environment is R 4.6.0 with the ten exact sources
in config/r_sources.csv. All ten R packages were rebuilt in a new isolated
library from verified CRAN archives, without using the prior project library.
Models were then refitted from copied, independently reconstructed cohorts.
The separate fresh-Python raw-to-display replay completed successfully and
reproduced 102 comparable output files exactly (gzip content compared after
decompression; elapsed-time and environment metadata excluded). After adding
the response-balance exporter to the common entry point, model preparation,
all model stages and all displays were run together again using the fresh
Python and R installations and the freshly reconstructed cohorts. All 102
comparable outputs matched the preceding clean-R validation. The final
diagnostic addition did not change raw processing, weights or fitted models.

The clean R model replay reproduced 390 columns in 40 result CSVs. Continuous
values used prespecified absolute and relative tolerances of 1e-8; integral
and categorical fields were required to match exactly. These tolerances are
stricter than the displayed manuscript rounding. The paired model inputs and
response weights had also been checked against retained study outputs.

The original 35 tables matched accepted-view cells in the manuscript/supplement,
including labels, counts and rounded numerical values. Six freshly rendered
figures and both table JSONs were byte-identical across the two R installations.
The six figures were visually inspected. They are not claimed pixel-identical
to historical Word graphics. Forty-one selected main-text numerical assertions
also match fresh exports; that assertion check supplements the table audit.

## Second critical review

The original density exporter used `tapply` without explicit bin levels. Empty
coordinate bins were dropped, so weighted mass was assigned to incorrect
plotted coordinates. Historical byte identity did not establish correctness.
The corrected exporter retains every bin on both axes. An independent
synthetic grid with gaps, unequal weights and boundary observations verifies
placement and mass conservation. The old Figure 1 must not be used.

Cox models can return finite coefficients and covariance after warning of an
infinite coefficient. All 16 production model-fit call sites now reject
non-estimability/convergence warnings and non-finite results. A separated
synthetic survey Cox model verifies this failure path; an ordinary fit
retains its exact coefficient and covariance. Unrelated warnings remain visible.

Output guards now protect every named output in multi-output R stages,
including partial earlier exports. Synthetic subprocess checks verify that
each density/flow output is protected before loading input cohorts.

After these changes, all model and display stages were replayed with the
clean Python and R environments from independently reconstructed cohorts.
Of 95 comparable CSV, JSON and PNG files, 93 were byte-identical. Only
the Figure 1 raster CSV and PNG changed, as expected. All fitted-model
CSV outputs, all 35 table exports and the other five figures were unchanged.
This replay did not repeat raw-minute processing because that code was unchanged.

Named contrasts additionally reject duplicate, unnamed, empty, unknown and
non-finite requests. Covariance axes are aligned to coefficient names rather
than assumed to share their ordering. A duplicate-term regression test failed
before the correction and passes afterward; a permuted-covariance test confirms
unchanged valid inference. A second complete model/display replay reproduces
all 95 compared CSV, JSON and PNG files exactly after this guard was added.

## Internal cross-checks and adversarial tests

36 Python tests cover weighted tied ranks, normalization and invalid inputs;
questionnaire gates, frequencies and missingness; nonwear/day boundaries,
participant chunk boundaries and malformed keys; source hashes and extraction;
stale package metadata; and complete, unique, finite sensitivity exports. They use synthetic inputs, not study participants.

`Rscript tests/test_inference.R` checks named linear contrasts using specified
coefficients, a non-diagonal covariance matrix and a tabulated t critical value.
It checks exponentiation, confidence intervals, design degrees of freedom and
a joint Wald calculation. `Rscript tests/test_balance.R` independently checks
target-population/Bernoulli standardization and zero-variance behavior.

`Rscript tests/test_density.R`, `tests/test_fit_checks.R` and
`tests/test_output_guards.R` test the three failure modes above. The sixth R script, `tests/test_no_transport.R`, checks same-reference reranking,
weighted ties, unchanged device scores and transport independence. All six
R test scripts and all 36 Python tests pass in the recorded environments.

The wrist response-balance export reproduces all eight columns of the retained
balance table within 1e-10. The target comprises 9,844 eligible records and 7,814
responders; maximum absolute standardized differences are 0.060310 before and
0.001599 after weighting. This validates the reported descriptive calculation,
not balance of unmeasured characteristics.

## Post-hoc transport sensitivity and final integration

The questionnaire audit found that transport occasions and duration on active
days have different reference units. The original hip score is therefore an
operational activity index, not a validated daily-duration estimate. The author
approved omission of transport, keeping the original reference and analysis
samples and reconstructing both the contrast and common component. The
specification was recorded before inspecting the new estimates.

An independent Python calculation using standard-library normal quantiles
checked every reconstructed score against the R implementation. Maximum
absolute differences were below 1.7e-12 for standardized components. Original
columns, device ranks, participant IDs and reference/analysis masks were retained.
All 45 original-model refit rows matched their previous numerical results.

The pooled mortality HR was 1.224 (95% CI 1.107–1.354), compared with 1.235
(1.127–1.352) originally, in the same 5,469 participants with 1,185 deaths.
The female estimate remained 1.130 but its interval became 0.996–1.281.
The pooled age contrast was 1.136 SD (1.018–1.254), versus 1.139 SD
(1.017–1.260), in the same 5,415 participants. These descriptive comparisons
do not validate the original score or establish differences between estimates.

A fresh model/export integration from the verified reconstructed cohorts
compared 111 CSV, JSON, gzip-CSV and PNG outputs with the preceding run:
106 were identical, including the new sensitivity fits and all six figures.
The five changed presentation files contain only the intended hip self-report
label corrections. The strengthened exporter reproduced S30/S31 exactly and
rejects missing, duplicated, unexpected, non-finite or inconsistent model rows.
All 37 generated tables match the revised DOCX candidate cells exactly.
Word revision authorship, retained comments/revisions and open reply threads
were checked separately; numerical parity is not a native-Word layout test.

## Differences investigated instead of concealed

- Two historical S7 counts reflected a CSV-reader precision difference. The
  current raw-value export preserves full precision; the working supplement
  was corrected to 891 unique values and 2,264 persons in tied groups for the
  2005–2006 self-report raw values. Actual log-input ties are exported separately.
- Four retained S15 processing sensitivities explicitly reconstruct historical
  12-significant-digit serialization before reranking. Current S12 uses full
  precision. No previous fitted coefficient is supplied as an analysis input.
- An older local environment contained conflicting pyreadstat distribution
  metadata. The imported code was 1.2.8; the clean environment and new preflight
  remove this ambiguity. Older audit records are not retrospectively rewritten.

## Interpretation and publication boundary

METHODS.md states limitations involving selection, measurement regimes,
generated-score/response-weight uncertainty, exploratory selection and model
assumptions. Reproducing an analysis does not eliminate those limitations.

Only reviewed code, configuration, documentation and synthetic tests belong
in the release archive. Local inputs, individual-level outputs, fitted RDS
objects, manuscripts, comment threads and caches are excluded. No Git history
is included. Private hosting under `profdrheld-eng/nhanes-activity-contrast` has been authorized. License selection and public visibility remain separate author decisions; private hosting does not establish public code availability.

## Reporting extension on 18 September 2026

The extension reuses a byte-identical copy of the previously validated revision_cohort.csv.gz and mortality_results.csv. Five new synthetic Python checks cover nonsequential counts and zero values, weighted empirical quantiles and boundary ties, invalid weights, unweighted person-years, and disease composites with partially unknown components. A fresh run in Python 3.9.6 and the recorded clean R 4.6.0 library exported 560 completeness cells, 112 criterion-failure cells, four full-precision cutpoint sets and six follow-up summaries. Every original quintile membership was reproduced.

Eighteen survey Cox fits were estimated on unchanged principal complete cases. Each of the six principal refits matched its original coefficient, SE, HR, interval and p value within 1e-8 for all three retained inference conventions; sample size, deaths and design degrees of freedom agreed exactly. This is a reporting-stage replay, not another raw-minute reconstruction or independent external validation. Fitted objects remain local; the code release contains no participant data.

## S36 time-dependence diagnostics (19 September 2026)

All 18 weighted Cox refits reproduced S36 coefficients within 1e-8. Counting-process refits without time interactions also reproduced coefficients and full survey covariance within 1e-8, preserving original membership, deaths, person-time and design degrees of freedom. Focused five-year interaction fits all converged with finite estimates/covariance. The synthetic R regression check covers deaths exactly at five years, late deaths, censoring, person-time conservation and excluded records with missing follow-up. Full raw-minute processing was not rerun for this diagnostic addition. Diagnostic findings and inferential limitations are described in METHODS.md.

## Portable S37 table export

On 2026-09-19, `export_reporting_diagnostics.py` was added after the diagnostic R stage. It formats all 18 models without refitting them and requires the diagnostic completion marker, complete model/contrast coverage, unique keys, finite valid probabilities and HR intervals, and matching sample support. Existing exports are protected against overwrite. Its CSV and JSON match all 95 cells of the current supplementary table, including the 72 result cells. The accompanying notes distinguish conventional Schoenfeld tests from survey-based step-time inference.

All 45 Python unit tests passed in the documented Python 3.9.6 project environment, including four new S37 tests. A supplementary run under Python 3.12 passed the S37 tests but exposed two existing mock-resolution errors in `test_environment.py`; this does not establish Python 3.12 compatibility. The R diagnostic calculations and fitted models were unchanged and were not rerun for this formatting addition.

## Source and presentation revision (19 September 2026)

Publication labels now use “HR ratio” for exponentiated coefficient differences; legacy model-output lookup keys are retained. PAXDAY documentation references now use their 2020 publication year. The updated main, supplement and no-transport exporters were executed against the retained validated model results. Cell-by-cell comparisons found 115 changed cells, exclusively the approved labels and one codebook-reference cell. No numerical result changed. All 45 tables in the revised manuscript and supplement matched their corresponding generated or retained reporting/checklist tables. This check did not refit models or repeat raw-minute reconstruction.

The accompanying methods documentation distinguishes relative ranks from measurement equivalence and paired baseline assessments from identical observation periods. These changes clarify interpretation; they do not validate the operational activity score or resolve the previously reported study limitations.

## Additional age-comparison integration (20 September 2026)

The new `age-comparison` CLI stage was executed with the recorded Python 3.9.6
and isolated R library, after verifying all 34 source hashes. It used a byte-copy
of previously validated model inputs; minute-level raw-data processing was not
repeated for this additive stage. No dependency or previous analysis changed.

All seven result CSVs reproduce the pre-integration additional analysis with zero
numerical difference: baseline contrasts (2 rows), coefficients (72), age contrasts
(27), curves (540), joint tests (9), domain support (6), and independent coefficient/
covariance checks (3). Both S40 panels match all cells in the current supplement.
The original period-specific age-80 estimates and sample sizes are checked before
fitting. Independent period fits reproduce the combined interaction coefficients
and covariance within 1e−10. The two-df test is also checked against `regTermTest`.

The Python suite passes 55 synthetic tests, including ten new S40 export tests.
Seven existing R test scripts pass. A repeated R-stage invocation fails explicitly
without changing existing output hashes. A separate display replay reproduces
both table panels and includes the explanatory note. These checks establish
computational agreement, not validation of every scientific assumption. The new
stage is wired into `all`; the entire minute-level `all` pipeline was not rerun
for this extension. No participant-level files or fitted objects are included.
