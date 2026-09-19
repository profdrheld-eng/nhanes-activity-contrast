# NHANES activity-contrast analyses

Code-only release candidate for the current manuscript's 2003–2006 hip and
2011–2014 wrist analyses. The second critical review and the authorized post-hoc transport sensitivity are documented in VALIDATION.md. Prepared for the private repository `profdrheld-eng/nhanes-activity-contrast`; public release is pending.
Participant-level data must be obtained separately from official NHANES/NCHS
sources and are not included in this package.

## What is reproduced

The pipeline reconstructs four cycle cohorts from 34 public input files,
fits the main and sensitivity models, refits response weights, exports
response-balance and other diagnostics, and generates 43 tables and 6 figures, including the post-hoc reporting additions (S32–S37). Table S37 is assembled automatically from the diagnostic CSV outputs.
It does not read historical coefficients or the Word manuscripts as inputs.

- [DATA_SOURCES.md](DATA_SOURCES.md): exact official inputs and folder layout.
- [METHODS.md](METHODS.md): operational definitions, domains, models and limitations.
- [RESULTS_MAP.md](RESULTS_MAP.md): every table, figure and central result's producer.
- [VALIDATION.md](VALIDATION.md): evidence, tolerances, differences and test scope.

## Install the recorded environment

The validated environment uses Python 3.9.6 and R 4.6.0. Other complete runtime
combinations have not been validated. First ensure `python3` refers to the
intended Python version, then create an isolated environment outside this folder:

```sh
python3 --version
python3 -m venv /your/local/nhanes-env
. /your/local/nhanes-env/bin/activate
python -m pip install -r requirements.txt
```

All seven Python dependencies, including transitive dependencies, are pinned.
The entry point checks that imported package versions agree with their
installed distribution metadata. Conflicting metadata produces an error.

`config/r_sources.csv` fixes the ten non-base R packages to exact versions,
official CRAN source URLs, file sizes and SHA-256 hashes. Download those named
archives into a separate local folder, then run:

```sh
python install_r_environment.py --archives /your/local/r-archives --library /your/local/nhanes-r-library
export R_LIBS=/your/local/nhanes-r-library
export R_LIBS_USER=/your/local/nhanes-r-library
export R_LIBS_SITE=NULL
```

The library directory must not already exist. The installer validates every
archive before installing anything and never changes a system library. Source
installation requires the C/C++ and Fortran toolchain appropriate to R 4.6.0.
Build logs remain in the new library. All ten exact source packages were
successfully rebuilt and used for the validation described in VALIDATION.md.
The shell commands shown are for macOS/Linux; Windows setup is not tested.

## Test and run

From this code directory, with the environments above active:

```sh
python -m unittest discover -s tests -v
Rscript tests/test_inference.R
Rscript tests/test_balance.R
Rscript tests/test_fit_checks.R
Rscript tests/test_density.R
Rscript tests/test_output_guards.R
Rscript tests/test_no_transport.R
python run_analysis.py verify-sources --data-dir /your/local/data
python run_analysis.py all --data-dir /your/local/data --work-dir /your/local/new-work --chunksize 2000000
```

Only synthetic data are used by unit tests. The data directory must follow
DATA_SOURCES.md and every file must match its manifest hash and size.
`all` requires a **new, nonexistent** work directory, separate from both the
raw-data directory and this code directory. It does not download or upload
study data. Original input files are not overwritten.

The explicit execution order is:

1. `process-hip`: extract and process both minute-level hip files, all five variants.
2. `build-cohort`: reconstruct all four cycles, eligibility and activity scores.
3. `prepare-models`: prepare common model inputs and descriptive summaries.
4. `fit-models`: refit response weights, balance, main/secondary/exploratory models and diagnostics.
5. `export-displays`: generate the original 37 tables and figures.
6. `reporting`: export domain-specific completeness, follow-up, cutpoints and same-case adjustment models, then generate Tables S32–S36.

For an intentional stage-by-stage run, replace `all` with the relevant stage
and reuse that run's work directory in the order above. Do not repeat stages
whose outputs already exist or mix partial outputs from different attempts.
`--cycle` is available only for cycle-specific source/cohort processing;
model stages require all four cycles. A successful `all` command writes
`pipeline_complete.json`. A failed or interrupted run does not.

## Outputs and resources

Numerical models and descriptive exports are in `work/models/`, response
fits and balance in `work/response/`, exploratory results in `work/exploratory/`,
and final table/figure presentations in `work/displays/`. The latter includes
`main_tables.json`, `supplement_tables.json`, `no_transport_tables.json`,
`reporting_tables.json` and six PNGs. New reporting CSVs and fitted objects are in
`work/reporting/`; never publish the fitted RDS objects.

The 34 source files total about 1.01 GB. A completed reference working directory
occupies about 5.3 GB, including extracted minute files and intermediates.
Allow additional space for installed environments and separate validation runs.
Approximately 147 million hip minute rows took about 25 minutes to process on
the validation machine at `--chunksize 2000000`. The buffer size can be changed,
but the command above reproduces the tested setting. A separately measured R
diagnostic stage used about 2.0 GiB peak resident memory and 11.6 seconds.
These are machine- and stage-specific observations, not a whole-pipeline
peak-memory measurement or guaranteed minimum system requirements.

Never upload the external work directory: CSV/CSV.GZ files and fitted RDS
objects can contain individual records. The release contains only reviewed
code, configuration, documentation and synthetic tests. It excludes study data,
model objects, manuscript comments, local caches and Git history.

## AI assistance and author responsibility

GPT-6 Astra (OpenAI), a large language model accessed through Codex, assisted with code development, debugging, formatting, explanatory code comments, and reproducibility checks, as well as the preparation and language editing of README files, execution instructions, and other repository documentation. Automated tests and comparisons against the reported results support computational reproducibility; they do not establish the validity of all methodological assumptions. The authors retain responsibility for the analytical decisions and the released code.

## Interpretation and publication

Numerical reproduction does not prove model assumptions or causal validity.
METHODS.md documents measurement differences, selection, exploratory choices,
and uncertainty not propagated through the complete scoring/weighting pipeline.
Historical serialization choices are explicitly retained where required to
reproduce earlier reported sensitivities; current analyses retain full precision.

The repository is prepared for private hosting under `profdrheld-eng/nhanes-activity-contrast`. Public visibility requires a separate author decision. The code license has not yet been selected; no open-source license is granted by this package. Before public release, select a license, review the complete tracked history, and record the exact paper-specific commit or release in the manuscript.

STROBE and RECORD checklists (Tables S38–S39) belong to the journal Supplementary File and are not separate code outputs or repository attachments.

The reporting stage can be added to a completed pre-reporting run using
`python run_analysis.py reporting --data-dir /your/local/data --work-dir /your/local/completed-work`.
It requires the existing model inputs and mortality results, refuses existing reporting
outputs, and checks that its principal refits match those results. The added stage
was tested separately on a byte-identical copy of the validated model input; the
minute-level pipeline was not rerun for this reporting-only extension.

The `reporting` stage also runs `R/run_reporting_diagnostics.R` after the S32–S36 exports. It writes Table S37 source values to `reporting_diagnostics/schoenfeld.csv`, `step_time.csv` and `support.csv`, together with 18 local diagnostic plots. To add diagnostics to an existing completed reporting run, use `Rscript R/run_reporting_diagnostics.R /path/to/work`. The output directory must not already exist; a completion marker is written only after all 18 models pass numerical checks. After successful diagnostics, `export_reporting_diagnostics.py` automatically writes `displays/table_s37.csv`, `reporting_diagnostic_tables.json` and `reporting_diagnostic_table_notes.json`. For an existing completed diagnostic run, use `python export_reporting_diagnostics.py --work-dir /path/to/work`; existing S37 exports are never overwritten. No Word file or local audit folder is required. These diagnostics do not replace the original models or establish proportionality.

Run the boundary/event/person-time regression check with `Rscript tests/test_reporting_time_split.R` in the same R environment.
