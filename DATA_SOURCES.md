# Public-use data sources

The package does not contain or redistribute participant-level data. Obtain
files directly from NHANES and NCHS. `config/sources.csv` lists the 34 files
required by the implemented cohort constructions, their official source URLs,
byte sizes, and SHA-256 checksums of the source snapshot used in this project.
All 34 source files were verified locally against these hashes. All download
links returned HTTP 200 with the expected file sizes during release preparation.
That link check does not replace verifying the downloaded content hashes.

Store each file in its cycle subdirectory, retaining the manifest filename:

```text
data/
  2003-2004/DEMO_C.XPT
  2003-2004/PAXRAW_C.zip
  ...
  2013-2014/NHANES_2013_2014_MORT_2019_PUBLIC.dat
```

The URLs identify NHANES questionnaire, examination, demographic, and monitor
files, and the NCHS 2019 public-use linked mortality release. Consult each
module's documentation rather than treating identical-looking variable names
as evidence of comparable measurement across cycles. For example:

- [NHANES 2003–2004 hip monitor documentation](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAXRAW_C.htm)
- [NHANES 2005–2006 individual leisure activity documentation](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAQIAF_D.htm)

A checksum mismatch stops processing. Do not replace the recorded checksum
merely to make a new download pass; determine whether the official release
changed, the download returned an HTML page, or a file was altered. No personal
account or project-local participant data should be necessary.

Raw files stay unchanged. Extracted monitor files and derived person-level
outputs belong in a separate local working directory, outside this repository.
Do not upload that working directory or saved model objects containing data.

## Reader precision

Minute-level hip data are read with ReadStat (`pyreadstat`). Some alternative
XPORT readers represent SAS numeric zero as a tiny positive floating-point
value; that is consequential for zero-run nonwear classification. Questionnaire
readers and their transformation parity are checked separately. Preserve
round-trip precision when writing/reading derived scores, as rounding raw
self-report inputs can create artificial rank ties.
