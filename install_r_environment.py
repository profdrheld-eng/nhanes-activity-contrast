#!/usr/bin/env python3
"""Install hash-verified CRAN source archives into a new, isolated R library.

Download the files listed in config/r_sources.csv into --archives first.
Requires R and its platform's C/C++/Fortran build tools. No system libraries
are modified. Installation logs remain beside the new library.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archives', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    library = args.library.resolve()
    if library == root or root in library.parents:
        parser.error('Install library outside this repository')
    if library.exists():
        parser.error('Use a new, nonexistent library directory')
    executable = shutil.which('R')
    if executable is None:
        parser.error('R must be installed and available on PATH')
    with (root/'config/r_sources.csv').open() as stream:
        records = list(csv.DictReader(stream))
    for record in records:
        path = args.archives/f"{record['package']}_{record['version']}.tar.gz"
        if path.stat().st_size != int(record['bytes']) or hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']:
            raise ValueError(f'Archive verification failed: {path.name}')
    library.mkdir(parents=True)
    env = dict(os.environ, R_LIBS=str(library), R_LIBS_USER=str(library), R_LIBS_SITE='NULL')
    for record in records:
        path = (args.archives/f"{record['package']}_{record['version']}.tar.gz").resolve()
        print('Installing '+record['package']+' '+record['version'], flush=True)
        with (library/f"install_{record['package']}.log").open('w') as log:
            subprocess.run([executable, 'CMD', 'INSTALL', '--library='+str(library), str(path)],
                           env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    (library/'installation_sources.json').write_text(json.dumps(records, indent=2)+'\n')


if __name__ == '__main__':
    main()
