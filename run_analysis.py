#!/usr/bin/env python3
"""Public analysis entry point for source, cohort, and current model stages.

Raw files are supplied locally in cycle subfolders according to config/sources.csv.
No data are downloaded or uploaded by this command. Use a fresh work directory.
"""
import argparse
import json
from pathlib import Path
import platform
import sys
import shutil
import subprocess
from nhanes_activity.sources import source_records, verify_file, extract_xport, sha256
from nhanes_activity.environment import analysis_packages


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['all','verify-sources','process-hip','build-cohort','prepare-models','fit-models','export-displays','reporting'])
    parser.add_argument('--data-dir',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path)
    parser.add_argument('--cycle',choices=['2003-2004','2005-2006','2011-2012','2013-2014'])
    parser.add_argument('--chunksize',type=int,default=5_000_000)
    args = parser.parse_args()
    packages = analysis_packages() if args.stage != 'verify-sources' else None
    root = Path(__file__).resolve().parent
    if args.stage in ['all','prepare-models','fit-models','export-displays','reporting'] and args.cycle:
        parser.error('Model stages require all four cycles; omit --cycle')
    if args.stage=='process-hip' and args.cycle in ['2011-2012','2013-2014']:
        parser.error('process-hip requires a hip cycle')
    rows = source_records(root/'config/sources.csv')
    if args.cycle:
        rows = [r for r in rows if r['cycle']==args.cycle]
    if args.stage!='verify-sources' and args.work_dir is None:
        parser.error('--work-dir is required for processing')
    for row in rows:
        verify_file(args.data_dir/row['relative_path'],row)
    print(f'Verified {len(rows)} source files against SHA-256 and byte size',flush=True)
    if args.stage=='verify-sources':
        return
    from nhanes_activity.pax import process_xport
    # Keep generated person-level files outside the code repository, including
    # when callers provide a relative path or an existing directory symlink.
    work = args.work_dir.resolve()
    if work == root or root in work.parents:
        parser.error('Work directory must be outside the code repository')
    data = args.data_dir.resolve()
    if work == data or data in work.parents or work in data.parents:
        parser.error('Work and raw-data directories must be disjoint')
    if args.stage=='all':
        if work.exists():
            parser.error('The all stage requires a new, nonexistent work directory')
        for stage in ['process-hip','build-cohort','prepare-models','fit-models','export-displays','reporting']:
            subprocess.run([sys.executable,str(root/'run_analysis.py'),stage,
                            '--data-dir',str(data),'--work-dir',str(work),
                            '--chunksize',str(args.chunksize)],check=True)
        (work/'pipeline_complete.json').write_text(json.dumps({
            'status':'complete','stages':['process-hip','build-cohort','prepare-models','fit-models','export-displays','reporting']},indent=2)+'\n')
        return
    work.mkdir(parents=True,exist_ok=True)
    if args.stage=='prepare-models':
        subprocess.run([sys.executable,str(root/'prepare_models.py'),'--data-dir',str(data),
                        '--work-dir',str(work)],check=True)
        return
    if args.stage=='fit-models':
        rscript=shutil.which('Rscript')
        if rscript is None:
            parser.error('Rscript must be installed and available on PATH')
        stages=['response_weights.R','response_balance.R','run_models.R','check_regime_nuisance.R',
                'check_age_support.R','run_additional_sensitivities.R','run_retained_sensitivities.R',
                'age_sex_overlay.R','run_diagnostics.R','run_exploratory.R','run_no_transport.R']
        for stage in stages:
            print('Running '+stage,flush=True)
            subprocess.run([rscript,str(root/'R'/stage),str(work)],check=True)
        for exporter in ['export_retained.py','export_score_diagnostics.py']:
            subprocess.run([sys.executable,str(root/exporter),'--work-dir',str(work)],check=True)
        return
    if args.stage=='export-displays':
        rscript=shutil.which('Rscript')
        if rscript is None:
            parser.error('Rscript must be installed and available on PATH')
        for stage in ['export_table1.R','plot_density.R','plot_age_main.R','plot_flow.R','plot_supplement.R']:
            print('Running '+stage,flush=True)
            subprocess.run([rscript,str(root/'R'/stage),str(work)],check=True)
        for exporter in ['export_main_tables.py','export_supplement_tables.py','export_no_transport.py']:
            subprocess.run([sys.executable,str(root/exporter),'--work-dir',str(work)],check=True)
        return
    if args.stage=='reporting':
        rscript=shutil.which('Rscript')
        if rscript is None:
            parser.error('Rscript must be installed and available on PATH')
        subprocess.run([sys.executable,str(root/'prepare_reporting.py'),'--work-dir',str(work)],check=True)
        subprocess.run([rscript,str(root/'R/run_reporting.R'),str(work)],check=True)
        subprocess.run([sys.executable,str(root/'export_reporting.py'),'--work-dir',str(work)],check=True)
        subprocess.run([rscript,str(root/'R/run_reporting_diagnostics.R'),str(work)],check=True)
        subprocess.run([sys.executable,str(root/'export_reporting_diagnostics.py'),'--work-dir',str(work)],check=True)
        return
    if args.stage=='build-cohort':
        from nhanes_activity.cohorts import CYCLES, write_cycle
        for cycle in ([args.cycle] if args.cycle else CYCLES):
            print(json.dumps(write_cycle(args.data_dir,work,cycle)),flush=True)
        return
    cycles = [args.cycle] if args.cycle else ['2003-2004','2005-2006']
    for cycle in cycles:
        destination = work/cycle
        destination.mkdir(exist_ok=False)
        row = next(r for r in rows if r['cycle']==cycle and r['filename'].endswith('.zip'))
        xport = destination/'minutes.xpt'
        extract_xport(args.data_dir/row['relative_path'],xport)
        environment={'python':sys.version,'platform':platform.platform(),
                     'packages':{p:r['loaded_version'] for p,r in packages.items()},
                     'package_details':packages,
                     'input_sha256':row['sha256'],'xport_sha256':sha256(xport),
                     'chunksize':args.chunksize}
        (destination/'environment.json').write_text(json.dumps(environment,indent=2)+'\n')
        print(json.dumps(process_xport(xport,cycle,destination/'pax',args.chunksize)),flush=True)

if __name__=='__main__':main()
