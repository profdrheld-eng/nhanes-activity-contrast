#!/usr/bin/env python3
"""Export tie counts and the centered value of equal activity ranks."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--work-dir',type=Path,required=True)
OUT=parser.parse_args().work_dir.resolve()/'models'
if any((OUT/name).exists() for name in ('rank_ties_and_bounds.csv','transform_input_ties_and_bounds.csv','contrast_zero_reference.csv')):
    raise FileExistsError('Refusing to overwrite score diagnostics')
CY=['2003-2004','2005-2006','2011-2012','2013-2014']
dist=pd.read_csv(OUT/'contrast_distribution.csv',float_precision='round_trip')
# Transparent finite-tail/tie checks and the location of equal ranks after centering.
cohort=pd.read_csv(OUT/'revision_cohort.csv.gz',low_memory=False,float_precision='round_trip')
tie=[]
for cy in CY:
    x=cohort[(cohort.cycle==cy)&cohort.reference_domain.astype(str).str.lower().isin(['true','1'])]
    w=x.WTMEC2YR.to_numpy()
    for name,v in [('SR_day',x.SR_day.to_numpy()),('device_raw',x.device_raw.to_numpy()),
                   ('log1p_SR_day',np.log1p(x.SR_day.to_numpy())),('device_log',x.device_log.to_numpy())]:
        o=np.argsort(v,kind='stable'); vs=v[o];ws=w[o]
        starts=np.r_[0,np.flatnonzero(np.diff(vs)!=0)+1];ends=np.r_[starts[1:],len(v)]
        ranks=np.empty(len(v));cum=0
        for s,e in zip(starts,ends):
            tw=ws[s:e].sum();ranks[s:e]=(cum+tw/2)/w.sum();cum+=tw
        tie.append(dict(cycle=cy,measure=name,n=len(x),unique_values=len(starts),
                        persons_in_tied_groups=int(sum(e-s for s,e in zip(starts,ends) if e-s>1)),
                        clipped_low=int(sum(ranks<.0001)),clipped_high=int(sum(ranks>.9999))))
ties=pd.DataFrame(tie)
ties[ties.measure.isin(['SR_day','device_raw'])].to_csv(OUT/'rank_ties_and_bounds.csv',index=False)
# Monotone transforms can collapse distinct floating-point inputs into ties.
# Preserve the manuscript's raw-value table and disclose actual rank inputs separately.
ties[ties.measure.isin(['log1p_SR_day','device_log'])].to_csv(OUT/'transform_input_ties_and_bounds.csv',index=False)
dd=dist[dist.weighting=='weighted'].copy();dd['equal_rank_location_standardized']=-dd['mean']/dd['sd']
dd[['cycle','mean','sd','equal_rank_location_standardized']].to_csv(OUT/'contrast_zero_reference.csv',index=False)
