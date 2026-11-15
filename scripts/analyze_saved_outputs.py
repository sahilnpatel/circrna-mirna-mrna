"""Modern descriptive reanalysis of archived outputs, not a full study replication."""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr, spearmanr, ttest_ind
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from data_io import aligned_rows,write_csv
from plots import scatter


def bh_adjust(p):
    p=np.asarray(p); order=np.argsort(p); ranked=p[order]*len(p)/np.arange(1,len(p)+1)
    out=np.empty(len(p));out[order]=np.minimum.accumulate(ranked[::-1])[::-1].clip(0,1)
    return out


def main():
    rows=aligned_rows(ROOT/'data/source');out=ROOT/'results';out.mkdir(exist_ok=True)
    stats=[]
    for cohort,rr in [('all_samples',rows),('primary_tumors',[r for r in rows if r['sample_type']=='01'])]:
        x=np.array([r['mir885_activity_proxy'] for r in rr]); y=np.array([r['MCM5'] for r in rr])
        pearson=pearsonr(x,y);spearman=spearmanr(x,y)
        stats.append({'cohort':cohort,'n':len(rr),'pearson_r':float(pearson.statistic),'pearson_p':float(pearson.pvalue),
                      'spearman_r':float(spearman.statistic),'spearman_p':float(spearman.pvalue)})
        scatter(out/f'{cohort}_correlation.svg',x,y,[{'01':0,'11':1,'02':2}[r['sample_type']] for r in rr],
                f'{cohort}: r = {pearson.statistic:.3f}','Inferred miR-885-3p activity (proxy)','MCM5 expression (saved scale)')
    write_csv(out/'correlations.csv',stats)
    rr=[r for r in rows if r['sample_type']=='01'];mcm=np.array([r['MCM5'] for r in rr]);cut=float(np.median(mcm))
    tests=[]
    for name in rr[0]:
        if not name.startswith(('GOBP_','REACTOME_')):continue
        val=np.array([r[name] for r in rr]);high=val[mcm>cut];low=val[mcm<=cut]
        test=ttest_ind(high,low,equal_var=False)
        tests.append({'pathway':name,'high_n':len(high),'low_n':len(low),'high_mean':float(high.mean()),
                      'low_mean':float(low.mean()),'difference':float(high.mean()-low.mean()),'p':float(test.pvalue)})
    for t,q in zip(tests,bh_adjust([t['p'] for t in tests])):t['bh_q']=float(q)
    write_csv(out/'pathway_comparisons.csv',tests)
    (out/'reanalysis.json').write_text(json.dumps({'median_mcm5_cutoff':cut,'cohort':'primary tumors',
        'test':'Welch t-test; Benjamini-Hochberg correction over 11 available saved pathways',
        'tie_rule':'MCM5 > median is high; ties enter low',
        'interpretation':'Association only; signatures can overlap and include MCM5. Not an exact reproduction of the paper.'},indent=2))
    print('Wrote correlations, descriptive pathway comparisons, and SVG figures.')


if __name__=='__main__':main()
