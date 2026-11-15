from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import norm
from math import ceil

ROOT=Path('/home/ubuntu/revision_work')
OUT=ROOT/'original_preserving_revision'/'analysis'
OUT.mkdir(parents=True,exist_ok=True)
rows=[]
for year in (2005,2009):
    df=pd.read_csv(ROOT/'input'/f'Book{year}.csv')
    xy=df[['Long','Lat']].to_numpy()
    _,idx=cKDTree(xy).query(xy,k=9)
    neighbors=idx[:,1:]
    y=df.PresenceAb.to_numpy(dtype=float)
    z=y-y.mean()
    denom=(z*z).sum()
    obs=(z*z[neighbors].mean(axis=1)).sum()/denom
    rng=np.random.default_rng(20261104)
    extremes=0
    for _ in range(9999):
        q=rng.permutation(z)
        stat=(q*q[neighbors].mean(axis=1)).sum()/denom
        extremes+=abs(stat)>=abs(obs)
    rows.append({'year':year,'N':len(y),'k':8,'weights':'directed 8-nearest-neighbor row-standardized','Moran_I':obs,'two_sided_permutation_p':(extremes+1)/10000,'permutations':9999,'seed':20261104})
pd.DataFrame(rows).to_csv(OUT/'reference_label_moran.csv',index=False)
models=[(2005,'RF',251,4,7,37),(2005,'SVM',248,7,6,38),(2009,'RF',253,2,4,40),(2009,'SVM',253,2,5,39)]
metrics=[]
for year,model,tn,fp,fn,tp in models:
    n=tn+fp+fn+tp
    po=(tp+tn)/n
    pe=((tp+fp)*(tp+fn)+(tn+fn)*(tn+fp))/n**2
    mcc=(tp*tn-fp*fn)/np.sqrt(float(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    metrics.append({'year':year,'model':model,'N':n,'TN':tn,'FP':fp,'FN':fn,'TP':tp,'OA':po,'Kappa':(po-pe)/(1-pe),'MCC':mcc,'Precision':tp/(tp+fp),'Recall':tp/(tp+fn),'Specificity':tn/(tn+fp),'F1':2*tp/(2*tp+fp+fn)})
pd.DataFrame(metrics).to_csv(OUT/'corrected_original_holdout_metrics.csv',index=False)
p=162/1021
print('Cochran n before rounding:',norm.ppf(.975)**2*p*(1-p)/.05**2,'ceiling:',ceil(norm.ppf(.975)**2*p*(1-p)/.05**2))
print(pd.DataFrame(rows).to_string(index=False))
print(pd.DataFrame(metrics).to_string(index=False))
