import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

df=pd.read_csv('./cyp71.txt',sep='\t')

df['pvalue_cyp71-2']=df.apply(lambda x: round(ttest_ind([x[2],x[3],x[4]],[x[5],x[6],x[7]],equal_var=False)[1],2),axis=1)
df['pvalue_cyp71-3']=df.apply(lambda x: round(ttest_ind([x[2],x[3],x[4]],[x[8],x[9],x[10]],equal_var=False)[1],2),axis=1)

df=df.fillna(1)

df.to_csv('./cyp71.txt',sep='\t',index=None)
