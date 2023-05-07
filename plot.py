#volcano plot

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df=pd.read_csv('./cyp71.txt',sep='\t')
df['cyp71-2/Col-0']=df['cyp71-2']/df['Col-0']
df['cyp71-3/Col-0']=df['cyp71-3']/df['Col-0']

##############cyp71-2######################
c2part1=df[(df['pvalue_cyp71-2']<=0.05)&(df['cyp71-2/Col-0']>=3/2)]#upregulated, red dot
c2part2=df[(df['pvalue_cyp71-2']<=0.05)&(df['cyp71-2/Col-0']<=2/3)]#downregulated, blue dot
c2part3=df[(df['pvalue_cyp71-2']>0.05)|((df['cyp71-2/Col-0']<3/2)&(df['cyp71-2/Col-0']>2/3))]#grey dot
plt.axhline(-np.log10(0.05),color='grey',linestyle='--')
plt.axvline(np.log2(3/2),color='grey',linestyle='--')
plt.axvline(np.log2(2/3),color='grey',linestyle='--')
plt.scatter(np.log2(c2part1['cyp71-2/Col-0']),-np.log10(c2part1['pvalue_cyp71-2']),color='r')
plt.scatter(np.log2(c2part2['cyp71-2/Col-0']),-np.log10(c2part2['pvalue_cyp71-2']),color='b')
plt.scatter(np.log2(c2part3['cyp71-2/Col-0']),-np.log10(c2part3['pvalue_cyp71-2']),color='grey')
plt.xlim((-4,4))
plt.ylim((0,3))
plt.savefig('volcano_cyp71-2.pdf')
plt.show()


##############cyp71-3######################
c3part1=df[(df['pvalue_cyp71-3']<=0.05)&(df['cyp71-3/Col-0']>=3/2)]#upregulated, red dot
c3part2=df[(df['pvalue_cyp71-3']<=0.05)&(df['cyp71-3/Col-0']<=2/3)]#downregulated, blue dot
c3part3=df[(df['pvalue_cyp71-3']>0.05)|((df['cyp71-3/Col-0']<3/2)&(df['cyp71-3/Col-0']>2/3))]#grey dot
plt.axhline(-np.log10(0.05),color='grey',linestyle='--')
plt.axvline(np.log2(3/2),color='grey',linestyle='--')
plt.axvline(np.log2(2/3),color='grey',linestyle='--')
plt.scatter(np.log2(c3part1['cyp71-3/Col-0']),-np.log10(c3part1['pvalue_cyp71-3']),color='r')
plt.scatter(np.log2(c3part2['cyp71-3/Col-0']),-np.log10(c3part2['pvalue_cyp71-3']),color='b')
plt.scatter(np.log2(c3part3['cyp71-3/Col-0']),-np.log10(c3part3['pvalue_cyp71-3']),color='grey')
plt.xlim((-4,4))
plt.ylim((0,3))
plt.savefig('volcano_cyp71-3.pdf')
plt.show()
