import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df=pd.read_csv('./cyp71.txt',sep='\t',index_col='miRNAid')
df['cyp71-2/Col0']=df['cyp71-2']/df['Col-0']
df['cyp71-3/Col0']=df['cyp71-3']/df['Col-0']
df['Col0']=np.log10(df['Col-0']+1)
df['cyp71-2']=np.log10(df['cyp71-2']+1)
df['cyp71-3']=np.log10(df['cyp71-3']+1)
df=df[['Col0','cyp71-2','cyp71-3','cyp71-2/Col0','cyp71-3/Col0','pvalue_cyp71-2','pvalue_cyp71-3']]

texts=['ath-miR159a','ath-miR159b','ath-miR160a','ath-miR165a','ath-miR167a','ath-miR171a','ath-miR172a']
v=1.5#set value

t=np.arange(0,5,0.2)
################plot cyp71-2################
#scatter
cyp712_up=df[df['cyp71-2/Col0']>=v]
cyp712_down=df[(df['cyp71-2/Col0']>0) & (df['cyp71-2/Col0']<=1/v)]
cyp712_down=cyp712_down[cyp712_down['pvalue_cyp71-2']<=0.05]
cyp712_middle=df[(df['cyp71-2/Col0']<v) & (df['cyp71-2/Col0']>1/v)]

plt.figure(figsize=(6,6),dpi=80)
labels=['cyp71-2/Col0>=%s'%v,'Col0/cyp71-2>=%s'%v,'Ratio<%s'%v]

up2=plt.scatter(cyp712_up['Col0'],cyp712_up['cyp71-2'],color='darkgrey',marker='.')
down2=plt.scatter(cyp712_down['Col0'],cyp712_down['cyp71-2'],color='red',marker='.')
middle2=plt.scatter(cyp712_middle['Col0'],cyp712_middle['cyp71-2'],color='darkgrey',marker='.')
plt.plot(t,t,'--',color='darkgrey')
for index in texts:
    plt.text(df.loc[index][0],df.loc[index][1],index,fontsize=4,horizontalalignment='left')

plt.xlim(0,5)
plt.ylim(0,5)
plt.xlabel('log10(Col0 expression+1)')
plt.ylabel('log10(cyp71-2 expression+1)')

plt.savefig('scatter_cyp71-2.pdf')
plt.close()



################plot cyp71-3################
#scatter
cyp713_up=df[df['cyp71-3/Col0']>=v]
cyp713_down=df[(df['cyp71-3/Col0']>0) & (df['cyp71-3/Col0']<=1/v)]
cyp713_down=cyp713_down[cyp713_down['pvalue_cyp71-3']<=0.05]
cyp713_middle=df[(df['cyp71-3/Col0']<v) & (df['cyp71-3/Col0']>1/v)]

plt.figure(figsize=(6,6),dpi=80)
labels=['cyp71-3/Col0>=%s'%v,'Col0/cyp71-3>=%s'%v,'Ratio<%s'%v]

plt.scatter(cyp713_up['Col0'],cyp713_up['cyp71-3'],color='darkgrey',marker='.')
plt.scatter(cyp713_down['Col0'],cyp713_down['cyp71-3'],color='red',marker='.')
plt.scatter(cyp713_middle['Col0'],cyp713_middle['cyp71-3'],color='darkgrey',marker='.')
plt.plot(t,t,'--',color='darkgrey')
for index in texts:
    plt.text(df.loc[index][0],df.loc[index][2],index,fontsize=4,horizontalalignment='left')

plt.xlim(0,5)
plt.ylim(0,5)
plt.xlabel('log10(Col0 expression+1)')
plt.ylabel('log10(cyp71-3 expression+1)')

plt.savefig('scatter_cyp71-3.pdf')
plt.close()


