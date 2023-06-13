import pandas as pd
import numpy as np
from scipy import stats
import os

def get_all_splice_events(samples):
	#merge all splice events together,remove no undefined strand splice event and standardize the data that suit for bedtools intersect to find overlap splice events.
	df=pd.read_csv('../%s/%sSJ.out.tab'%(samples[0],samples[0]),sep='\t',header=None)
	df=df[df[3]!=0]
	df=df[df[6]>=5]
	df=df[[0,1,2,3]]
	for sample in samples[1:]:
		tmp=pd.read_csv('../%s/%sSJ.out.tab'%(sample,sample),sep='\t',header=None)
		tmp=tmp[tmp[3]!=0]
		tmp=tmp[tmp[6]>=5]#only splice site with enough reads will be record
		tmp=tmp[[0,1,2,3]]
		df=pd.concat([df,tmp])
	df=df.drop_duplicates()
	df=df.sort_values([0,1,2])
	df[3]=df.apply(lambda x: ['+' if x[3]==1 else '-'][0],axis=1)
	df=df.reset_index()
	df['splice']=df.index +1
	df['splice']=df['splice'].astype('str')
	df['splice']='splice'+df['splice']
	df[5]='.'
	df=df[[0,1,2,'splice',5,3]]
	df.to_csv('all.splice.bed',sep='\t',index=None,header=None)

def build_splice_database():
	#site1:(1)splice event1 (2)splice event2 ....
	def rankit(x):
		x=x.reset_index()
		x[9]=x.index+1
		x[9]=x[9].astype('str')
		return x
	os.system('bedtools intersect -wa -wb -s -F 0.5 -a ath_gene.bed -b all.splice.bed >> all.splice.gene.bed')
	os.remove('all.splice.bed')
	df=pd.read_csv('all.splice.gene.bed',sep='\t',header=None)
	df=df.drop([0,1,2,4,5],axis=1)
	df=df.drop_duplicates(subset=[6,7,8,9,10,11])
	df[6]=df[3]+df[6]
	del df[3]
	df=df.sort_values([6,7,8])
	df.to_csv('all.splice.gene.bed',sep='\t',header=None,index=None)
	os.system('bedtools merge -S + -d 0 -i all.splice.gene.bed >> tmp1.data')#chain:+ merge
	os.system('bedtools merge -S - -d 0 -i all.splice.gene.bed >> tmp2.data')#chain:- merge
	df1=pd.read_csv('tmp1.data',sep='\t',header=None)
	df2=pd.read_csv('tmp2.data',sep='\t',header=None)
	df1[5]='+'#chain:+
	df2[5]='-'#chain:-
	df=pd.concat([df1,df2])
	df=df.sort_values([0,1,2])
	df=df.reset_index()
	df[3]=df.index+1
	df[3]=df[3].astype('str')
	df[3]='site'+df[3]
	df[4]='.'
	df=df[[0,1,2,3,4,5]]
	df.to_csv('tmp.bed',sep='\t',header=None,index=None)
	os.system('bedtools intersect -wa -wb -s -a tmp.bed -b all.splice.gene.bed >> site.data')
	os.remove('tmp1.data')
	os.remove('tmp2.data')
	os.remove('tmp.bed')
	fa=open('site.data')
	fb=open('splice.database','w')
	fb.write('Chr\tstart\tend\tsplice_id\tchain\n')
	line0=fa.readline()
	seq0=line0.rstrip().split('\t')
	site0=seq0[3].lstrip('site')
	fb.write(seq0[6]+'\t'+seq0[7]+'\t'+seq0[8]+'\tsplice'+site0+'.1\t'+seq0[11]+'\n')
	i=2
	for line in fa:
		seq=line.rstrip().split('\t')
		site=seq[3].lstrip('site')
		if site==site0:
			fb.write(seq[6]+'\t'+seq[7]+'\t'+seq[8]+'\tsplice'+site+'.'+str(i)+'\t'+seq[11]+'\n')
			i+=1
		else:
			fb.write(seq[6]+'\t'+seq[7]+'\t'+seq[8]+'\tsplice'+site+'.1\t'+seq[11]+'\n')
			site0=site
			i=2
	fa.close()
	fb.close()
	os.remove('site.data')
	os.remove('all.splice.gene.bed')
		
def calculate_splice_sites(samples):
	def cal_percent(tmp):
		if len(tmp)==1:
			s=['\t'.join(tmp[0])+'\t1.00\t1.00\t1.00\t1.00\t1.00\t1.00\t1.00\t1.00\t1.00\n']
			return s
		cal=[e[6:] for e in tmp]
		cal=np.array(cal)
		cal=cal.astype(int)
		per=cal/(0.00001+cal.sum(axis=0))
		per=per.round(2)
		for i in range(len(cal)):
			tmp[i]+=[str(e) for e in per[i]]
		s=['\t'.join(e)+'\n' for e in tmp]
		return s
	df=pd.read_csv('splice.database',sep='\t')
	df['gene']=df.apply(lambda x: x['Chr'][:-4],axis=1)
	df['Chr']=df.apply(lambda x: x['Chr'][-4:],axis=1)
	df=df[['gene','Chr','start','end','splice_id','chain']]
	for sample in samples:
		tmp=pd.read_csv('../%s/%sSJ.out.tab'%(sample,sample),sep='\t',header=None)
		tmp=tmp[tmp[3]!=0]
		tmp[3]=tmp.apply(lambda x: ['+' if x[3]==1 else '-'][0],axis=1)
		tmp=tmp[[0,1,2,3,6]]
		tmp=tmp.rename(columns={0:'Chr',1:'start',2:'end',3:'chain',6:sample})
		df=pd.merge(df,tmp,on=['Chr','start','end','chain'],how='left')#must use left in how
		df=df.fillna(0)
		df[sample]=df[sample].astype('int')
	df.to_csv('all.splice.result',sep='\t',index=None)
	fa=open('all.splice.result')
	fb=open('all.splice.percent.result','w')
	fb.write('gene\tChr\tstart\tend\tsplice\tchain\tCol0_rep1\tCol0_rep2\tCol0_rep3\tcyp71-3_rep1\tcyp71-3_rep2\tcyp71-3_rep3\tse-1_rep1\tse-1_rep2\tse-1_rep3\tpCol0_rep1\tpCol0_rep2\tpCol0_rep3\tpcyp71-3_rep1\tpcyp71-3_rep2\tpcyp71-3_rep3\tpse-1_rep1\tpse-1_rep2\tpse-1_rep3\n')
	next(fa)
	line0=fa.readline()
	seq0=line0.rstrip().split('\t')
	site0=seq0[4].split('.')[0]
	tmp=[seq0]
	for line in fa:
		seq=line.rstrip().split('\t')
		site=seq[4].split('.')[0]
		if site==site0:
			tmp.append(seq)
		else:
			s=cal_percent(tmp)
			fb.writelines(s)
			tmp=[seq]
			site0=site
	s=cal_percent(tmp)
	fb.writelines(s)
	fa.close()
	fb.close()
	os.remove('all.splice.result')

def calculate_pvalue():
	def pvalue(m,n,a,b):
		if sum(m)<=30 and sum(n)<=30:
			return 1
		if sum(a)==3 and sum(b)==3:#only 1 splice site
			return 1
		else:
			p=stats.ttest_ind(a,b)[1]
			return round(p,3)
	df=pd.read_csv('all.splice.percent.result',sep='\t')
	df['pCol0']=round((df['pCol0_rep1']+df['pCol0_rep2']+df['pCol0_rep3'])/3,3)
	df['pcyp71-3']=round((df['pcyp71-3_rep1']+df['pcyp71-3_rep2']+df['pcyp71-3_rep3'])/3,3)
	df['pse-1']=round((df['pse-1_rep1']+df['pse-1_rep2']+df['pse-1_rep3'])/3,3)
	df['fold_change(cyp71-3 vs Col-0)']=round(df['pcyp71-3']/df['pCol0'],2)
	df['pvalue(cyp71-3 vs Col-0)']=df.apply(lambda x:pvalue([x['Col0_rep1'],x['Col0_rep2'],x['Col0_rep3']],[x['cyp71-3_rep1'],x['cyp71-3_rep2'],x['cyp71-3_rep3']],[x['pCol0_rep1'],x['pCol0_rep2'],x['pCol0_rep3']],[x['pcyp71-3_rep1'],x['pcyp71-3_rep2'],x['pcyp71-3_rep3']]),axis=1)
	df['fold_change(se-1 vs Col-0)']=round(df['pse-1']/df['pCol0'],2)
	df['pvalue(se-1 vs Col-0)']=df.apply(lambda x:pvalue([x['Col0_rep1'],x['Col0_rep2'],x['Col0_rep3']],[x['se-1_rep1'],x['se-1_rep2'],x['se-1_rep3']],[x['pCol0_rep1'],x['pCol0_rep2'],x['pCol0_rep3']],[x['pse-1_rep1'],x['pse-1_rep2'],x['pse-1_rep3']]),axis=1)
	df.to_csv('all.splice.pvalue.bed',sep='\t',index=None)
	os.remove('all.splice.percent.result')

def AS_find():
	df=pd.read_csv('all.splice.pvalue.bed',sep='\t')
	df1=df[((df['pCol0']>=0.1)|(df['pcyp71-3']>0.1))&((df['fold_change(cyp71-3 vs Col-0)']>=1.25)|(df['fold_change(cyp71-3 vs Col-0)']<=0.8))&(df['pvalue(cyp71-3 vs Col-0)']<=0.05)]
	df2=df[((df['pCol0']>=0.1)|(df['pse-1']>0.1))&((df['fold_change(se-1 vs Col-0)']>=1.25)|(df['fold_change(se-1 vs Col-0)']<=0.8))&(df['pvalue(se-1 vs Col-0)']<=0.05)]
	df1['site']=df1.apply(lambda x: x['splice'].split('.')[0],axis=1)
	df2['site']=df2.apply(lambda x: x['splice'].split('.')[0],axis=1)
	df1=df1[['site']]
	df2=df2[['site']]
	df['site']=df.apply(lambda x: x['splice'].split('.')[0],axis=1)
	df1=pd.merge(df1,df,on='site',how='left')
	df2=pd.merge(df2,df,on='site',how='left')
	del df1['site']
	del df2['site']
	df1=df1.drop_duplicates()
	df2=df2.drop_duplicates()
	df1=df1.sort_values(['Chr','start','end','splice'])
	df2=df2.sort_values(['Chr','start','end','splice'])
	df1.to_csv('cyp71_AS.result',sep='\t',index=None)
	df2.to_csv('se_AS.result',sep='\t',index=None)
	fa=open('cyp71_AS.result')
	fb=open('se_AS.result')
	fc=open('AS_in_cyp71_and_se.result','w')
	lines1=fa.readlines()
	lines2=fb.readlines()
	list0=[line for line in lines1 if line in lines2]
	fc.writelines(list0)
	fa.close()
	fb.close()
	fc.close()

fa=open('../fastq/mRNA_samples.txt')
samples=[line.rstrip().split('\t')[1] for line in fa]
fa.close()
get_all_splice_events(samples)
build_splice_database()
calculate_splice_sites(samples)
calculate_pvalue()
AS_find()
