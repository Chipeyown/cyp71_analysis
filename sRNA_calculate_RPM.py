import sys
def reversea(s):
	s1=''
	for i in s[::-1]:
		if i=='A':
			s1+='T'
		elif i=='T':
			s1+='A'
		elif i=='G':
			s1+='C'
		elif i=='C':
			s1+='G'
	return s1

pre=sys.argv[1]

dict={}
fa=open('%s.sam'%pre)
for line in fa:
	if line.startswith('@'):
		continue
	seq=line.split('\t')
	if seq[1]=='4':
		continue
	if seq[1]=='0':
		sequence=seq[9]
	else:
		sequence=reversea(seq[9])
	try:
		dict[sequence]+=1
	except:
		dict[sequence]=1
fa.close()

with open('shells/%s.mapped.log'%pre) as f:
	n=int(f.readlines()[-1].split()[1])

fb=open('%s.normalized.txt'%pre,'w')
fb.write('Sequence\tRPM\n')
for key,value in dict.items():
	fb.write(key+'\t'+str(round(value*1000000/n,3))+'\n')
fb.close()

