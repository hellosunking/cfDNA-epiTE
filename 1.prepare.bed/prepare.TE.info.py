# -- coding:UTF-8 --<code>
import sys
import gzip
#f=open('stat.txt','w')
out=open(sys.argv[1],'w')
for i in range(len(sys.argv)):
	filename=sys.argv[i]
	if (i>1):
		file = gzip.open(filename,'rt')
		num=0
		size=0
		for line in file :
			num=num+1
			s=int(line.strip().split('\t')[2])-int(line.strip().split('\t')[1])
			size=size+s
		out.write(filename.split('.')[0]+"\t"+filename.split('.')[1]+"\t"+str(num)+"\t"+str(size)+'\n')
		file.close()
out.close()


