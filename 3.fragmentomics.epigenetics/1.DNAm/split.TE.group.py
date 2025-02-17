# -- coding:UTF-8 --<code>
import sys
#f=open('stat.txt','w')
out=open(sys.argv[1],'w')
for i in range(len(sys.argv)):
	filename=sys.argv[i]
	if (i>1):
		file = open(filename,'r')
		num=0
		size=0
		for line in file :
			num=num+1
			s=int(line.strip().split('\t')[2])-int(line.strip().split('\t')[1])
			size=size+s
		out.write(filename.replace(".bed","")+"\t"+filename+"\t"+str(num)+"\t"+str(size)+'\n')
		file.close()
out.close()


