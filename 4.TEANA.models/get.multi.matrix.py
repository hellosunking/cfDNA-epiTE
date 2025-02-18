from sys import argv
sum=0
dic={}
type={}
out=open(argv[1],'w')
f=open(argv[2],"r")
lable=["Bile","Breast","Colorectal","Gastric","Lung","Ovarian","Pancreatic"]
for line in f:
	if  line.strip().split("\t")[1] != "Control" :
		dic[line.strip().split("\t")[0]]=[]
		type[line.strip().split("\t")[0]]=line.strip().split("\t")[1]
for i in range(7):
	fh=open(argv[i+3],'r')
	#fh.readline()
	for line in fh:
		ling=lable[i]+"pred"
		line=line.replace('pred',ling)
		sid=line.strip().split(" ")[0]
		print(sid)
		len=(line.strip().split("\t")[-1])
		#print(len)
		if sid in dic.keys():
			dic[sid].append(len)
		#len=line.strip().split("\t")[1:]
		#print(sid)
		# sum=sum+line.split("\t")[2]
#print(len)
for  (x,y) in dic.items() :
	if x=="Sid":
		out.write(x+"\t"+type[x]+"\t"+str("\t".join(y))+"\n")
for  (x,y) in dic.items() :
	if x!="Sid":
		out.write(x+"\t"+type[x]+"\t"+str("\t".join(y))+"\n")
fh.close()
f.close()
out.close()
