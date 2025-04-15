from sys import argv
sum=0
dic={}
type={}
out=open(argv[1],'w')
f=open(argv[2],"r")
label=["Bile","Breast","Colorectal","Gastric","Lung","Ovarian","Pancreatic"]

for line in f:
	if line.strip().split("\t")[1] != "Control" :
		dic[line.strip().split("\t")[0]]=[]
		type[line.strip().split("\t")[0]]=line.strip().split("\t")[1]

for i in range(7):
	fh=open(argv[i+3],'r')
	#fh.readline()
	for line in fh:
		ling=label[i]+"pred"
		line=line.replace('pred',ling)
		sid=line.strip().split(" ")[0]
		len=(line.strip().split("\t")[-1])
		if sid in dic.keys():
			dic[sid].append(len)

out.write("Sid\tType\t"+str("\t".join(label))+"\n")

for (x,y) in dic.items() :
	if x != "sid" :
		out.write(x+"\t"+type[x]+"\t"+str("\t".join(y))+"\n")

fh.close()
f.close()
out.close()
