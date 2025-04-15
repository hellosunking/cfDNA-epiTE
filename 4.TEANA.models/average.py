from sys import argv
sum=0
dic={}
out=open(argv[101],'w')
for i in range(100):
    fh=open(argv[i+1],'r')
    fh.readline()
    for line in fh:
        sid="   ".join(line.split("\t")[0:2])
        if sid in dic.keys() :
            dic[sid]=dic[sid]+float(line.strip().split("\t")[2])
        else:
            dic[sid]=float(line.strip().split("\t")[2])
out.write("Sid"+"\t"+"Type"+"\t"+"pred"+"\n")
for  (x,y) in dic.items() :
    n=float(y/10)
    out.write(x+"\t"+str(n)+"\n")
fh.close()
out.close()

