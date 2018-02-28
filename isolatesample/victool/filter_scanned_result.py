import os,sys

f = open("goodrecohandscan_1mu1p.csv")
ll = f.readlines()
f.close()

for l in ll[2:202]:
    l = l.strip()
    info = l.split(",")
    r = int(info[0])
    s = int(info[1])
    e = int(info[2])
    good = int(info[5])
    mode = info[-2]
    
    if good==0 or mode!="1mu1p":
        continue
    #if good==0:
    #    continue

    #print info
    #print r,s,e,good,mode
    print "%d,%d,%d"%(r,s,e)


