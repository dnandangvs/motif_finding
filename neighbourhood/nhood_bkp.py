#!/usr/bin/python3
fhot=open('hotspots','r')
fcold=open('coldspots','r')
fsand=open('flanked','w')
fnon=open('unflanked','w')
hotspots=dict()
for i in range(1,23):
	hotspots["chr"+str(i)]=[]
coldspots=dict()
for i in range(1,23):
	coldspots["chr"+str(i)]=[]
hcount=0
ccount=0
count=0
hranks=dict()
cranks=dict()

for line in fhot:
	words=line.split()
	hotspots[words[0]].append((int(words[1]),int(words[2])))
	hcount=hcount+1
	hranks[line.strip()]=hcount
	
	

for line in fcold:
	words=line.split()
	coldspots[words[0]].append((int(words[1]),int(words[2])))
	ccount=ccount+1
	cranks[line.strip()]=ccount
for i in range(1,23):
	chrom="chr"+str(i)
	
	hotspots[chrom].sort()
	coldspots[chrom].sort()
	start=0

	for y in range(len(coldspots[chrom])-1):
		dummy,cend=coldspots[chrom][y]
		cstart,dummy=coldspots[chrom][y+1]
		tmp=0
		for x in range(start,len(hotspots[chrom])):
			hstart,hend=hotspots[chrom][x]
			if(hstart>=cend and hend<=cstart):
				tmp=tmp+1
			else:
				
				if(tmp==1):
					count=count+1
					hpos=chrom+'\t'+str(hstart)+'\t'+str(hend)
					fsand.write(hpos+'\t'+str(hranks[string])+'\n')
					break
				if(hstart>cstart):
					start=x
					break
			if(x==len(hotspots[chrom])-1):
				if(tmp==1):
					count=count+1
			
					
print(count)
