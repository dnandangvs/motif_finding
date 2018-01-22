#!/usr/bin/python3
fhot=open('hotspots','r')
fcold=open('coldspots','r')
fsand=open('flanked','w')
fnon=open('unflanked','w')
fhmap=open('HOTSPOTS_MAP.csv','r')
pop=dict()
val=''

for line in fhmap:
	words=line.split()
	temp=words[0]+'\t'+words[1]+'\t'+words[2]
	
	# for i in range(3,len(words)):
	# 	if(i==len(words)-1):
	# 		val=val+words[i]
	# 	else:
	# 		val=val+words[i]+'\t'
	pop[temp]=words[2]+'\t'+words[3]+'\t'+words[4]+'\t'+words[5]+'\t'+words[6]+'\t'+words[7]+'\t'+words[8]+'\t'+words[9]+'\t'+words[10]+'\t'+words[11]+'\t'+words[12]+'\t'+words[13]+'\t'+words[14]+'\t'+words[15]+'\t'+words[16]+'\t'+words[17]+'\t'+words[18]+'\t'+words[19]+'\t'+words[20]+'\t'+words[21]+'\t'+words[22]+'\t'+words[23]+'\t'+words[24]+'\t'+words[25]+'\t'+words[26]+'\t'+words[27]

print('reading HOSTP finished')

fsand.write('chr'+'\t'+'HS_start'+'\t'+'HS_end'+'\t'+'HS_rank'+'\t'+'pCS_start'+'\t'+'pCS_end'+'\t'+'pCS_rank'+'\t'+'pCS_end-HS_start'+'\t'+'fCS_start-HS_end'+'\t'+'Mean'+'\t'+'fCS_start'+'\t'+'fCS_end'+'\t'+'fCS_rank'+'\t'+'Mean CS_rank'+'\t'+'HS_rank-pCS_rank'+'\t'+'HS_rank-fCS_rank'+'\t'+'HS_rank-meanCS_rank'+'\t'+'width'+'\t'+'Validated'+'\t'+'Conserved'+'\t'+'COS_rate'+'\t'+'CEU_rate'+'\t'+'TSI_rate'+'\t'+'FIN_rate'+'\t'+'GBR_rate'+'\t'+'IBS_rate'+'\t'+'YRI_rate'+'\t'+'LWK_rate'+'\t'+'ASW_rate'+'\t'+'ACB_rate'+'\t'+'MKK_rate'+'\t'+'CHB_rate'+'\t'+'JPT_rate'+'\t'+'CHS_rate'+'\t'+'CDX_rate'+'\t'+'KHV_rate'+'\t'+'MXL_rate'+'\t'+'PUR_rate'+'\t'+'CLM_rate'+'\t'+'PEL_rate'+'\t'+'GIH_rate'+'\t'+'REFINED_rate'+'\n')
fnon.write('chr'+'\t'+'HS_start'+'\t'+'HS_end'+'\t'+'HS_rank'+'\t'+'pCS_start'+'\t'+'pCS_end'+'\t'+'pCS_rank'+'\t'+'pCS_end-HS_start'+'\t'+'fCS_start-HS_end'+'\t'+'Mean'+'\t'+'fCS_start'+'\t'+'fCS_end'+'\t'+'fCS_rank'+'\t'+'Mean CS_rank'+'\t'+'HS_rank-pCS_rank'+'\t'+'HS_rank-fCS_rank'+'\t'+'HS_rank-meanCS_rank'+'\t'+'width'+'\t'+'Validated'+'\t'+'Conserved'+'\t'+'COS_rate'+'\t'+'CEU_rate'+'\t'+'TSI_rate'+'\t'+'FIN_rate'+'\t'+'GBR_rate'+'\t'+'IBS_rate'+'\t'+'YRI_rate'+'\t'+'LWK_rate'+'\t'+'ASW_rate'+'\t'+'ACB_rate'+'\t'+'MKK_rate'+'\t'+'CHB_rate'+'\t'+'JPT_rate'+'\t'+'CHS_rate'+'\t'+'CDX_rate'+'\t'+'KHV_rate'+'\t'+'MXL_rate'+'\t'+'PUR_rate'+'\t'+'CLM_rate'+'\t'+'PEL_rate'+'\t'+'GIH_rate'+'\t'+'REFINED_rate'+'\n')

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
fhot.close()
fcold.close()
test=0

for i in range(1,23):
	chrom="chr"+str(i)
	
	hotspots[chrom].sort()
	coldspots[chrom].sort()
	start=0
	
	for y in range(len(coldspots[chrom])-1):
		pcstart,pcend=coldspots[chrom][y]
		fcstart,fcend=coldspots[chrom][y+1]
		tmp=0
		for x in range(start,len(hotspots[chrom])):
			hstart,hend=hotspots[chrom][x]
			if(y==0):
				if(hend<=pcstart):
					test=test+1
					# print(chrom+'\t'+str(hstart)+'\t'+str(hend))
					# print(hotspots[chrom][x+1])
					# print(pcstart,pcend)
					hs=chrom+'\t'+str(hstart)+'\t'+str(hend)
					pcs=chrom+'\t'+str(pcstart)+'\t'+str(pcend)
					fcs=chrom+'\t'+str(fcstart)+'\t'+str(fcend)
					pdis=pcend-hstart
					fdis=fcstart-hend
					fnon.write(hs+'\t'+str(hranks[hs])+'\t'+str(pcstart)+'\t'+str(pcend)+'\t'+str(cranks[pcs])+'\t'+str(pdis)+'\t'+str(fdis)+'\t'+str((pdis+fdis)/2)+'\t'+str(fcstart)+'\t'+str(fcend)+'\t'+str(cranks[fcs])+'\t'+str((cranks[pcs]+cranks[fcs])/2)+'\t'+str(hranks[hs]-cranks[pcs])+'\t'+str(hranks[hs]-cranks[fcs])+'\t'+str(hranks[hs]-((cranks[pcs]+cranks[fcs])/2))+pop[hs]+'\n')

			if(y==len(coldspots[chrom])-2):
				if(hstart>=fcstart):
					test=test+1
					# print(chrom+'\t'+str(hstart)+'\t'+str(hend))
					# print(fcstart,fcend)
					hs=chrom+'\t'+str(hstart)+'\t'+str(hend)
					pcs=chrom+'\t'+str(pcstart)+'\t'+str(pcend)
					fcs=chrom+'\t'+str(fcstart)+'\t'+str(fcend)
					pdis=pcend-hstart
					fdis=fcstart-hend
					fnon.write(hs+'\t'+str(hranks[hs])+'\t'+str(pcstart)+'\t'+str(pcend)+'\t'+str(cranks[pcs])+'\t'+str(pdis)+'\t'+str(fdis)+'\t'+str((pdis+fdis)/2)+'\t'+str(fcstart)+'\t'+str(fcend)+'\t'+str(cranks[fcs])+'\t'+str((cranks[pcs]+cranks[fcs])/2)+'\t'+str(hranks[hs]-cranks[pcs])+'\t'+str(hranks[hs]-cranks[fcs])+'\t'+str(hranks[hs]-((cranks[pcs]+cranks[fcs])/2))+pop[hs]+'\n')

			if(hstart>=pcend and hend<=fcstart):
				tmp=tmp+1
			else:
				
				if(tmp==1):
					count=count+1
					hstart,hend=hotspots[chrom][x-1]
					hs=chrom+'\t'+str(hstart)+'\t'+str(hend)
					pcs=chrom+'\t'+str(pcstart)+'\t'+str(pcend)
					fcs=chrom+'\t'+str(fcstart)+'\t'+str(fcend)
					pdis=pcend-hstart
					fdis=fcstart-hend
					fsand.write(hs+'\t'+str(hranks[hs])+'\t'+str(pcstart)+'\t'+str(pcend)+'\t'+str(cranks[pcs])+'\t'+str(pdis)+'\t'+str(fdis)+'\t'+str((pdis+fdis)/2)+'\t'+str(fcstart)+'\t'+str(fcend)+'\t'+str(cranks[fcs])+'\t'+str((cranks[pcs]+cranks[fcs])/2)+'\t'+str(hranks[hs]-cranks[pcs])+'\t'+str(hranks[hs]-cranks[fcs])+'\t'+str(hranks[hs]-((cranks[pcs]+cranks[fcs])/2))+pop[hs]+'\n')
					
				if(tmp>1):
					for pos in range(tmp):

						hstart,hend=hotspots[chrom][x-pos-1]
						hs=chrom+'\t'+str(hstart)+'\t'+str(hend)
						pcs=chrom+'\t'+str(pcstart)+'\t'+str(pcend)
						fcs=chrom+'\t'+str(fcstart)+'\t'+str(fcend)
						pdis=pcend-hstart
						fdis=fcstart-hend
						fnon.write(hs+'\t'+str(hranks[hs])+'\t'+str(pcstart)+'\t'+str(pcend)+'\t'+str(cranks[pcs])+'\t'+str(pdis)+'\t'+str(fdis)+'\t'+str((pdis+fdis)/2)+'\t'+str(fcstart)+'\t'+str(fcend)+'\t'+str(cranks[fcs])+'\t'+str((cranks[pcs]+cranks[fcs])/2)+'\t'+str(hranks[hs]-cranks[pcs])+'\t'+str(hranks[hs]-cranks[fcs])+'\t'+str(hranks[hs]-((cranks[pcs]+cranks[fcs])/2))+pop[hs]+'\n')
				hstart,hend=hotspots[chrom][x]
					
				if(hstart>fcstart):
					start=x
					break
			if(x==len(hotspots[chrom])-1):
				if(tmp==1):
					count=count+1
					fsand.write(hs+'\t'+str(hranks[hs])+'\t'+str(pcstart)+'\t'+str(pcend)+'\t'+str(cranks[pcs])+'\t'+str(pdis)+'\t'+str(fdis)+'\t'+str((pdis+fdis)/2)+'\t'+str(fcstart)+'\t'+str(fcend)+'\t'+str(cranks[fcs])+'\t'+str((cranks[pcs]+cranks[fcs])/2)+'\t'+str(hranks[hs]-cranks[pcs])+'\t'+str(hranks[hs]-cranks[fcs])+'\t'+str(hranks[hs]-((cranks[pcs]+cranks[fcs])/2))+pop[hs]+'\n')
			
					
print(count,test)
print(hcount,ccount)

