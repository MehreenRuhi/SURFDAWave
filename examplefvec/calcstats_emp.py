import sys
import os
fname=sys.argv[1]

sites=[]
data=[]
def hetfunc(datasamp):
        heter=0
        for tt in range(len(datasamp[0])):
               cou=0
               for ttt in range(len(datasamp)):
                        cou+=float(datasamp[ttt][tt])
               p=cou/float(len(datasamp))


               heter+=(2*p)*(1-p)
        return heter




with open(fname,'r') as i:
	for j in i:
		popinfor=[]
		if j[0]!='#':
			geneline=j.rstrip('\n').split('\t')[9:]
			for eachind in geneline:
				for hap in eachind.split('|'):
					popinfor.append(hap)
		if len(set(popinfor))==2:
			sites.append(j.split('\t')[1])
			data.append(''.join(popinfor))




f2=open(fname+'.stats','a')
f3=open(fname+'.sites','a')

from collections import Counter

import numpy as np
popsize=len(data[1])
allstats=[]
nsites=[]
for t in range(0,(len(sites)/5)):
	statslist=[]
        nsites.append(sites[t*5])
     

        d= data[t*5:(t*5)+10]

        d=[map(None, dd) for dd in d]

        hapmat=np.matrix(d)

        hapmat2=hapmat.T

        haplis=np.array(hapmat2).tolist()
        haplist=[''.join(y) for y in haplis]
        count=Counter(haplist)
        counts=[]

        for l in set(haplist):
                counts.append(count[l])

        avepi=hetfunc(haplist)/len(haplist[0])
	
     	statslist.append(avepi)
	
        sortedcount= sorted(counts,reverse=True)
        hlist=[]
        for each in sortedcount:
                hlist.append(each)
        hlist2=[eachH/(float(len(haplist))) for eachH in hlist]
        lenhlist=len(hlist)
        hapnum=0
        hapfreq=[]
        while hapnum<5 and hapnum<lenhlist:
                hapfreq.append(hlist[hapnum]/float(popsize))
                hapnum+=1
        hsqlist=[j ** 2 for j in hlist2]

        h2sq=hsqlist[1:]
        h2sum=float(sum(h2sq))
        h1=float(sum(hsqlist))
        h21=h2sum/h1
        one=hlist2[0]
        two=hlist2[1]
        h12sqlist=h2sq[1:]
        h12part1=(one+two)*(one+two)
        h12part2=sum(h12sqlist)
        h12=h12part1+h12part2
	statslist.append(h1)
	statslist.append(h12)
	statslist.append(h21)
        for eachhapfreq in hapfreq:
              		statslist.append(eachhapfreq)
        if lenhlist<5:
                numzeros=5-lenhlist
                for eachzero in range(numzeros):
                       
			statslist.append(0.0)
	allstats.append(statslist)
	

for eachstat in range(len(allstats)-128):
        
        haps=allstats[eachstat:eachstat+128]
	flathaps=[item for sublist in haps for item in sublist]

        for d in range(len(flathaps)):
                f2.write(str(flathaps[d])+',')
        f2.write('\n')
	f3.write(str(nsites[eachstat+64])+'\n')

f2.close()
f3.close()







			
