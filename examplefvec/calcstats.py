from collections import Counter
import sys

def hetfunc(datasamp):
	heter=0
        for tt in range(len(datasamp[0])):
 	       cou=0
               for ttt in range(len(datasamp)):
			cou+=float(datasamp[ttt][tt])
               p=cou/float(len(datasamp))


               heter+=(2*p)*(1-p)
	return heter

file1=sys.argv[1]
f2=open(file1+'.stats','a')
f=open(file1,'r')
g=f.readlines()
popsize=int(g[0].split(' ')[1])
pop=popsize+2
for i in range(len(g)):
	segs=[]
	data=[]
	if g[i][0:3]=='seg':
		
		segs=g[i+1].rstrip(' \n').lstrip('positions: ').split(' ')
		segs=map(float,segs)
		seglistlen=len(segs)
		if seglistlen>660:
			mid= seglistlen/2
			data=g[i+2:i+pop]
			oneten=seglistlen/10
			ranges=range(-320,330,5)	
			for win in range(128):
				midsnp=[]
				winsegs=[]
				for eachline in data:
					midsnp.append(eachline[mid])
					winsegs.append(eachline[mid+ranges[win]:mid+ranges[win+2]])
					
				pi=hetfunc(winsegs)
				avepi=pi/len(winsegs[0])
				f2.write(str(avepi))
				f2.write(',')
					
		                count=Counter(winsegs)
                                counts=[]

                                for l in set(winsegs):
                                        counts.append(count[l])
                                sortedcount= sorted(counts,reverse=True)
				hlist=[]
                                for each in sortedcount:
                                        hlist.append(each)
				lenhlist=len(hlist)
				hapnum=0
				hapfreq=[]
				
				while hapnum<5 and hapnum<lenhlist:
					hapfreq.append(hlist[hapnum]/float(popsize))
					hapnum+=1
				
					
 					
				hlist2=[eachH/(float(len(winsegs))) for eachH in hlist]
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
                                f2.write(str(h1))
                                f2.write(',')
                                f2.write(str(h12))
                                f2.write(',')
                                f2.write(str(h21))
                                f2.write(',')
				for eachhapfreq in hapfreq:
					f2.write(str(eachhapfreq)+',')
				if lenhlist<5:
					numzeros=5-lenhlist
					for eachzero in range(numzeros):
						f2.write(str(0.0)+',')
					
			f2.write('\n')
		else:
			pass

		

f.close()	
f2.close()
