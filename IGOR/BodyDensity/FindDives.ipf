#pragma rtGlobals=1		// Use modern global access method and strict wave access.

////Extract dives
/// set variables (DiveDef, D_DiveDef, fs) before running the macro
///Updated by Tomoko Narazaki (19 May 2016)

Function FindDives (D, diveDef, d_diveDef, fs)
wave D
variable diveDef, d_diveDef, fs
variable p,  numPoints, nd, numDives
numPoints=numpnts(D)
variable Tdep, dep		//Tdep: time od deepest point in dive; dep: dive depth

//variable DiveDef
//DiveDef = 2			//DiveDef: threshold of a dive
print "Dives are defined as any submergences  >", DiveDef, "m"
//variable D_DiveDef
//D_DiveDef = 10			//D_DiveDef: consider as 'deep dives' if dives are deeper than D_DiveDef
print "Find deep dives only of which max depth >", D_DiveDef, "m"
//variable fs
//fs = 4					//fs: sampling frequency of DEPTH data
print "sampling freq of depth data =", fs, "Hz"

Make/O/D/N = (numpnts(D), 6) Dive
Duplicate/O D D1
Dive = nan
D1 = 0

//Discriminate "dive" or "not dive"
p=1
do
	if(D[p] > DiveDef)
		D1[p] = 1		//D1=1: when it is during dive; D1=0: when it is NOT during dive
	endif
p+=1
while(p<=numPoints)

//Search the start of each dive
p=1
nd =0
do
	if(D1[p] - D1[p-1]==1)	
		Dive[nd][0] = p*(1/fs) + leftx(D)	//Dive[][0]: start time of dive (i.e. # of seconds from 1904/1/1)	
		nd+=1
	endif
p+=1
while(p<=numPoints)

print nd
numDives = nd-1			//counting the total no. of dives
						//print numDives						
DeletePoints/M = 0 numDives+1, numPoints-numDives, Dive
						//delete the excess parts from wave "Dive"
						
//Searching of end of each dive
p=0
nd=0
do
	if(D1[p+1] - D1[p]==-1)		
		Dive[nd][1] = p*(1/fs) + leftx(D)		//Dive[][1]: end time of dive (i.e. # of seconds from 1904/1/1)	
	nd+=1
	endif
p+=1
while(nd<=numDives)				

//get dive duration, post-dive surface duration, max dive depth & the deepest point of dive
nd=0
do
	Dive[nd][2] = Dive[nd][1]-Dive[nd][0] +1*(1/fs)		//Dive[][2]: Dive duration in s
	Dive[nd][3] = Dive[nd+1][0] - Dive[nd][1]				//Dive[][3]: Post-dive surface durationi in s
	
	Tdep = Dive[nd][0]
	dep = DiveDef
	do
		if(D(Tdep)>=dep)
			Dive[nd][4] = Tdep				//Dive[][4]: Time of the deepest point od the dive
											//(i.e. # of seconds from 1904/1/1)
			dep = D(Tdep)
			Dive[nd][5] = dep				//Dive[][5]: max dive depth
		endif
	Tdep+=1
	while(Tdep<=Dive[nd][1])
nd+=1
while(nd<numDives)					

//get  dive duraion etc of the last dive
Dive[numDives][2] = Dive[numDives][1]-Dive[numDives][0]+1*(1/fs)
Dive[numDives][3] = numPoints*deltax(D) + leftx(D) - Dive[numDives][1]

Tdep=Dive[numDives][0]
dep = DiveDef
do
	if(D(Tdep)>=dep)
		Dive[numDives][4] = Tdep
		
		dep = D(Tdep)
		Dive[numDives][5] = dep
	endif
Tdep+=1
while(Tdep<=Dive[numDives][1])


//Remove SHALLOW dives
nd = numDives
do
	if(Dive[nd][5]<=D_DiveDef)
		DeletePoints/M=0 nd, 1, Dive
	endif
nd-=1
while(nd>=0)



//Make waves for further analysis
Make/D/O/N = (numpnts(Dive)/6) Start, Fin, DiveNumber, MaxD, MaxTime, StartT, FinT, Duration, StartD, FinD
nd=0
do
	Start[nd] = (Dive[nd][0] - leftx(D)	)*fs		
	Fin[nd] = (Dive[nd][1] - leftx(D))*fs
	StartT[nd] = Dive[nd][0]
	FinT[nd] = Dive[nd][1]
	StartD[nd] = D(Dive[nd][0])
	FinD[nd] = D(Dive[nd][1])
	MaxD[nd] = Dive[nd][5]
	MaxTime[nd]=Dive[nd][4]
	DiveNumber[nd] = nd +1
	Duration[nd] = Dive[nd][2]
nd+=1
while(nd<numpnts(Dive)/6)
killwaves D1
End


