#pragma rtGlobals=1		

////Separate a dive into descent, bottom and ascent phase
/// The end of descent is searched from the start of dive - the first point where animals' pitch change to positive
/// The start of ascent is searched from the end of dive - the last point where animal's pitch is negative

//Set variables (D_DiveDef, Pdef, fs) before running the macro
///Updated by Tomoko Narazaki (19 May 2016)


Function FindBottom (D, pitch, Dive, d_diveDef, pitchDef, fs)
wave D, pitch, Dive
variable d_diveDef, pitchDef, fs
variable p, nd, numDives
numDives = numpnts(Dive)/6-1
variable Tdep, dep, depB, TdepB
//variable D_DiveDef
//D_DiveDef = 10
print "Deep dive definition =", D_DiveDef, "m"
//variable PitchDef, fs
//PitchDef = 1				//Pdef: pitch definition
print "pitch definition  =", PitchDef, "degree"
//fs = 4				//fs: sampling frequency of depth data (Hz)
print "sampling frequency of depth data =", fs, "Hz"

Make/O/D/N=(numpnts(Dive)/6, 4) Bottom

//Search the start of bottom phase
nd = 0
do
	if (Dive[nd][5] < D_DiveDef)
		Bottom[nd][0] = NaN			//  in case of shallow dive
		Bottom[nd][1] = NaN			//  in case of shallow dive
	else
		depB = D_DiveDef		//bottom phase should start deeper than deepDiveDef
		TdepB = Dive[nd][0]
		do
			if (D(TdepB) >= depB)
				break
			endif
		TdepB+=1
		while(TdepB<=Dive[nd][4])								
			
		Tdep = TdepB		//search from >10 m
		do
			if(pitch(Tdep) >=-(PitchDef))
				break
			endif	
		Tdep += 1
		while (Tdep <= Dive[nd][4])	
		Bottom[nd][0] = Tdep		// Bottom[nd][0]: time of bottom start
		dep = D(Tdep)
		Bottom[nd][1] = dep			// Bottom[nd][1]: depth at bottom start
	endif				
	if (Bottom[nd][0] == 0)
		Bottom[nd][0] = NaN
		Bottom[nd][1] = NaN
	endif
nd += 1		
while (nd <= numDives)


//Search the end of bottom phase
nd = 0
do
	if (Dive[nd][5] < D_DiveDef)
		Bottom[nd][2] = 0			
		Bottom[nd][3] = 0			
	else	
	depB = D_DiveDef			//bottom phase should end deeper than DeepDiveDef
	TdepB = Dive[nd][1]
	do
		if (D(TdepB) >= depB)
			break
		endif
	TdepB -= 1
	while(TdepB>=Dive[nd][4])

	Tdep = TdepB
	do
		if(pitch(Tdep)<=PitchDef)
			break
		endif
	Tdep -= 1
	while (Tdep >= DIve[nd][4])			
	Bottom[nd][2] = Tdep	// Bottom[nd][2]: time of bottom end
	dep = D(Tdep)
	Bottom[nd][3] = dep		// Bottom[nd][3]: depth at bottom end
	endif		
	if (Bottom[nd][2] == 0)
		Bottom[nd][2] = NaN
		Bottom[nd][3] = NaN
	endif
nd += 1		
while (nd <= numDives)
	
duplicate/O D Phase
phase = nan				//phase[]=nan: not during dive
nd=0
do
	phase[(Dive[nd][0] - leftx(D))*fs, (Bottom[nd][0]- leftx(D))*fs] = -1		//Phase[]=-1: descent
	phase[(Bottom[nd][0]- leftx(D))*fs, (Bottom[nd][2]- leftx(D))*fs] = 0		//phase[] = 0: bottom
	phase[(Bottom[nd][2]- leftx(D))*fs, (Dive[nd][1]- leftx(D))*fs]= 1		//phase[] = 1: ascent
nd+=1
while(nd<=numDives)		
End
 


