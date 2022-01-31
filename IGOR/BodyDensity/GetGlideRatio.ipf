#pragma rtGlobals=1		

/// Get Glide ratio during descent and ascent phase
/// Descent glide ratio is obtained for each dive (from the time at the animal reached 'DeepDef' to  the start  of bottom phase)
/// Ascent glide ratio is obtained for each dive (from the end of bottom phase to the time at the animal reached 'DeepDef')

// NOTE: sampling frequency must be the same for D, SGtype, pitch
//Set variables (DeepDef, fs)

//Updated by Tomoko Narazaki (May 20, 2016)


Function GetGlideRatio (Dive, Bottom, SGtype, pitch, D, deepDef, fs)
wave  Dive, Bottom, SGType, pitch, D
variable deepDef, fs
variable nd , p, dep, Tdep
//variable DeepDef, fs
//DeepDef = 100			//DeepDef: get glide ratio deeper than this value
print "get glide ratio for >", DeepDef, "m"
//fs = 4					//fs: sampling frequency in Hz
print "sampling frequency of D and SGtype", fs, "Hz"

Make/O/D/N = (numpnts(DiveNumber), 10) G_ratio
G_ratio = nan
Duplicate/O D deepG
deepG=0

nd=0
do
///Descent phase
	if(Bottom[nd][1] >=DeepDef)			//get descent glide ratio only if bottom start depth >= DeepDef
		Tdep=Dive[nd][0]
		do
			if(D(Tdep)>=DeepDef)
				break
			endif
		Tdep+=1
		while(Tdep<=Bottom[nd][0])		
		deepG[(Tdep-leftx(D))*fs, (Bottom[nd][0]-leftx(D))*fs] =1
		dep=D(Tdep)
		wavestats/Q /R = (Tdep, Bottom[nd][0]) SGType
		G_ratio[nd][0] = V_npnts/fs								//G_ratio[][0]: total duration during descent phase(s)
		G_ratio[nd][1] = V_sum/fs									//G_ratio[][1]: glide duration during descent (s)
		G_ratio[nd][2] = G_ratio[nd][1]/G_ratio[nd][0]				//G_ratio[][2]: glide ratio during descent		
		wavestats/Q /R = (Tdep, Bottom[nd][0]) pitch
		G_ratio[nd][3] = V_avg										//G_ratio[][3]: mean pitch during descent (degree) 
		G_ratio[nd][4] = (Bottom[nd][1] - dep)/G_ratio[nd][0]		//G_ratio[][4]: descent rate (m/s)				
	endif
///Ascent phase
	if(Bottom[nd][3] >=DeepDef)		//get ascent glide ratio only if bottom end depth >= DeepDef
		Tdep=Dive[nd][1]
		do
			if(D(Tdep)>=DeepDef)
				break
			endif
		Tdep-=1
		while(Tdep>=Bottom[nd][2])	
		deepG[(Bottom[nd][2]-leftx(D))*fs, (Tdep-leftx(D))*fs] =1	
		dep=D(Tdep)
		wavestats/Q /R = (Bottom[nd][2], Tdep) SGtype
		G_ratio[nd][5] = V_npnts/fs								//G_ratio[][5]: total duration during ascent phase(s)
		G_ratio[nd][6] = V_sum/fs									//G_ratio[][6]: glide duration during ascent (s)
		G_ratio[nd][7] = G_ratio[nd][6]/G_ratio[nd][5]				//G_ratio[][7]: glide ratio during ascent		
		wavestats/Q /R = (Bottom[nd][2], Tdep) pitch
		G_ratio[nd][8] = V_avg										//G_ratio[][8]: mean pitch during ascent (degree)	
		G_ratio[nd][9] = (dep - Bottom[nd][3])/G_ratio[nd][5]		//G_ratio[][9]: ascent rate (m/s)					
	endif	
nd+=1
while(nd<=numpnts(DiveNumber))
killwaves DeepG

END
		