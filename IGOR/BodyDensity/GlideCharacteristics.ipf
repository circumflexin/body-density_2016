#pragma rtGlobals=1		

// All glides are segmented to sub-glides (i.e. all sub-glides have the same duration)
// Then, characteristics of each sub-glides are summarized for further analysis
// NOTE: sampling frequency must be the same for D, Temp, SwimSP, Dsw,pitch, roll, head
//Set variables (dur, Ddef, D_DiveDef, fs)

// If your data does not have temperature, do the following before running the macro:
// Duplicate D Temp
// Temp = Nan

//Updated by Tomoko Narazaki (May 20, 2016)

Function SubGlides_summary (SGtype, D, SwimSP, Temp, Dive, DiveNumber, Start, Fin, Dsw, pitch, roll, head, sub_dur, Ddef, fs)
wave SGtype, D, SwimSP, Temp, Dive, DiveNumber, Start, Fin, Dsw, pitch, roll, head
variable sub_dur, Ddef, fs
variable p, ng, nd
variable numPoints, numGlides, numDives
numPoints = numpnts(D)
numDives = numpnts(DiveNumber)
variable VP, dep, XX, YY

variable D_DiveDef
//sub_dur = 5			//sub_dur: duration of sub-glides in second
print "Glides are segmented into ", sub_dur, "-second sub-glides"
//Ddef = 0				//Ddef: any glides shallower than this is excluded
print "Extract  glides only when  >", Ddef,"m deep"
//D_DiveDef = 1000			//D_DiveDef: Deep dive definition
//print "consider deep dives when max depth >", DeepDiveDef, "m"
//fs = 4				//fs: sampling frequency

Make/O/D/N = (numpnts(D), 28) Glide
Glide = nan

//search start of each sub=glide
duplicate/O D G
G = 0
p=1
do
	if(G[p-1] ==0)			//when p-1 was not glide
	if(SwimSp[p]>0)			// speed data is available
	wavestats/Q/R = [p, p+sub_dur*fs-1] SGtype
		if(V_sum  > sub_dur*fs-1)				
			G[p, p+sub_dur*fs*2-1] = 1			//G = 1: sub-glide
			Glide[ng][0] = p					//Glide[][0]: sub-glide start point
			Glide[ng][1] = p+sub_dur*fs-1	//Glide[][1]: sub-glide end point
		ng+=1
		endif
	endif
	endif
p+=1
while(p<=numPoints)

numGlides = ng -1		//counting the total no. of sub-glides
						
DeletePoints/M = 0 numGlides+1, numPoints-numGlides, Glide
					//delete the excess parts of "Glide"

////sub-glide duration (for checking)
G=0
ng = 0
do
	Glide[ng][2] = (Glide[ng][1] - Glide[ng][0]+1)/fs
	G[Glide[ng][0], Glide[ng][1]] = 1
ng+=1
while(ng<=numGlides)

///////////mean depth, speed, temp, pitch (&sd), temp & seawater density
ng =0
do
	wavestats/Q/R = [Glide[ng][0], Glide[ng][1]] D
		Glide[ng][3] = V_avg					//Glide[][3]: mean depth
	wavestats/Q/R = [Glide[ng][0], Glide[ng][1]] SwimSp
		Glide[ng][5] = V_avg					//Glide[][5]: mean SwimSp
	wavestats/Q/R = [Glide[ng][0], Glide[ng][1]] pitch
		Glide[ng][6] = V_avg					//Glide[][6]: mean pitch 
		Glide[ng][7] = sin(V_avg*pi/180)		//Glide[][7]: sin of pitch
		Glide[ng][8] = V_sdev				//Glide[][8]: sd of pitch

//added by taiki//////////////////////////////////////////////////
		Glide[ng][24] = V_max - V_min	// //Glide[][25]: changes in pitch
	wavestats/Q/R = [Glide[ng][0], Glide[ng][1]] roll
		Glide[ng][25] = V_max - V_min			//Glide[][25]: changes in roll
		Glide[ng][26] = V_avg			//Glide[][26]: mean roll
////////////////////////////////////////////////////////////////
		
	wavestats/Q/R = [Glide[ng][0], Glide[ng][1]] temp
		Glide[ng][9] = V_avg					//Glide[][9]: mean temp
	wavestats/Q/R = [Glide[ng][0], Glide[ng][1]] Dsw
		Glide[ng][10] = V_avg				//Glide[][10]: mean seawater density	
ng+=1
while(ng<=numGlides)

//////depth change during glides
ng = 0
do
	Glide[ng][4] = (D[Glide[ng][1]] - D[Glide[ng][0]])		//Glide[][4]: total depth change during the glide (m)
ng+=1
while(ng<=numGlides)

///// acceleration during the glide
Make/D/O/N = (2) W_sigma, W_coef
ng = 0
do
	duplicate/R=[Glide[ng][0], Glide[ng][1]]SwimSp vG
	wavestats/Q vG
	if(V_numNans <=0)				//get acceleration only when SwimSp is available throuought the glide
		CurveFit/Q/NTHR=0 line  vG /D
		VP=abs(V_pr)
		Glide[ng][12]=V_Pr				//Glide[][12]:The linear coefficient r (Pearson's r) 
		Glide[ng][13] = W_sigma[1]		//Glide[][13]: SE of the gradient
                                                               
		if(VP>0)					
			Glide[ng][11]= W_coef[1] //K1			//Glide[][11]:acceleration (i.e. the gradient of the regression line)
		else
			Glide[ng][11]= 0				
		endif
	else
		Glide[ng][11] = nan			// no acceleration calculated when SwimSp data is missing during the glide
		Glide[ng][12] = nan			// NAN for r coefficient, too.		
	endif                                                                
killwaves vG
ng+=1
while(ng<=numGlides)

//descent, bottom or ascent phase?
ng=0
do
	wavestats/Q/R = [Glide[ng][0], Glide[ng][1]] Phase
	if(V_numNans ==0)
		if(V_Sum>0)
			Glide[ng][14] = 1	//Glide[][14] = 1: ascent glide
		elseif(V_sum==0)
			Glide[ng][14] = 0	//Glide[][14] = 0: bottom glide
		elseif(V_sum < 0)
			Glide[ng][14] = -1	//Glide[][14] = -1: descent glide
		endif
	elseif(V_numNans > 0)
		Glide[ng][13] = nan		// Glide[][14] = nan : the glide was not during dives
	endif
ng+=1
while(ng<=numGlides)

// Dive number, depth, durtaion
duplicate D num
num = 0
duplicate num dp dr
nd = 0
do
	num[Start[nd], Fin[nd]] = DiveNumber[nd]
	dp[Start[nd], Fin[nd]] = Dive[nd][5]
	dr[Start[nd], Fin[nd]] = DIve[nd][2]
nd+=1
while(nd<=numDIves)

ng=0
do
	wavestats/Q/R = [Glide[ng][0], Glide[ng][1]] num
	Glide[ng][15] = V_max				//Glide[][15]: Dive Number in wihch the glide occoured
	wavestats/Q/R = [Glide[ng][0], Glide[ng][1]] dp
	Glide[ng][16] = V_max				//Glide[][16]: Max dive depth (m) of the dive
	wavestats/Q/R = [Glide[ng][0], Glide[ng][1]] dr
	Glide[ng][17] = V_max				//Glide[][17]: Dive duration (s) of the dive
ng+=1
while(ng<=numGlides)

///////mean and variations in pitch using circular statistics
Make/D/N = (sub_dur*fs) angle
angle = nan
duplicate angle sinAngle cosAngle
ng=0
do
	p=0
	do
		angle[p] = Pitch[Glide[ng][0]+p]
		sinAngle[p] = sin(angle[p]*pi/180)
		cosAngle[p] = cos(angle[p]*pi/180)
	p+=1
	while(p<=numpnts(angle))	
	wavestats/Q sinAngle
	YY = V_Sum/(sub_dur*fs)
	wavestats/Q cosAngle
	XX = V_Sum/(sub_dur*fs)
	Glide[ng][19] = sqrt(XX^2+YY^2)			
					//Glide[][19]: measure of concentration (r) of pitch angle during the sub-glide
	if(YY>=0 && XX >=0)
		Glide[ng][18] = asin(YY/sqrt(XX^2+YY^2))*180/pi		//Glide[][18]: mean pitch calculated using circular stats
	elseif(YY>=0 && XX < 0)
		Glide[ng][18] = 180 - asin(YY/sqrt(XX^2+YY^2))*180/pi
	elseif(YY< 0 && XX >=0)
		Glide[ng][18] = asin(YY/sqrt(XX^2+YY^2))*180/pi
	elseif(YY < 0 && XX < 0)
		Glide[ng][18] =  -180 - asin(YY/sqrt(XX^2+YY^2))*180/pi 
	endif
ng+=1
while(ng<=numGlides)

///////mean and variations in roll using circular statistics
angle = nan
ng=0
do
	p=0
	do
		angle[p] = roll[Glide[ng][0]+p]
		sinAngle[p] = sin(angle[p]*pi/180)
		cosAngle[p] = cos(angle[p]*pi/180)
	p+=1
	while(p<=numpnts(angle))	
	wavestats/Q sinAngle
	YY = V_Sum/(sub_dur*fs)
	wavestats/Q cosAngle
	XX = V_Sum/(sub_dur*fs)
	Glide[ng][21] = sqrt(XX^2+YY^2)			
				//Glide[][21]: measure of concentration (r) of roll angle during the sub-glide
	if(YY>=0 && XX >=0)
		Glide[ng][20] = asin(YY/sqrt(XX^2+YY^2))*180/pi		//Glide[][20]: mean roll
	elseif(YY>=0 && XX < 0)
		Glide[ng][20] = 180 - asin(YY/sqrt(XX^2+YY^2))*180/pi
	elseif(YY< 0 && XX >=0)
		Glide[ng][20] = asin(YY/sqrt(XX^2+YY^2))*180/pi
	elseif(YY < 0 && XX < 0)
		Glide[ng][20] =  -180 - asin(YY/sqrt(XX^2+YY^2))*180/pi 
	endif
ng+=1
while(ng<=numGlides)

///////mean and variations in heading using circular statistics
angle = nan
ng=0
do
	p=0
	do
		angle[p] = head[Glide[ng][0]+p]
		sinAngle[p] = sin(angle[p]*pi/180)
		cosAngle[p] = cos(angle[p]*pi/180)
	p+=1
	while(p<=numpnts(angle))	
	wavestats/Q sinAngle
	YY = V_Sum/(sub_dur*fs)
	wavestats/Q cosAngle
	XX = V_Sum/(sub_dur*fs)	
	Glide[ng][23] = sqrt(XX^2+YY^2)		
			//Glide[][23]: measure of concentration (r) of heading angle during the sub-glide
	if(YY>=0 && XX >=0)
		Glide[ng][22] = asin(YY/sqrt(XX^2+YY^2))*180/pi		//Glide[][22]: mean heading
	elseif(YY>=0 && XX < 0)
		Glide[ng][22] = 180 - asin(YY/sqrt(XX^2+YY^2))*180/pi
	elseif(YY< 0 && XX >=0)
		Glide[ng][22] = asin(YY/sqrt(XX^2+YY^2))*180/pi
	elseif(YY < 0 && XX < 0)
		Glide[ng][22] =  -180 - asin(YY/sqrt(XX^2+YY^2))*180/pi 
	endif
ng+=1
while(ng<=numGlides)

//Remove very shallow glides
ng=numGlides
do
	if(Glide[ng][3]< Ddef)
		DeletePoints/M=0 ng, 1, Glide
	endif
ng-=1
while(ng>=0)
killwaves num dp dr angle sinAngle cosAngle

//Glide[][0]: sub-glide start pt
//Glide[][1]: sub-glide end pt
//Glide[][2]:sub-glide duration (s)
//Glide[][3]: mean depth (m) during the sub-glide
//Glide[][4]: total depth change during the sub-glide (m)
//Glide[][5]: mean SwimSp (m/s)
//Glide[][6]: mean pitch (degree)
//Glide[][7]: SIN of pitch
//Glide[][8]: SD of pitch
//Glide[][9]: mean temperature
//Glide[][10]: mean seawater density (Dsw)
//Glide[][11]: acceleration during the sub-glide
//Glide[][12]: R-value for the regression swim speed vs time
//Glide[][13]: SE of the gradient for the regression swim speed vs time
//Glide[][14]: 0 for bottom phase, -1 for descent phase, 1 for ascent phase
//Glide[][15]: dive number in which the sub-glide occured
//Glide[][16]: max dive depth(m) of the dive
//Glide[][17]: dive duration (s) of the dive
//Glide[][18]: mean pitch (deg) calculated using circular statistics
//Glide[][19]: Measure of concentration (r) of pitch durin the sub-glide
//Glide[][20]: mean roll (deg) calculated using circular statistics
//Glide[][21]: Measure of concentration (r) of roll durin the sub-glide
//Glide[][22]: mean heading (deg) calculated using circular statistics
//Glide[][23]: Measure of concentration (r) of heading durin the sub-glide

//added by taiki on 2019/03/19
//Glide[][24]: changes in pitch
//Glide[][25]: changes in roll
//Glide[][26]: mean roll (degree)

//added by taiki on 2019/05/08
//Glide[][27]: stable glide or not: stable = 1, not = 0, calculated later of this function

Glide[][27] = 0

End
	
Function 	SelectStableGlide(Glide, D)
	wave Glide, D
	variable noUnstableGlide = 0
	//make DeltaDepth
	DepthToDeltaDepth(D)
	print "not stable glide is at Glide ID:"
	variable i
	for(i = 0; i < DimSize(Glide, 0); i+=1) // i for Glide
		wavestats/q/r=[Glide[i][0], Glide[i][1]] DeltaDepth
		if(V_min * V_max > 0) // i.e. continuous ascent or descent during 1 sec
			Glide[i][27] = 1
		else
			noUnstableGlide+=1
			Glide[i][27] = 0
		endif
	endfor
	print "No. of all glides = " && 	DimSize(Glide, 0)
	print "No. of unstable glides = " && noUnstableGlide
End

Function DepthToDeltaDepth(D)
	wave D
	variable i
	Duplicate/O D, DeltaDepth
	DeltaDepth = NaN
	for(i = 0; i < DimSize(D, 0); i+=1) // i for D
		DeltaDepth[i] = D[i+1] - D[i]
	endfor
End