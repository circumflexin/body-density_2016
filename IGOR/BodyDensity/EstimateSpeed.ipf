#pragma rtGlobals=1


//Estimate speed 
//NOTE: sampling freq of pitch & D must be the same to use EstimateSpeed_1
//Updated by Tomoko Narazaki (May19, 2016)

Function EstimateSpeed (pitch, D, fc, fs, pitchDef)		//use pitch in degree
wave pitch, D
variable fc, fs, pitchDef
variable p//, fc
//fc = 5				//fc: smoothing parameter for vertical speed. default is 5s
//variable fs
//fs = 4				//fs: sampling freq of depth data (Hz)
print "sampling frequency of depth data =", fs, "Hz"
//variable pitchDef
//pitchDef = 30		//pitchDef: threshold to estimate swim speed
print "estimate swim speed only when pitch >", pitchDef, "degree"

duplicate/O D vertSp SwimSp
vertSp = nan
p=1
do
	vertSp[p] = (D[p] - D[p-1])*fs		//vertSp: vertical speed in m/s
p+=1
while(p<=numpnts(D))
Smooth fs*fc, vertSp
duplicate/O pitch sinP
sinP = sin(pitch*pi/180)
p=0
do
	if(abs(pitch[p]) > pitchDef)
		SwimSp[p] = -vertSp[p]/sinP[p]
	else
		SwimSp[p] = nan
	endif
p+=1
while(p<=numpnts(D))
//remove negative speed 
p=0
do
	if(SwimSp[p]<0)
		SwimSp[p] = nan
	endif
p+=1
while(p<=numpnts(SwimSp))
killwaves sinP
END

