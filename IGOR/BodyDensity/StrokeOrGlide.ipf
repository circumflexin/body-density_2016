#pragma rtGlobals=1		// Use modern global access method and strict wave access.


//Make a wave indicating whether during stroke or glide (SGtype)
	//SGtype = 1; stroking occured during a certain interval
	//SGtype = 0; stroking NOT occured within a certain interval
/// set variable (fs) before running the macro
///Updated by Tomoko Narazaki (19 May 2016)

Function StrokeOrGlide (Stroke, interval, fs)
wave Stroke
variable interval, fs
variable p, q
variable fs_D

//interval = 1		//interval: examine whether stroke occured during this period
print "check if stroke occured during +/-", interval, "seconds"

//fs = 4				//fs: sampling frequency of acceleration data
print "sampling freq of accelration =", fs, "Hz"

duplicate/O stroke StrokeOrNot SGtype
StrokeOrNot = nan
SGtype = nan
 
 p=0
 do
 	if(numtype(stroke[p])==2)
 		StrokeOrNot[p] = 0
 	elseif(numtype(stroke[p])==0)
 		StrokeOrNot[p] = 1
 	endif
 p+=1
 while(p<=numpnts(stroke))
 
 p=interval
 do
 	wavestats/Q/R = [(p-interval*fs), (p+interval*fs)]StrokeOrNot
 	if(V_sum==0)
 		SGtype[p] = 1
 	elseif(V_sum > 0)
 		SGtype[p] = 0
 	endif
p+=1
while (p<=numpnts(SGtype))
killwaves StrokeOrNot
END





