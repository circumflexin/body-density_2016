#pragma rtGlobals=1		// Use modern global access method and strict wave access.

/////Macros to extract strokes: ExtractStroke_1 & ExtractStroke_2
/// use in order


//Extract all positive and negative peaks in heave specific acceleration
//then, make histogram to find out threshold for further analysis

/// set variables (fs, s_dur) before running the macro
///Updated by Tomoko Narazaki (19 May 2016)

Function ExtractStroke_1(heaveSA, fs, s_dur)
wave heaveSA
variable fs, s_dur
//variable fs			
//fs = 4		//fs: sampling frequency (Hz)
print "sampling freq of specific acceleration =", fs, "Hz"
//variable s_dur		//s_dur: minimum stroke duration in seconds
//s_dur = 2
print "minimum stroke duration is set as", s_dur, "seconds"
variable p, numPoints
numPoints = numpnts(heaveSA)
duplicate/O heaveSA diffSA stroke peaks
diffSA[0] = 0
stroke = nan
peaks = nan
p=1
do
diffSA[p] = heaveSA [p] - heaveSA[p-1]
p+=1
while(p<=numPoints)

//find peaks
p=1
do
	if(diffSA[p]*diffSA[p+1]<0)
		peaks[p] = heaveSA[p]
	else
		peaks[p] = nan
	endif
p+=1
while(p<=numPoints)

//remove peaks that exist within minimum stroke duration
p=1
do
	wavestats/Q/R = [p -s_dur*fs/2, p+s_dur*fs/2] heaveSA
	if((peaks[p]<V_max)&&(peaks[p]>0))				//remove errors in positive peaks
		peaks[p] = nan
	elseif((peaks[p]>V_min)&&(peaks[p]<0))			//remove errors in negative peaks
		peaks[p] = nan
	endif
p+=1
while(p<=numPoints)

Duplicate/O peaks stroke
killwaves diffSA

//make histogram
peaks = abs(peaks)
Make/N=500/O peaks_Hist;DelayUpdate
Histogram/P/B=4 peaks,peaks_Hist
Display peaks_Hist
ModifyGraph mode=5

END




///Extract strokkes using the threshold
/// set variable (stroke_def) before running the macro
///Updated by Tomoko Narazaki (19 May 2016)

Function ExtractStroke_2(stroke, stroke_def)
wave stroke		
variable stroke_def
variable p
//stroke_def = 0.05			//stroke_def: threshold for stroking
print "stroke definition =", stroke_def, "G"
p=0
do
	if(abs(stroke[p])< stroke_def)
		stroke[p] = nan
	endif
p+=1
while(p<=numpnts(stroke))
END
