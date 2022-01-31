#pragma rtGlobals=1

//Estimate seawater desity from CTD data
// This macro requires the following waves: 
	//DPT (measured CTD depth), TMP (measured CTD temperature) & SAL (measured CTD salinity)
// Set the variable 'Deepest' before running the macro. 
//'Deepest' must be deeper than animal's max dive depth.
// updated by Tomoko Narazaki (May 20, 2016)

Function SWdensityFromCTD(DPT, TMP, SL, deepest)
wave DPT, TMP, SL		//data from CTD measurement
variable deepest	
variable p
variable MaxDep			//CTD's deepest measurement
wavestats/Q DPT
MaxDep = trunc(V_max)
print "deepest measurement in CTD was", MaxDep, "m"

//variable Deepest
//Deepest = 2000			//Deepest: must be deeper than animal's max dive depth
print "create CTD plot up to", Deepest, "m"
Make/D/O/N = (MaxDep) dep
dep = nan
duplicate/O dep sal tempr

// calculate average temperature and salinity for each depth (per metre)
p=1
do
	extract/O DPT, dest, DPT>=p&&DPT<p+1
	if(numpnts(dest)>0)				//when there is CTD measurement at this depth range
		wavestats/Q dest
		dep[p] = V_avg				///dep: average depth per 1m interval
	
		extract/O TMP, dest, DPT>=p&&DPT<p+1
		wavestats/Q dest
		tempr[p] = V_avg			//tempr: average temp per 1m interval
	
		extract/O SL, dest, DPT>=p&&DPT<p+1
		wavestats/Q dest
		sal[p] = V_avg				//sal: average salinity per 1 m interval
	else
		dep[p] = nan
		tempr[p] = nan
		sal[p] = nan
	endif	
p+=1
while(p<=MaxDep)
killwaves  dest

// Assume temp and salinity is constant at deep depth
// if animals dive deeper than CTD measurement range.
if (Deepest > MaxDep)
	insertpoints numpnts(dep), 1, dep, sal, tempr
endif
dep[MaxDep+1] = Deepest
sal[MaxDep+1] = sal[MaxDep-1]
tempr[MaxDep+1] = tempr[MaxDep-1]

/////calculate seawater density using equation of state of seawater
 /// based on UNESCO (1983)
duplicate/O dep F3 G3 H3 I3 N3 AA BB P3 Q3 S3 T3 SWdensity
Q3=2*(9.8*dep+1)*10^-5
F3 = (1+0.1*dep)*1.01325
G3 =999.842594+(6.793952*10^-2)*tempr+(-9.09529*10^-3)*tempr^2+(1.001685*10^-4)*tempr^3+(-1.120083*10^-6)*tempr^4+(6.536322*10^-9)*tempr^5
H3 = G3 + sal*(0.824493+(-4.0899*10^-3)*tempr+(7.6438*10^-5)*tempr^2+(-8.2467 *10^-7)*tempr^3+(5.3875*10^-9)*tempr^4)+sal^1.5*(-5.72466*10^-3+(1.0227*10^-4)*tempr +(-1.6546*10^-6)*tempr^2)+sal^2*(4.8314*10^-4)
I3 = 19652.21+148.4206*tempr+(-2.327105)*tempr^2+(1.360477*10^-2)*tempr^3+(-5.155288*10^-5)*tempr^4
N3 = I3 + sal*(54.6746-0.603459*tempr+(1.09987*10^-2)*tempr^2+(-6.167*10^-5)*tempr^3)+sal^1.5*(7.944*10^-2+(1.6483*10^-2)*tempr+(-5.3009*10^-4)*tempr^2)
AA = 3.239908+(1.43713*10^-3)*tempr+(1.160921*10^-4)*tempr^2+(-5.77905*10^-7)*tempr^3 + ((2.2838*10^-3)+(-1.0981*10^-5)*tempr+(-1.6078*10^-6)*tempr^2)*sal+(1.91075*10^-4)*sal^1.5
BB = 8.50935*10^-5+(-6.12293*10^-6)*tempr+(5.2787*10^-8)*tempr^2 + ((-9.9348*10^-7)+(2.0816*10^-8)*tempr +(9.1697*10^-10)*tempr^2)*sal
P3 = N3 + AA*F3 + BB*F3^2
S3=-4*H3*(P3*(9.8*dep+1)*10^-5)
T3=SQRT(P3^2+S3)

SWdensity =(P3-T3)/Q3
killwaves  F3 G3 H3 I3 N3 AA BB P3 Q3 S3 T3
//p = F3; ƒÏw =G3; ƒÏ(S,T,0) = H3; Kw = I3; K(S, T, 0) =N3;  K(S,T,p) = P3
duplicate/O dep depCTD
killwaves dep sal tempr

/// Remove any NANs from depCTD and SWdensity
p = numpnts(depCTD)
do
	if(numtype(SWdensity[p])==2)
	DeletePoints/M=0 p, 1, SWdensity
	DeletePoints/M=0 p, 1, depCTD	
	endif
p-=1
while(p>=0)		

display depCTD vs SWdensity
SetAxis/A/R left
ModifyGraph mode=4,marker=19,msize=2,rgb=(0,0,0)
Label left "Depth";DelayUpdate
Label bottom "SWdensity"

End	

