#pragma rtGlobals=1

// Estimate seawater density around the tagged animal
// from seawater density profile obtained from CTD measurement
// Updated by Tomoko Narazaki (May 20, 2016)

Function EstimateDsw (SWdensity, depCTD, D)
wave SWdensity, depCTD, D		//SWdensity: seawater density from CTD measurement
								//depCTD: CTD's depth data where SWdensity was recorded
								//D: animal's depth data
variable p
Duplicate/O D Dsw
p=0
do
	Dsw[p] = interp(D[p], depCTD, SWdensity)	//Dsw: seawater density around the tagged animal
p+=1
while(p<=numpnts(D))
END

