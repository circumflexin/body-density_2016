function Dsw=EstimateDsw(SWdensity, depCTD, p)

%Dsw=EstimateDsw(SWdensity, depCTD, p)
%SWdensity: seawater density from CTD measurement
%depCTD: CTD's depth data where SWdensity was recorded
%p: animal's depth data
xi=p;
x= depCTD;
y= SWdensity;

Dsw = interp1(x,y,xi, 'linear','extrap');

end