function tempw=Estimate_temp_ES(temp, depCTD, p)
% Estimate the temperature around the whale
% Based on function: Dsw=EstimateDsw(SWdensity, depCTD, p)
% temp: output from function SWdensityFromCTD_ES
% depCTD: CTD's depth data where SWdensity was recorded
% p: animal's depth data
xi=p;
x= depCTD;
y= temp;

tempw = interp1(x,y,xi,'linear','extrap');

end