function[SWdensity,depCTD]=SWdensityFromCTD(DPT, TMP, SL,D);

newdepth=0:max(DPT);
temp=NaN(size(newdepth));
sali=NaN(size(newdepth));
temp(round(DPT))=TMP;
sali(round(DPT))=SL;
%linear interpret between the NaNs 
temp=fixgaps(temp);
sali=fixgaps(sali);
dens = sw_dens0(sali,temp);
figure(5); clf;
plot(dens,newdepth);
ylabel('Depth (m)')
xlabel('Density (kg/m^3)')
if max (D(:,6))> max(DPT)
    newdepth2=0:max(D(:,6));
    sali=[sali,(sali(end)*ones((length(newdepth2)-length(newdepth)),1))'];
    temp=[temp,(temp(end)*ones((length(newdepth2)-length(newdepth)),1))'];
    dens = sw_dens0(sali,temp);
    figure(5); clf;
    plot(dens,newdepth2)
    ylabel('Depth (m)')
    xlabel('Density (kg/m^3)')
end 
depCTD=newdepth2;
SWdensity=dens;
end