function[SWdensity,depCTD, temp]=SWdensityFromCTD(DPT, TMP, SL,D); % ES - made temp an output
%lat = 70.82611;                      % CORRECT for CTD; latitude in decimal degrees

newdepth=0:max(DPT); % ES - changed from 0 to 1
temp=NaN(size(newdepth));
sali=NaN(size(newdepth));
temp(round(DPT))=TMP;
sali(round(DPT))=SL;
%linear interpret between the NaNs 
temp=fixgaps(temp);
sali=fixgaps(sali);
P = gsw_p_from_z_ES(-newdepth,lat);  % ES - Added this line to calculate pressure from depth
dens = sw_dens(sali,temp,P);         % ES - Changed from sw_dens0 to sw_dens so as to include pressure (depth)
figure(5); clf;
plot(dens,newdepth);
ylabel('Depth (m)')
xlabel('Density (kg/m^3)')
if max (D(:,6))> max(DPT)
    newdepth2=0:max(D(:,6));  % ES - changed from 0 to 1
    sali=[sali,(sali(end)*ones((length(newdepth2)-length(newdepth)),1))'];
    temp=[temp,(temp(end)*ones((length(newdepth2)-length(newdepth)),1))'];
    dens = sw_dens0(sali,temp);
    %    dens = sw_dens(sali,temp,P);     % ES - Changed from sw_dens0 to sw_dens so as to include pressure (depth)
    figure(5); clf;
    plot(dens,newdepth2)
    ylabel('Depth (m)')
    xlabel('Density (kg/m^3)')
end 
depCTD=newdepth2;
SWdensity=dens;
end