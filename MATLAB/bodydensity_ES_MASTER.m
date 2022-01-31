%%%Body density estimation. Japan workshop May 2016
%Lucia Martina Martin Lopez lmml2@st-andrews.ac.uk
% Edited by Eilidh Siegal, July 2017 (es250@st-andrews.ac.uk)
% All edits by ES are commented, with comments preceded by "ES - "

clc             % ES - added this line                 
clear           % ES - added this line   
addpath(genpath('Y:\Eilidh\PhD\Chapter3\Part_1_glide_extraction\MatLab\BodyDensity_MATLAB\bodydensity_tools'))  % ES - added this line   
addpath(genpath('Y:\Eilidh\PhD\Chapter3\Part_1_glide_extraction\MatLab\BodyDensity_MATLAB'))                    % ES - added this line   
addpath(genpath('Y:\Eilidh\PhD\Chapter3\Part_1_glide_extraction\MatLab\BodyDensity_MATLAB\d2matlab_tools'))     % ES - added this line   
addpath(genpath('Y:\Eilidh\PhD\Chapter3\Part_1_glide_extraction\MatLab\BodyDensity_MATLAB\d3matlab_tools'))     % ES - added this line   

%% 1_LOAD DATA
% Load prh
tag='';                 % Insert deployment here
settagpath('prh','Y:\tag\tag2\metadata\prh\not25')
loadprh(tag);           


%% 2_DEFINE DIVES and make a summary table describing the characteristics of each dive.

% ES - used 2m as the threshold of the dive (surface) and 10m for dive definition (i.e. deep dive if below 10m)
% as in the IgorPro FindDives.ipf provided by Tomoko Narazaki (edited last 19 May 2016)

k=1:length(p);
mindivedef = 10;                        % ES - edited dive defintiion to keep same as in Miller et al. 2016
surface=2;                              % ES - edited dive defintiion to keep same as in Miller et al. 2016
T=finddives(p,fs,mindivedef,surface); 
D=[];                                   % [start_time(s) end_time(s) dive_duration(s) max_depth(s) max_depth(m) ID(n)]
D(:,1)=T(:,1);                          % start time of each dive in seconds since the tag on time
D(:,2)=T(:,2);                          % end time of each dive in seconds since the tag on time
D(:,3)=T(:,2)-T(:,1);                   % dive duration in seconds
D(:,4)=[T(2:end,1)-T(1:end-1,2);NaN];   % post-dive surface duration in seconds
D(:,5)=T(:,4);                          % time of the deepest point of each dive
D(:,6)=T(:,3);                          % maximum dive depth of each dive
D(:,7)=1*(1:size(T,1))';                % dive ID number, starting from 1


%% 3_SEPARATE LOW AND HIGH ACCELERATION SIGNALS 
% 3.1_QUICK SEPARATION OF DESCENT AND ASCENT PHASES

BOTTOM=[];
DES=[];
ASC=[];
for dive=1:size(T,1)
    kk=round(fs*T(dive,1)):round(fs*T(dive,2));                 % it is selecting the whole dive
    enddes=(find((pitch(kk)*180/pi)>0,1,'first')+T(dive,1)*fs); % search for the first point at which pitch is positive
    startasc=(find((pitch(kk)*180/pi)<0,1,'last')+T(dive,1)*fs);% search for the last point at which the pitch is negative
   BOTTOM(:,1)=(enddes)/fs;                                     % time in seconds at the start of bottom phase (end of descent)
   BOTTOM(:,3)=(startasc)/fs;                                   % time in seconds at the end of bottom phase (start of descent)
   des=round(fs*T(dive,1)):round(enddes);                       % selects the whole descent phase
   asc=round(startasc):round(fs*T(dive,2));                     % selects the whole ascent phase
   DES=[DES,des];
   ASC=[ASC,asc];
end


%% 3.2_SEPARATE LOW AND HIGH ACCELERATION SIGNALS 

%3.2.1_select periods of analysis.
% The power spectra is calculated of the longitudinal and 
% dorso-ventral accelerometer signals during descents and ascents to 
% determine the dominant stroke frequency for each animal in each phase with
% a fft of 512 and a sampling rate of fs. 
% Output: S is the amount of power in each particular frequency (f)

% During whole deployment 
[S,f]=speclev(Aw(k,:),512,fs);plot(f,S);
% During all descents and ascents phases where mainly steady swimming
% occurs. When calculated for the whole dive it may be difficult to 
% differenciate the peak at which stroking rate occurs as there is other 
% kind of movements than only steady swimming
[S,f]=speclev(Aw(DES,:),512,fs);plot(f,S); % During descent (ES added this line)
[S,f]=speclev(Aw(ASC,:),512,fs);plot(f,S); % During ascent (ES added this line)

%3.2.2_Determine FR fluking rate and cut-off frequency
% calculate the power spectrum of the accelerometer data at the whale frame
% plot PSD for the whole deployment, all descents and all ascent phases.
FR = nan(6,1); FRmag = FR;
FRl = nan(6,1); FRlmag = FRl;
FRn = 1;
cs = 'bgr';
va = {'top';'';'bottom'};

% Subplot 1: whole deployment
figure (1); clf; 
[S,f]=speclev(Aw(k,:),512,fs), grid;
ax(1)=subplot(311);
for nn = [1 3]
    [peakloc, peakmag] = peakfinder(S(:,nn),0.1); % make 0.1 smaller if it is not finding the peak you expect
    peakloc(1) = []; peakmag(1) = [];
    smoothS = runmean(S(:,nn),10);
    peakmag = peakmag - smoothS(peakloc);
    [~,peak] = max(peakmag);
    peak = peakloc(peak);
    plot(f(peak),S(peak,nn),'go','markersize',10,'linewidth',2);
    [minf,bf] = min(S(1:peak,nn));
    FRl(FRn) = f(bf);
    hold on
    plot([f(bf) f(bf)],[min(S(:,nn)) minf],'k--','linewidth',2);
    text(f(bf),min(min(S(:,[1 3]))),['f = ' num2str(round(f(bf)*100)/100)],'horizontalalignment','right','color',cs(nn),'verticalalignment',va{nn});
    FR(FRn) = f(peak);
    FRmag(FRn) = S(peak,nn);
    FRn = FRn+1;
end
[~,b] = max(FRmag(1:2));
if b == 2; v1 = 'top'; v2 = 'bottom'; else v1 = 'bottom'; v2 = 'top'; end
text(FR(1),FRmag(1)+2*sign(length(v1)-length(v2)),num2str(round(FR(1)*100)/100),'verticalalignment',v1,'horizontalalignment','center');
text(FR(2),FRmag(2)-2*sign(length(v1)-length(v2)),num2str(round(FR(2)*100)/100),'verticalalignment',v2,'horizontalalignment','center');
b = plot(f,S(:,1),'b'); grid;
r = plot(f,S(:,3),'r');
grid on;
ax1=gca;
set(get(ax1,'Xlabel'),'String', [{'$\bf\ Frequency \hspace{1mm} (Hz) $'}],'interpreter','latex','FontSize',8,'FontName','Arial')
ys = get(ax1,'ylim');
hp1 = text(0.02,diff(ys)*.92+min(ys),'Whole deployment','FontSize',12,'FontWeight','bold','horizontalalignment','left') ;% to write the name of the pannel
legend([b r],'HPF acc x axis (surge)','HPF acc z axis (heave)')

% Subplot 2: descent
ax(2)=subplot(312);
[S,f]=speclev(Aw(DES,:),512,fs);
for nn = [1 3]
    [peakloc, peakmag] = peakfinder(S(:,nn),0.1); % make 0.1 smaller if it is not finding the peak you expect
    peakloc(1) = []; peakmag(1) = [];
    smoothS = runmean(S(:,nn),10);
    peakmag = peakmag - smoothS(peakloc);
    [~,peak] = max(peakmag);
    peak = peakloc(peak);
    plot(f(peak),S(peak,nn),'go','markersize',10,'linewidth',2);
    [minf,bf] = min(S(1:peak,nn));
    FRl(FRn) = f(bf);
    hold on
    plot([f(bf) f(bf)],[min(S(:,nn)) minf],'k--','linewidth',2);
    text(f(bf),min(min(S(:,[1 3]))),['f = ' num2str(round(f(bf)*100)/100)],'horizontalalignment','right','color',cs(nn),'verticalalignment',va{nn});
    FR(FRn) = f(peak);
    FRmag(FRn) = S(peak,nn);
    FRn = FRn+1;
end
[~,b] = max(FRmag(1:2));
if b == 2; v1 = 'top'; v2 = 'bottom'; else v1 = 'bottom'; v2 = 'top'; end
text(FR(1),FRmag(3)+2*sign(length(v1)-length(v2)),num2str(round(FR(3)*100)/100),'verticalalignment',v1,'horizontalalignment','center');
text(FR(2),FRmag(4)-2*sign(length(v1)-length(v2)),num2str(round(FR(4)*100)/100),'verticalalignment',v2,'horizontalalignment','center');
b = plot(f,S(:,1),'b'); grid;
r = plot(f,S(:,3),'r');
grid on;
ax1=gca;
set(get(ax1,'Xlabel'),'String', [{'$\bf\ Frequency \hspace{1mm} (Hz) $'}],'interpreter','latex','FontSize',8,'FontName','Arial')
ys = get(ax1,'ylim');
hp1 = text(0.02,diff(ys)*.92+min(ys),'Descents','FontSize',12,'FontWeight','bold','horizontalalignment','left') ;% to write the name of the pannel
legend([b r],'HPF acc x axis (surge)','HPF acc z axis (heave)')

% Subplot 3: ascent
ax(3)=subplot(313);
[S,f]=speclev(Aw(ASC,:),512,fs);
for nn = [1 3]
[peakloc, peakmag] = peakfinder(S(:,nn),0.1);
peakloc(1) = []; peakmag(1) = [];
smoothS = runmean(S(:,nn),10);
peakmag = peakmag - smoothS(peakloc);
[~,peak] = max(peakmag);
peak = peakloc(peak);
plot(f(peak),S(peak,nn),'go','markersize',10,'linewidth',2)
[minf,bf] = min(S(1:peak,nn));
FRl(FRn) = f(bf);
hold on
plot([f(bf) f(bf)],[min(S(:,nn)) minf],'k--','linewidth',2);
text(f(bf),min(min(S(:,[1 3]))),['f = ' num2str(round(f(bf)*100)/100)],'horizontalalignment','right','color',cs(nn),'verticalalignment',va{nn});
FR(FRn) = f(peak);
FRmag(FRn) = S(peak,nn);
FRn = FRn+1;
end
[~,b] = max(FRmag(3:4));
if b == 2; v1 = 'top'; v2 = 'bottom'; else v1 = 'bottom'; v2 = 'top'; end
text(FR(3),FRmag(5)+2*sign(length(v1)-length(v2)),num2str(round(FR(5)*100)/100),'verticalalignment',v1,'horizontalalignment','center');   
text(FR(4),FRmag(6)-2*sign(length(v1)-length(v2)),num2str(round(FR(6)*100)/100),'verticalalignment',v2,'horizontalalignment','center'); 
b = plot(f,S(:,1),'b'); grid;
r = plot(f,S(:,3),'r');
grid on;
ax2=gca;
set(get(ax2,'Xlabel'),'String', [{'$\bf\ Frequency \hspace{1mm} (Hz) $'}],'interpreter','latex','FontSize',8,'FontName','Arial')
ys = get(ax2,'ylim');
hp1 = text(0.02,diff(ys)*.92+min(ys),'Ascents','FontSize',12,'FontWeight','bold','horizontalalignment','left') ;% to write the name of the pannel
legend([b r],'HPF acc x axis (surge)','HPF acc z axis (heave)')
linkaxes(ax, 'x'); % links x axes of the subplots for zoom/pan


% f = number that multiplied by the FR gives the cut-off frequency fl, of 
% the low pass filter. f is a fraction of FR. You can set default value to
% 0.4 if not, otherwise set f as fl (frequency at the negative peak in the
% power spectral density plot)/FR.
FR = mean(FR(1:2)); 
%f = mean(FRl)/FR;  % from last graph
f=0.23/FR;          % ES - general comment added: Miller et al (2016) cut-off was between 0.19-0.25Hz
alpha=25;
n=1;
k=1:length(p);
J = 0.1;            % dummy variable
tmax = 1/FR;        % dummy variable

%3.2.3_Separate low and high pass filtered acceleration signal using the
%parameters defined earlier and the function Ahf_Alnf.
[Anlf,Ahf,GL,KK] = Ahf_Anlf(Aw,fs,FR,f,n,k,[],[]) ;

%3.2.4_ calculate the smooth pitch from the low pass filter acceleration signal
% to avoid incorporating signals above the stroking periods
[smoothpitch,smoothroll]=a2pr(Anlf(k,:));
%check the difference between pitch and smoothpitch
figure (2); clf;
ax1 = subplot(2,1,1);
plott(p,fs)
ax2 = subplot(2,1,2);
plot((1:length(p))/fs,pitch*180/pi)
hold on
plot((1:length(p))/fs,smoothpitch*180/pi,'r')
legend('pitch','smoothpitch')
linkaxes([ax1 ax2], 'x');


%% 4_DEFINE PRECISE DESCENT AND ASCENT PHASES 
Bottom=[];
Phase(1:length(p))=NaN;

isdive = false(size(p));
for i = 1:size(T,1); isdive(round(T(i,1)*fs):round(T(i,2)*fs)) = true; end
pd = p; pnd = p; pd(~isdive) = nan; pnd(isdive) = nan;
I = 1:length(p);
figure(4); clf; 
sp1 = subplot(5,1,1:4);
plot(I,pd,'g'); hold on; plot(I,pnd,'r'); set(gca,'ydir','rev','ylim',[min(p) max(p)]);
tit = title('Click within the last dive you want to use');
legend('Dive','Not Dive','location','southeast');
sp2=subplot(5,1,5);
plot(I,pitch,'g',I,smoothpitch,'r');
legend('pitch','smoothpitch')
x = ginput(1);
nnn = find(T(:,1)<x(1)/fs,1,'last');  % nn is dive number of the last dive you want to include % ES - changed the name of variable from "nn" to "nnn", as another variable already named nn


% ES - to choose dive phases manually must:
% 1. Comment out the initial enddes and startasc lines (lines 248 and 249) that use smooth pitch to select points where pitch first and last crosses 0
% 2. Uncomment lines from "figure" to "asc=round(x(2))/fs+T(dive,1);"
% 3. Create new "startasc" and "enddes" (lines 259 and 260), which uses
% lines 257 and 258 and multiplies by fs to get in correct time frame

for dive=1:nnn                                  % ES - Line added so as to only include until dive where clicked previously on dive profile    
    kk=round(fs*T(dive,1)):round(fs*T(dive,2)); % it is selecting the whole dive
    kkI = false(size(p)); kkI(kk) = true;
    enddes=round(find((smoothpitch(kkI&p>mindivedef*.75)*180/pi)>0,1,'first')+T(dive,1)*fs);% search for the first point after diving below mindivedef at which pitch is positive
    startasc=round(find((smoothpitch(kkI&p>mindivedef*.75)*180/pi)<0,1,'last')+T(dive,1)*fs);%search for the last point beroe diving above mindivedef at which the pitch is negative
    %if you want to do it manually as some times there is a small ascent
    %phase during the descent and a small descent phase during the ascent.
    %     figure
    %     ax(1)=subplot(211); plott(p(kk),fs)	% plott plots sensor data against a time axis
    %     ax(2)=subplot(212); plott(pitch(kk)*180/pi,fs,0)
    %     linkaxes(ax, 'x'); % links x axes of the subplots for zoom/pan
    %     [x,y]=ginput(2);	% click on where the pitch angle first goes to zero in the descent and last goes to zero in the ascent
%     %     des=round(x(1))/fs+T(dive,1);       % ES - Original code
%     %     asc=round(x(2))/fs+T(dive,1);       % ES - Original code
    enddes=round(x(1))*fs+T(dive,1)*fs;         % ES - edited to correct the two lines above
    startasc=round(x(2))*fs+T(dive,1)*fs;       % ES - changed output names and multiples by fs to get correct time frame
    Phase(kk(kk<enddes))=-1;
    Phase(kk(kk<startasc & kk>enddes)) = 0 ;
    Phase(kk(kk>startasc)) = 1 ;
    Bottom(dive,1)=(enddes)/fs; %Time in seconds at the start of bottom phase (end of descent)
    Bottom(dive,2)= p(enddes);% Depth in m at the start of the bottom phase (end of descent phase)
    Bottom(dive,3)=(startasc)/fs;%Time in seconds at the end of bottom phase (start of descent)
    Bottom(dive,4)= p(startasc);% Depth in m at the end of the bottom phase (start of descent phase)
end
pasc = p; pdes = p;
pasc(Phase<1 | isnan(Phase)) = nan;
pdes(Phase>-1 | isnan(Phase)) = nan;
plot(sp1,pasc,'k'); plot(sp1,pdes,'b');
legend('Dive','Not Dive','Ascent','Descent','location','southeast');
delete(tit)
linkaxes([sp1 sp2], 'x');


%% 5_ESTIMATE SWIM SPEED
thdeg = 30; %degree threshold above which speed can be estimated
[SwimSp] = inst_speed(p,smoothpitch,fs,FR,f,k,thdeg);
%help inst_speed
figure(4);
sp2= subplot(5,1,5);
plot(k,SwimSp,'g'); ylim([0 max(SwimSp)]);
legend('speed')
linkaxes([sp1 sp2], 'x'); % links x axes of the subplots for zoom/pan


%% 6_ESTIMATE SEAWATER DESNSITY AROUND THE TAGGED WHALES

% ES - edited script here to include pressure (i.e. depth) in the SWdensity estimation.
% Involved changing sw_dens0(sali,temp) to sw_dens(sali,temp,P). 
% 1st step required calculating pressure from depth (because pressure, not depth, is input into sw_dens fucntion)
% Downloaded gsw_matlab_v3_05_7 tools from: http://www.teos-10.org/software.htm
% Used gsw_p_from_z function (edited to gsw_p_from_z_ES) to get pressure from depth
% Downloaded Mixing Oceanographic toolbbox from: https://uk.mathworks.com/matlabcentral/fileexchange/47595-mixing--mx--oceanographic-toolbox-for-em-apex-float-data/content/mixing_library/private1/seawater/sw_dens.m
% Used sw_dens tool within SWdesnityfromCTD_ES to calculate SW density with pressure (i.e. depth) considered
% lat = 70.82611; latitude of CTD taken from Miller et al 2016 and  converted to decimal degrees needed for gsw_matlab_v3_05_7

% Home > Import > CTD data as matrix
DTS = CTDJanMayen2013;  % ES - insert correct matrix name; make sure downcast only

DPT=DTS(:,1);           % ES - make sure correct columns called
TMP=DTS(:,3);
SL=DTS(:,2);

[SWdensity,depCTD,temp]=SWdensityFromCTD_ES(DPT, TMP, SL,D);        % ES - edited function to make temp an output and to include depth (pressure) in seawater density calculation
Dsw=EstimateDsw(SWdensity, depCTD, p);
tempw = Estimate_temp_ES(temp, depCTD, p);                          % ES - added line and created function to calculate temperature around the whale


%% 7_EXTRACT STROKES AND GLIDES
% Can be done using the body rotations (pry) estimated using the
% magnetometer method, or it can be done using the dorso-ventral axis of the
% high-pass filtered acceleration signal 
% Using whichever method, tmax and J need to be determined. 
% tmax is the maximum duration allowable for a fluke stroke in seconds, it can be set as 1/FR
% J is the magnitude threshold for detecting a fluke stroke in m /s^2

% 7.1 Using the heave high pass filtered acceleration signal,(use n=3)
% [Anlf,Ahf,GL,KK] = Ahf_Anlf(Aw,fs,FR,f,n,k,J,tmax)
% units of J are in m/s2; set J and tmax [] until determined.
n=3;                % ES - added this line (n otherwise = 1)
[Anlf,Ahf,GL,KK] = Ahf_Anlf(Aw,fs,FR,f,n,k,[],[]);

figure (5); clf;
sp1 = subplot(1,4,1:3);
ax2 = plotyy(1:length(p),p,1:length(pitch),Ahf(:,3));
set(ax2(1),'ydir','rev','ylim',[0 max(p)]);
set(ax2(2),'nextplot','add');
plot(ax2(2),1:length(pitch),Ahf(:,1),'m');
maxy = max(max(Ahf(round(T(1,1)*fs):round(T(nn,2)*fs),[1 3])));
set(ax2(2),'ylim',[-2*maxy 2*maxy],'ytick',round(-2*maxy*10)/10:0.1:2*maxy)
legend('Depth','HPF acc z axis','HPF acc x axis');
title('Zoom in to find appropriate thresholds for fluking, then enter it for J');
linkaxes([ax1 ax2], 'x');

flukeAcc1A = Ahf(ASC,1); % hpf-x acc ascent 
flukeAcc1D = Ahf(DES,1); % hpf-x acc descent 
flukeAcc3A = Ahf(ASC,3); % hpf-z acc ascent
flukeAcc3D = Ahf(DES,3); % hpf-z acc descent

sp2 = subplot(2,4,4);
TOTAL = abs([flukeAcc1A; flukeAcc1D]);  % to look at heave, change to 3
Y=buffer(TOTAL,2*round(1/FR*fs));
hist(max(Y),100); set(sp2,'xlim',[0 max(max(Y))]);
title('hpf-x acc');

sp3 = subplot(2,4,8);
TOTAL = abs([flukeAcc3A; flukeAcc3D]);
Y=buffer(TOTAL,2*round(1/FR*fs));
hist(max(Y),100,'FaceColor','g'); set(sp3,'xlim',[0 max(max(Y))]);
title('hpf-z acc');
% Choose a value for J based on the histogram for:
%   hpf-x, then when detecting glides in the next step use Ahf_Anlf
%   function with n=1
%   hpf-z then when detecting glides in the next step use Ahf_Anlf
%   funct ion with n=3
J= 0.27;    % in m/s2  
tmax=1/FR;  % in seconds

%[Anlf,Ahf,GL,KK] = Ahf_Anlf(Aw,fs,FR,f,n,k,J,tmax);

%in case you want to check the performance of both methods in the following
%figure define glides and strokes obtained from the high pass filtered
%acceleration as GLa and KKa respectively
n=3;
[Anlf,Ahf,GLa,KKa] = Ahf_Anlf(Aw,fs,FR,f,n,k,J,tmax); 

% 7.2 Using the body rotations (pry) estimated using the magnetometer method (use n=1)
figure (6); clf;
sp1 = subplot(1,4,1:3);
n=1;                    % ES - added this to ensure correct n
[MagAcc,pry,Sa,GL,KK] = magnet_rot_sa(Aw,Mw,fs,FR,f,alpha,1,k,J,tmax);
ax1 = plotyy(1:length(p),p,1:length(pry(:,1)),pry(:,1));
set(ax1(1),'ydir','rev','ylim',[0 max(p)]);
maxy = max(pry(round(T(1,1)*fs):round(T(nn,2)*fs)));
set(ax1(2),'ylim',[-2*maxy 2*maxy],'ytick',round(-2*maxy*10)/10:0.1:2*maxy)
set(ax1(2),'nextplot','add');
legend('Depth','rotations in y axis');
title('Zoom in to find appropriate thresholds for fluking, then enter it for J');
sp2 = subplot(1,4,4);
flukePA = pry(ASC,1); % body rotations ascent in radians
flukePD =pry(DES,1);
TOTAL = abs([flukePA; flukePD]);
Y=buffer(TOTAL,2*round(1/FR*fs));
hist(max(Y),100); set(sp2,'xlim',[0 max(max(Y))]);
% Choose a value for J based on the histogram for:
%   pry(:,1), then when detecting glides in the next step use magnet_rot_sa
%   function with n=1

J = 0.4;% in radians 
tmax=1/FR;%in seconds
[MagAcc,pry,Sa,GL,KK] = magnet_rot_sa(Aw,Mw,fs,FR,f,alpha,1,k,J,tmax);

%KK = matrix of cues to zero crossings in seconds (1st column) and 
%zero-crossing directions (2nd column). +1 means a positive-going 
%zero-crossing. Times are in seconds.
%this is already ensuring that all glides are longer than tmax/2

%Check glides duration and positive and negative zero crossings (KK) based 
% on selected J and tmax
figure (7); clf;
sp1 = subplot(5,1,1:4);
plot(k,pd,'g'); hold on; plot(k,pnd,'r','linewidth',2); set(gca,'ydir','rev','ylim',[min(p) max(p)]);
legend('Dive','Not Dive','location','southeast');
plot(sp1,pasc,'k','linewidth',2); plot(sp1,pdes,'b','linewidth',2);
sp2= subplot(5,1,5);
ax2= plotyy(k,SwimSp,k,pry(:,1));
set(ax2(2),'nextplot','add')
%plot(ax2(2),KK(:,1)*fs,pry(round(KK(:,1)*fs),1),'r*')
plot(ax2(2),KKa(:,1)*fs,0,'r*')         % ES - changed from KK to KKa when using accel. method
set(ax2(1),'ylim',[0 max(SwimSp)]);
set(ax2(2),'ylim',[min(pry([ASC DES])) max(pry([ASC DES]))]);
legend('speed','body rotations')
linkaxes([sp1 sp2 ax2], 'x'); % links x axes 
GLdura= GLa(:,2)-GLa (:,1);
GLTa=[GLa(:,1),GLdura];
GLdur= GL(:,2)-GL (:,1);
GLT=[GL(:,1),GLdur]; 
t=(0:length(pry(:,1))-1)/fs;
glk=eventon(GLT,t);
glka=eventon(GLTa,t);
pgl=p;
pgla=p;
pgl(glk==0)=NaN;
pgla(glka==0)=NaN;
h3=plot(sp1,t*fs,pgl,'m:','Linewidth',3);% glides detected with the body rotations (pry)
hold on
h4=plot(sp1,t*fs,pgla,'y:','Linewidth',3);% glides detected with the hpf acc
legend(sp1,'Bottom','Not Dive','Ascent','Descent','Glide','location','southeast');

%SGtype indicates whether it is stroking (1) or gliding(0)
glk1(glka==0)=1;                    % ES - changed from glk to glka when using accel. method
glk1(glka<0)=0;                     % ES - changed from glk to glka when using accel. method
SGtype=glk1;


%% 8 MAKE 5SEC SUB-GLIDES
dur=5;
SGL=splitGL(dur,GLa);               % ES - changed from GL to GLa when using accel. method
SGLT=[SGL(:,1),SGL(:,2)-SGL(:,1)];

%%check that all subglides have a duration of 5 seconds
% rangeSGLT=[min(SGLT(:,2)),max(SGLT(:,2))];
% delete(h3);
% figure(7);
% glk=eventon(SGLT,t);
% pgl=p;
% pgl(glk==0)=NaN;
% h3=plot(sp1,t*fs,pgl,'m','Linewidth',3);


%% 9 Create summary table required for body density 

% FR=0.4;       % ES - commented this as already calculated
% f=0.4;        % ES - commented this as already calculated
alpha=25;
n=1;
k=1:length(p);
J=2/180*pi;
tmax=1/FR;

[MagAcc,pry,Sa,GL,KK] = magnet_rot_sa(Aw,Mw,fs,FR,f,alpha,n,k,J,tmax);
[smoothhead] = m2h(MagAcc.Mnlf(k,:),smoothpitch,smoothroll);
smoothpitchdeg=smoothpitch*180/pi;
smoothheaddeg=smoothhead*180/pi;  % ES - added this line to convert from radians to degrees
smoothrolldeg=smoothroll*180/pi;  % ES - added this line to convert from radians to degrees

% ES - commented below is original Glide variable creation 
% Glide = zeros(length(SGL),24) ;
% for i=1:length(SGL)
%             cue1=SGL(i,1)*fs;
%             cue2=SGL(i,2)*fs;
%     Glide(i,1)=SGL(i,1);%sub-glide start point in seconds
%     Glide(i,2)=SGL(i,2);%sub-glide end point in seconds
%     Glide(i,3)=SGL(i,2)-SGL(i,1);%sub-glide duration
%     Glide(i,4)=mean(p(round(cue1):round(cue2)));%mean depth(m)during sub-glide
%     Glide(i,5)=abs(p(round(cue1))-(p (round(cue2))));%total depth(m)change during sub-glide
%     Glide(i,6)=mean(SwimSp(round(cue1):round(cue2)));%mean swim speed during the sub-glide, only given if pitch>30 degrees
%     Glide(i,7)=mean((smoothpitchdeg(round(cue1):round(cue2))));%mean pitch during the sub-glide
%     Glide(i,8)=sin(mean(smoothpitch(round(cue1):round(cue2))));%mean sin pitch during the sub-glide 
%     Glide(i,9)=std(smoothpitch(round(cue1):round(cue2)));%SD of pitch during the sub-glide  
%     Glide(i,10)=mean(tempw(round(cue1):round(cue2)));%mean temperature during the sub-glide 
%     Glide(i,11)=mean(Dsw(round(cue1):round(cue2)));%mean seawater density (kg/m^3) during the sub-glide
%     try  
%     xpoly=(round(cue1):round(cue2))';
%         ypoly=SwimSp(round(cue1):round(cue2));
%         [B,BINT,R,RINT,STATS] = regress(ypoly,[xpoly ones(length(ypoly),1)]);
%    
%     Glide(i,12)=B(1);%mean acceleration during the sub-glide
%     Glide(i,13)=STATS(1);%R2-value for the regression swim speed vs. time during the sub-glide
%     Glide(i,14)=STATS(4);%SE of the gradient for the regression swim speed vs. time during the sub-glide
%     catch
%     Glide(i,12)=NaN;%mean acceleration during the sub-glide
%     Glide(i,13)=NaN;%R2-value for the regression swim speed vs. time during the sub-glide
%     Glide(i,14)=NaN;%SE of the gradient for the regression swim speed vs. time during the sub-glide
%     end
%     sumphase=sum(Phase(round(cue1):round(cue2)));  
%         sp=NaN;
%         sp(sumphase<0)=-1;
%         sp(sumphase==0)= 0 ;
%         sp(sumphase>0)=1 ;
%     Glide(i,15)=sp;%Dive phase:0 bottom, -1 descent, 1 ascent, NaN not dive phase       
%         Dinf=D(find((D(:,1)*fs)<cue1 &(D(:,2)*fs)>cue2),:);
%         if isempty(Dinf),
%           Dinf=NaN(size(D,1),size(D,2)) ;
%         end
%     Glide(i,16)= Dinf(7);%Dive number in which the sub-glide recorded
%     Glide(i,17)=Dinf(6);%Maximum dive depth (m) of the dive
%     Glide(i,18)=Dinf(3);%Dive duration (s) of the dive
%     Glide(i,19)=circ_mean(smoothpitchdeg(round(cue1):round(cue2)));%Mean pitch(deg) calculated using circular statistics 
%     Glide(i,20)=1-(circ_var(smoothpitch(round(cue1):round(cue2))));%Measure of concentration (r) of pitch during the sub-glide (i.e. 0 for random direction, 1 for unidirectional)
%     Glide(i,21)=circ_mean(smoothrolldeg(round(cue1):round(cue2)));%Mean roll (deg) calculated using circular statistics
%     Glide(i,22)=1-(circ_var(smoothroll(round(cue1):round(cue2))));%Measure of concentration (r) of roll during the sub-glide
%     Glide(i,23)=circ_mean(smoothheaddeg(round(cue1):round(cue2)));%Mean heading (deg) calculated using circular statistics 
%     Glide(i,24)=1-(circ_var(smoothhead(round(cue1):round(cue2))));%Measure of concentration (r) of heading during the sub-glide
% end

% ES - new Glide creation code, including:
% 1. Fix for acceleration values (i.e. in m/s2 not m/s/#samples/s) 
% 2. Fix to ensure acceleration calculated only if speed values available throughout glide (as in Miller 2016 Igor Pro procedures)
% 3. Fix for correct units of pitch, roll and heading
Glide = zeros(length(SGL),24) ;
for i=1:length(SGL)
            cue1=SGL(i,1)*fs;
            cue2=SGL(i,2)*fs;
    Glide(i,1)=SGL(i,1);%sub-Glide start point in seconds
    Glide(i,2)=SGL(i,2);%sub-Glide end point in seconds
    Glide(i,3)=SGL(i,2)-SGL(i,1);%sub-Glide duration
    Glide(i,4)=mean(p(round(cue1):round(cue2)));%mean depth(m)during sub-Glide
    Glide(i,5)=abs(p(round(cue1))-(p (round(cue2))));%total depth(m)change during sub-Glide
    Glide(i,6)=mean(SwimSp(round(cue1):round(cue2)));%mean swim speed during the sub-Glide, only given if pitch>30 degrees
    Glide(i,7)=mean((smoothpitchdeg(round(cue1):round(cue2))));%mean pitch during the sub-Glide
    Glide(i,8)=sin(mean(smoothpitch(round(cue1):round(cue2))));     % ES - changed this from "smoothpitchdeg" to "smoothpitch",i.e. changed from degrees to radians %mean sin pitch during the sub-Glide
    Glide(i,9)=std(smoothpitch(round(cue1):round(cue2)));           % ES - changed this from "smoothpitchdeg" to "smoothpitch",i.e. changed from degrees to radians  %SD of pitch during the sub-Glide
    Glide(i,10)=mean(tempw(round(cue1):round(cue2)));               % ES - changed this from temp to tempw so refers to temperature around the whale not CTD temp %mean temperature during the sub-Glide 
    Glide(i,11)=mean(Dsw(round(cue1):round(cue2)));%mean seawater density (kg/m^3) during the sub-Glide
    
    if all(isnan(SwimSp(round(cue1):round(cue2))))==0;              % ES - added this line, if all speed values within the glide are not NaN
    if isnan(Glide(i,6)) ==0                                        % ES - added this line, so that only calcuate accelertion if there is a mean swim speed (otherwise accel given as 0 and following errors occur: "> In regress (line 84) Warning: X is rank deficient to within machine precision."
    try  
    xpoly=(round(cue1):round(cue2))';
    % xpoly=(SGL(i,1):0.2:SGL(i,2))'; 
        ypoly=SwimSp(round(cue1):round(cue2));
        [B,BINT,R,RINT,STATS] = regress(ypoly,[xpoly ones(length(ypoly),1)]);
   
    Glide(i,12)=B(1)*fs;    % ES - changed this to ensure accel in m/s2, not m/fs %mean acceleration during the sub-Glide
    Glide(i,13)=STATS(1);   % ES - general note: this does not change with different gradients and SE %R2-value for the regression swim speed vs. time during the sub-Glide
    Glide(i,14)=STATS(4)*fs;% ES - changed this to correct for slope in s, in line with new gradient in ms2 %SE of the gradient for the regression swim speed vs. time during the sub-Glide; 
    catch
    Glide(i,12)=NaN;%mean acceleration during the sub-Glide
    Glide(i,13)=NaN;%R2-value for the regression swim speed vs. time during the sub-Glide
    Glide(i,14)=NaN;%SE of the gradient for the regression swim speed vs. time during the sub-Glide
    end
    else    % ES - added these lines (to next end) to complete the if statement re. mean speeds not NaN
    Glide(i,12)=NaN;%mean acceleration during the sub-Glide
    Glide(i,13)=NaN;%R2-value for the regression swim speed vs. time during the sub-Glide
    Glide(i,14)=NaN;%SE of the gradient for the regression swim speed vs. time during the sub-Glide
    end
    else    % ES - added these lines (to next end) to complete the if statement re. SwimSp nan
    Glide(i,12)=NaN;%mean acceleration during the sub-Glide
    Glide(i,13)=NaN;%R2-value for the regression swim speed vs. time during the sub-Glide
    Glide(i,14)=NaN;%SE of the gradient for the regression swim speed vs. time during the sub-Glide
    end     % ES - added this to complete if statement
    sumphase=sum(Phase(round(cue1):round(cue2)));  
        sp=NaN;
        sp(sumphase<0)=-1;
        sp(sumphase==0)= 0 ;
        sp(sumphase>0)=1 ;
    Glide(i,15)=sp;%Dive phase:0 bottom, -1 descent, 1 ascent, NaN not dive phase       
        Dinf=D(find((D(:,1)*fs)<cue1 &(D(:,2)*fs)>cue2),:);
        if isempty(Dinf),
          Dinf=NaN(size(D,1),size(D,2)) ;
        end
    Glide(i,16)= Dinf(7);%Dive number in which the sub-Glide recorded
    Glide(i,17)=Dinf(6);%Maximum dive depth (m) of the dive
    Glide(i,18)=Dinf(3);%Dive duration (s) of the dive
    Glide(i,19)=circ_mean(smoothpitchdeg(round(cue1):round(cue2)));     % ES - changed this from "smoothpitch" to "smoothpitchdeg" %Mean pitch(deg) calculated using circular statistics  
    Glide(i,20)=1-(circ_var(smoothpitch(round(cue1):round(cue2))));         %Measure of concentration (r) of pitch during the sub-Glide (i.e. 0 for random direction, 1 for unidirectional)
    Glide(i,21)=circ_mean(smoothrolldeg(round(cue1):round(cue2)));      % ES - changed this from "smoothroll" to "smoothrolldeg" %Mean roll (deg) calculated using circular statistics  
    Glide(i,22)=1-(circ_var(smoothroll(round(cue1):round(cue2))));          %Measure of concentration (r) of roll during the sub-Glide
    Glide(i,23)=circ_mean(smoothheaddeg(round(cue1):round(cue2)));      % ES - changed this from "smoothhead" to "smoothheaddeg" %Mean heading (deg) calculated using circular statistics
    Glide(i,24)=1-(circ_var(smoothhead(round(cue1):round(cue2))));          %Measure of concentration (r) of heading during the sub-Glide
end


%% 10. Calculate glide ratio
G_ratio = zeros(nnn,10);            % ES - changed this from "G_ratio = zeros(size(T,1),10) ;" so as not to include last dive
for dive=1:nnn                      % ES - changed this from "for dive=1:size(T,1)" so as not to include last dive
    
    kkdes=round(fs*T(dive,1)):round(fs*Bottom(dive,1)); % it is selecting the whole dive
    kkas=round(fs*Bottom(dive,3)):round(fs*T(dive,2)); % it is selecting the whole dive
    G_ratio(dive,1)=length(kkdes)/fs;%total duration of the descet phase (s)
    G_ratio(dive,2)=length(find(SGtype(kkdes)==0))/fs;%total glide duration during the descet phase (s)
    G_ratio(dive,3)=G_ratio(dive,2)/G_ratio(dive,1);%glide ratio during the descet phase
    G_ratio(dive,4)=mean(smoothpitchdeg(kkdes)*180/pi);%mean pitch during the descet phase(degrees)     % ES - changed this from "smoothpitch" to "smoothpitchdeg"
    G_ratio(dive,5)=Bottom(dive,2)/G_ratio(dive,1);%descent rate (m/s)
    G_ratio(dive,6)=length(kkas)/fs;%total duration of the ascet phase (s)
    G_ratio(dive,7)=length(find(SGtype(kkas)==0))/fs;%total glide duration during the ascet phase (s)
    G_ratio(dive,8)=G_ratio(dive,7)/G_ratio(dive,6);%glide ratio during the ascet phasennn
    G_ratio(dive,9)=mean(smoothpitchdeg(kkas)*180/pi);%mean pitch during the ascet phase(degrees)       % ES - changed this from "smoothpitch" to "smoothpitchdeg"
    G_ratio(dive,10)=Bottom(dive,3)/G_ratio(dive,6);%ascent rate (m/s)
end



%% 11. Export csvs
% ES - Changed from csvwrtie to dlmwrite. Do not use "csvwrite" as this reduces precision of results
%%%%%%%%%%% INSERT CORECT TAG NAME
dlmwrite('Y:\Eilidh\PhD\Chapter3\Part_1_glide_extraction\Glide_summaries\Run10\TAGNAME.csv',Glide,'precision','%.7f')
dlmwrite('Y:\Eilidh\PhD\Chapter3\Part_1_glide_extraction\MatLab\BodyDensity_MATLAB\glide_ratio_tables\Run10\TAGNAME.csv',G_ratio,'precision','%.7f')
