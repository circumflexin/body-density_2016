
function SGL=splitGL(dur,GL);
%   SGL=splitGL(dur,GL);
%   Make x duration subglides.

%    INPUT:
%       dur = desired duration of glides
%       GL = matrix containing the start time (first column) and end time
%           (2nd column) of any glides.Times are in seconds.
%
%
%    OUTPUT: 
%          SGL = matrix containing the start time (first column) and end 
%           time(2nd column) of the generated subglides.All glides must
%           have duration equal to the given dur value.Times are in seconds.
%
%   Lucia Martina Martin Lopez (May 2016)
%   lmml2@st-andrews.ac.uk
%   Updated Feb 2020 Alec Burslem
%   acb35@st-andrews.ac.uk


SUM=[];
GLINF = [GL,GL(:,2)-GL(:,1)];
for i=1:length(GL)
    sum=[];
    ng = floor(GLINF(i,3)/dur); % AB- floor division drops the partial glide here
    if ng > 0
        STARTGL=[];
        ENDGL=[];
        for k = 1:ng
            startglide = GLINF(i,1)+(k-1)*dur;
            endglide = startglide+dur;
            STARTGL=[STARTGL;startglide];
            ENDGL=[ENDGL;endglide];
        end
        sum=[STARTGL,ENDGL];
    end
    SUM=[SUM;sum];
    SGL=[SUM];
end


        
% Original:
% SUM=[];
% for i=1:length(GL)
%     sum=[];
%     GLINF=[GL,GL(:,2)-GL(:,1)];
%     ng=(GLINF(i,3)/(dur+1));
%     if abs(GLINF(i,3))>dur %If this duration is bigger than 5
%         STARTGL=[];
%         ENDGL=[];
%         for k=1:round(ng)
%             v=((1:max(round(ng)))-1)*6; %should be *dur?
%             startglide1=GLINF(i,1)+v(k);
%             endglide1=startglide1+5;
%             STARTGL=[STARTGL;startglide1];
%             ENDGL=[ENDGL;endglide1];
%         end
%         sum=[STARTGL,ENDGL];
%         
%     end
%     SUM=[SUM;sum];
%     SGL2=[SUM];
% end
% 
% 
% 







