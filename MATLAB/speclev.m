function    [SL,f]=speclev(x,nfft,fs,w,nov)
%
%    [SL,f]=speclev(x,nfft,fs,w,nov) S is the amount of power in each particular frequency (f)
%  x is the signal from which the speclev (power spectra) is going to be
%  calculated 
%fft  is the fourier transform length
%%fs is the frequency sample
%%w is that if there is less than 4 inputs then w is equal to the nfft 
%% nov is that if there is less than 5 inputs then nov which is the overlap will be half of the size of the nfft
%    mark johnson, WHOI
%    mjohnson@whoi.edu
if nargin<2,
   nfft = 512 ;
end

if nargin<3,
   fs = 1 ;
end

if nargin<4 | isempty(w),
   w = nfft ;
end

if nargin<5,
   nov = nfft/2 ; %%%overlap of half of the size of the fft
end

if length(w)==1,
   w = hanning(w) ;
end

P = zeros(nfft/2,size(x,2)) ;
for k=1:size(x,2),
   [X,z] = buffer(x(:,k),length(w),nov,'nodelay') ;
   X = detrend(X).*repmat(w,1,size(X,2)) ;
   F = abs(fft(X,nfft)).^2 ;
   P(:,k) = sum(F(1:nfft/2,:),2) ;
end

ndt = size(X,2) ;

% these two lines give correct output for randn input
% SL of randn should be -10*log10(fs/2)

slc = 3-10*log10(fs/nfft)-10*log10(sum(w.^2)/nfft) ;
SL = 10*log10(P)-10*log10(ndt)-20*log10(nfft)+slc ;
f = (0:nfft/2-1)/nfft*fs ;
