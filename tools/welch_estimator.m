function phi=welch_estimator(x,v,K,L)
% computes the power spectral density via Welch method based on
% P. Welch, The use of fast Fourier transform for the estimation of power spectra: A method based on time
% averaging over short, modified periodograms, IEEE Trans. Audio Electroacoust. 15, 70 (1967).
%
% x: the signal vector (n x 1) 
% v: window function    
% K: overlap length      
% L: size of the padded FFT grid


N=length(x);     % total data length
M=length(v);     % length of each subsample


S=floor((N-M+K)/K); % number of subsamples
Pv=mean(v.^2);     % the power of the window vector v

% function handle for computing psd for each subsample
psd_fft = @(y) (abs(fft((y(1:M).*v(:)),L)).^2)/M;

% compute the windowed fft and average over subsamples
phi=zeros(L,1);
for i = 1 : S,
   phi=phi + psd_fft(x((i-1)*K+1:(i-1)*K+M));
end

phi=phi/S/Pv;


end

