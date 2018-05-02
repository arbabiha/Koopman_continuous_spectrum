function phi=bartlett_estimator(x,M,L)
% computes the power spectral density via Bartlett method based on
% M. S. Bartlett, Periodogram analysis and continuous spectra, Biometrika 37, 1 (1950).
%
% x: the signal vector (n x 1) 
% M: length of each subsample    
% L: size of the padded FFT grid

phi=welch_estimator(x,ones(M,1),M,L);   % bartlett is a special case of Welch. :-)
end