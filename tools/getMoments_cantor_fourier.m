% Returns fourier moments exp(2*pi*i*k), k=0,...,N of the cantor measure on [0,1]

function y = getMoments_cantor_fourier( N )
t = 2*pi*[0:N];
M = [0:10000]';
size(cos( (1./(3.^M)) * t ))
y = exp(1i*t/2).* prod(cos( (1./(3.^M)) * t ));
y = y';
end

