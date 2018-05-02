function y = getLebesgueMoments_fourier( N, a, b, doubleSided )

if(~exist('doubleSided','var'))
    doubleSided = 0;
end

k = [1:N]';
y = ( exp(1i*2*pi*k*b) - exp(1i*2*pi*k*a) ) ./ (2*pi*1i*k);
if(~doubleSided)
    y = [b - a; y];
else
    y = [conj(flipud(y)) ;  b - a; y];
end

end

