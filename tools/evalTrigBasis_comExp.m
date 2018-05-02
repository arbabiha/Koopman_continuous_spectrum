% Evalueates exp(1i*2*pi*x(j)*k), k=0...N

% INPUT:
% x - 1 x K vector of points in R
% N - degree

%OUTPUT
% N x K matrix V

function  V = evalTrigBasis_comExp( x,N, doubleSided )
if(~exist('doubleSided','var'))
    doubleSided = 0;
end
x = x(:); x = x';
if(~doubleSided)
    V = exp(1i*2*pi*[0:N]'*x);
else
    V = exp(1i*2*pi*[-N:N]'*x); 
end

end

