% x - 1 x K vector of points to evaulate the Christoffel-Daruboux kernel
% invM - inverse of the moment matrix
% OUTPUT: c - values of the CD kernel at x
function c = evalChristoffelPol( invM, x )
N = size(invM,1)-1;
V = evalTrigBasis_comExp( x, N, 0 );
tmp = invM*V;
c = zeros(1,numel(x));
for i = 1:numel(x)
    c(i) = V(:,i)'*tmp(:,i);
end
end

