% getDiracMoments_fourier( x0, N)
% Returns the fourier moments of a dirac measure support on a point x0 \in
% [0,1], i.e., y(k) = exp(i*2*pi*x0*k), k=0,1,..N
% getDiracMoments_fourier( x0, N, w) returns the moments of the weighted
% sum of diracs at points x0(1)..x0(m)

function y = getDiracMoments_fourier( x0, N, w )
    if(~exist('w','var') || isempty(w))
        w = ones(numel(x0),1);
    end
    w = w(:);
    y = exp(1i*2*pi*[0:N]'*x0)*w;
end

