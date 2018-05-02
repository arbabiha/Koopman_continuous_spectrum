% y -  vector of Fourier moments of lenghth:
%      2*N+1 if bothSidesProvided = 1 (default) ordered from negative to positive frequencies
%           i.e., y(mid + k) = \int_0^1 exp(2*pi*x*k)\, dmu(x), k = -N...N, where mid = (N+1)/2
%           just accounts from matlab's indexing from 1 and not from -N

%       N+1 if bothSidesProvided = 0, only the nonnegative frquencies, i.e,
%           y(k) = \int_0^1 exp(2*pi*x*k)\, dmu(x), k = 0,..,N
%           In this case, we assume the measure is real and hence the
%           negative moments are just complex conjugates of the positve ones


% Returns the matrix

% |m0 \bar{m1} \bar{m2} ...              |
% |m1 m0       \bar{m1} ...              | 
% |m2  m1       m0      ...              |  


function M = momentMat_fourier( y, bothSidesProvided )
if(mod(numel(y),2) == 0)
    error('Number of fourier moments must be odd');
end
y = y(:);

if(nargin == 1)
    bothSidesProvided = 1;
end
if(bothSidesProvided == 0)
    y = [conj(fliplr(y(2:end)')') ; y];
end

N = (numel(y) - 1) / 2;
mid = (numel(y)+1)/2;
M = zeros(N+1,N+1);
for i = 1:N+1
    for j = 1:N+1
        M(i,j) = y(mid + i-j);
    end
end

end

