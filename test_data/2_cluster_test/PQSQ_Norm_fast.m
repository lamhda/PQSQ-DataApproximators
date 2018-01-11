function norm = PQSQ_Norm_fast( x, potentialFunction, ret )
%Calculate the PQSQ Error for one or more points vector.
%  x is n-by-m matrix with rows which correspond to data points.
%  intervals is m-by-K matrix with K values of thresholds which corresponds
%       to each coordinate. intervals(:,K) MUST be equal to Inf.
%   A and B are matrices of coefficients of quadratic functions
%   ret is type of returned values:
%       0 or omitted means returning the whole matrix of potentials
%       1 means returning the row vector of sums of potentials for each
%           coordinate
%       2 means returning the column vector of sums of potentials for each
%           data point.

    [n, m] = size(x);
    norm = zeros(n,m);

    for k=1:m
        norm(:,k) = PQSQ_fast(x(:,k),potentialFunction,k);
    end

    if nargin>2 
        if ret>0
            norm = sum(norm,ret);
        end
    end
end

