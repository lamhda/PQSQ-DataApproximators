function c = PQSQ_Mean_fast(x, potentialFunction, varargin)
%PQSQ_Mean calculates PQSQ central point: point which provide the minimum
%of PQSQ error sum for x-c.
%   x is n-by-m matrix
%   potentialFunction is structure, created by definePotentialFunction.

    m = size(x, 2);
    
    %Calculate results for three 'typical' initial points
    center_zero = PQSQ_Mean_Prime_fast(x,potentialFunction,zeros(1,m),varargin);
    center_mean = PQSQ_Mean_Prime_fast(x,potentialFunction,mean(x),varargin);
    center_median = PQSQ_Mean_Prime_fast(x,potentialFunction,median(x),varargin);

    %Calculate norms for comparison
    norm = [PQSQ_Norm_fast(bsxfun(@minus,x,center_zero), potentialFunction,1);...
        PQSQ_Norm_fast(bsxfun(@minus,x,center_mean), potentialFunction,1);...
        PQSQ_Norm_fast(bsxfun(@minus,x,center_median), potentialFunction,1)];
    
    %Select the minimal values
    [~,ind] = min(norm);
    
    %Form output vector
    c = center_zero;
    ii = ind==1;
    c(ii) = center_mean(ii);
    ii = ind==2;
    c(ii) = center_median(ii);
end

