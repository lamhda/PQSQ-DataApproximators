function c = PQSQR_Mean(x, potentialFunction, eps, c0)
%PQSQR_Mean calculates PQSQ central point: point which provide the minimum
%of PQSQ error sum for x-c.
%   x is n-by-m matrix of data
%   potentialFunction is structure which define PQSQ potential function
%       (see definePotentialFunction.m for details)
%   eps is tolerance level for iterations convergence (optional)
%   c0 is the 1-by-m row vector with prototype of centre (optional)
%
%Examples
%   Generate data
%   x = rand(100,5);
%   Compare different centroids:
%   cMean = mean(x,1)
%   cMedian  = median(x,1)
%   Create potential function for L1
%   potFun = PQSQR_definePotentialFunction( x, 5, @L1 );
%   cPQSQL1 = PQSQR_Mean(x, potFun)
%   Create potential function for L2
%   potFun = PQSQR_definePotentialFunction( x, 5, @L2 );
%   cPQSQL2 = PQSQR_Mean(x, potFun)
%   Create potential function for L1.5
%   potFun = PQSQR_definePotentialFunction( x, 5, @L1_5 );
%   cPQSQL1_5 = PQSQR_Mean(x, potFun)
%
%   Result
% cMean =
%     0.5344    0.4675    0.5014    0.5115    0.5122
% cMedian =
%     0.5358    0.4079    0.4866    0.5293    0.5452
% cPQSQL1 =
%     0.5462    0.4697    0.4953    0.4984    0.5153
% cPQSQL2 =
%     0.5344    0.4675    0.5014    0.5115    0.5122
% cPQSQL1_5 =
%     0.5391    0.4684    0.4989    0.5063    0.5134
% we can see that cPQSQL2 is the same as cMean. It is expected result.





    if nargin<3
        eps = 0.001;
    end

    m = size(x, 2);
    
    %Calculate results for three 'typical' initial points
    if nargin>3 
        %Use prototype instead of zeros
        center_zero = PQSQR_Mean_Prime(x,potentialFunction,c0,eps);
    else
        center_zero = PQSQR_Mean_Prime(x,potentialFunction,zeros(1,m),eps);
    end
    center_mean = PQSQR_Mean_Prime(x,potentialFunction,mean(x),eps);
    center_median = PQSQR_Mean_Prime(x,potentialFunction,median(x),eps);

    %Calculate norms for comparison
    norm = [sum(PQSQR(bsxfun(@minus,x,center_zero), potentialFunction));...
            sum(PQSQR(bsxfun(@minus,x,center_mean), potentialFunction));...
            sum(PQSQR(bsxfun(@minus,x,center_median), potentialFunction))];
    
    %Select the minimal values
    [~,ind] = min(norm);
    
    %Form output vector
    if ind == 1
        c = center_zero;
    elseif ind == 2
        c = center_mean;
    else
        c = center_median;
    end
end

