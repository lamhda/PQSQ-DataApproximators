function c = PQSQ_Mean(x, potentialFunction, eps)
%PQSQ_Mean calculates PQSQ central point: point which provide the minimum
%of PQSQ error sum for x-c.
%   x is n-by-m matrix of data
%   potentialFunction is structure which define PQSQ potential function
%       (see definePotentialFunction.m for details)
%   eps is tolerance level for iterations convergence
%
%Examples
%   Generate data
%   x = rand(100,5);
%   Compare different centroids:
%   cMean = mean(x)
%   cMedian  = median(x)
%   Create potential function for L1
%   potFun = definePotentialFunction( x, 5, @L1 );
%   cPQSQL1 = PQSQ_Mean(x, potFun);
%   Create potential function for L2
%   potFun = definePotentialFunction( x, 5, @L2 );
%   cPQSQL2 = PQSQ_Mean(x, potFun);
%   Create potential function for L1.5
%   potFun = definePotentialFunction( x, 5, @L1_5 );
%   cPQSQL1_5 = PQSQ_Mean(x, potFun);
%
%   Result
%cMean =
%    0.4928    0.4749    0.5234    0.4981    0.4944
%cMedian =
%    0.4702    0.4514    0.5490    0.4680    0.5150
%cPQSQL1 =
%    0.4740    0.4576    0.5597    0.4775    0.5228
%cPQSQL2 =
%    0.4928    0.4749    0.5234    0.4981    0.4944
%cPQSQL1_5 =
%    0.4893    0.4710    0.5305    0.4919    0.4972
% we can see that cPQSQL2 is the same as cMean. It is expected result.
% cPQSQL1 is the cloasest for the cMedian. 
% cPQSQL1_5 is intermediate.


    if nargin<3
        eps = 0.001;
    end

    m = size(x, 2);
    
    %Calculate results for three 'typical' initial points
    center_zero = PQSQ_Mean_Prime(x,potentialFunction,zeros(1,m),eps);
    center_mean = PQSQ_Mean_Prime(x,potentialFunction,mean(x,1),eps);
    center_median = PQSQ_Mean_Prime(x,potentialFunction,median(x,1),eps);

    %Calculate norms for comparison
    norm = [PQSQ_Norm(bsxfun(@minus,x,center_zero), potentialFunction,1);...
        PQSQ_Norm(bsxfun(@minus,x,center_mean), potentialFunction,1);...
        PQSQ_Norm(bsxfun(@minus,x,center_median), potentialFunction,1)];
    
    %Select the minimal values
    [~,ind] = min(norm);
    
    %Form output vector
    c = center_zero;
    ii = ind==2;
    c(ii) = center_mean(ii);
    ii = ind==3;
    c(ii) = center_median(ii);
end

