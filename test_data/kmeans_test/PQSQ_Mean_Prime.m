function c = PQSQ_Mean_Prime(x, potentialFunction, mn, eps)
%PQSQ_Mean calculates PQSQ central point: point which provide the minimum
%of PQSQ error sum for x-c.
%   x is n-by-m matrix
%   potentialFunction is structure which define PQSQ potential function
%       (see definePotentialFunction.m for details)
%   eps is tolerance level for iterations convergence

    m = size(x, 2);
    c = mn;

    for k=1:m
        count=0;
        while(count<100)
            m0 = c(k);
            inds = identifyIntervals(abs(x(:,k)-c(k)),potentialFunction.intervals(k,:));
            as = potentialFunction.A(k,inds)';
            x2 = sum(as);
            x1 = sum(as.*x(:,k));
            if x2~=0
                c(k) = x1/x2;
            else
                c(k) = 0;
            end
            
            count = count+1;
            delta = abs(c(k)-m0)/abs(c(k)+0.001);
            if(delta<eps)
                break;
            end;
        end;
    end;
end