function c = PQSQ_Mean_Prime_fast(x, potentialFunction, mn, varargin)
%PQSQ_Mean calculates PQSQ central point: point which provide the minimum
%of PQSQ error sum for x-c.
%   x is n-by-m matrix
%   potentialFunction is structure, created by definePotentialFunction.

    verbose=0;
    eps=0.001;
    
    m = size(x, 2);
    c = mn;

    for i=1:length(varargin)
        if strcmp(varargin{i},'verbose')
            verbose = varargin{i+1};
        end
        if strcmp(varargin{i},'eps')
            eps = varargin{i+1};
        end
    end
    
    for k=1:m
        
        if(verbose)
            display(sprintf('Coordinate %i:',k));
        end
        
        count=0;
        while(count<100)
            m0 = c(k);
            inds = splitPointIntervalsFast1(abs(x(:,k)-c(k)),potentialFunction.intervals(k,:));
            as = potentialFunction.A(k,inds)';
            x2 = sum(as);
            x1 = sum(as.*x(:,k));
            if x2~=0
                c(k)=x1/x2;
            else
                c(k)=0;
            end
            
            count=count+1;
            delta = abs(c(k)-m0)/abs(c(k)+0.001);
            if(delta<eps)
                break;
            end;
        end;
    end;
end