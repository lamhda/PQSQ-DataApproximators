function [V,U] = firstPrincipalComponentPQSQ_fast(x, potentialFunction, varargin)
%firstPrincipalComponentPQSQ calculates the first principal component for
%specified data matrix x, intervals and potential_function_handle.
% x is n-by-m matrix with n observations and m attributes.
% intervals is m-by-(p+1) matrix with (p+1) boundaries of p intervals for
%   each dimension of data in row.
% potential_function_handle is handler for energy potential function.
% Name, Value optionas:
% 'verbose', 1 is indicator of detailed output (debuging).
% 'eps', positive number is tolerance level for component calculation.
% 'init', -1/0/positive number is method of component initialization:
%               -1 means uniform random generation in interval 0, stdev for
%                   each coordiate
%               0 means initialization by the first L2 principal component
%                   for data matrix x.
%               positive number is the number of data point to use it as
%                   initialization.
% 'optimize', 
%       0 means never use additional optimization
%       1 means additional optimization at the end of calculations
%       2 means additional optimization at the each iteration and at the
%           end.
%
%Output arguments
%V is m-by-1 vector of found principal component.
%U is n-by-1 vector of length of projections onto found principal
%   component.

    verbose=0;
    eps=0.001;
    optimizeProjections = 1;
    initiatilization = 0; % 0 - PC1, -1 - random vector, >0 - ith data vector
%     meanDefined = 0;

    for i=1:length(varargin)
        if strcmp(varargin{i},'verbose')
            verbose = varargin{i+1};
        elseif strcmp(varargin{i},'eps')
            eps = varargin{i+1};
        elseif strcmp(varargin{i},'init')
            initiatilization = varargin{i+1};
        elseif strcmp(varargin{i},'optimize')
            optimizeProjections = varargin{i+1};
        end
    end

    %define dimensions of vectors and matrices
    [n, m] = size(x);

    RS = zeros(n, m);

    %%%%%%%%%%%%%%%%%%%%%
    % initialize V and U
    %%%%%%%%%%%%%%%%%%%%%
    if initiatilization==-1
        stdev = std(x)';
        V = rand(m,1).*stdev;
    elseif initiatilization==0
        [V,~] = eigs(x'*x,1);
    elseif initiatilization>0
        V=x(initiatilization,:)';
    end

    %Renormalize V
    V = V/norm(V);
    %Initialize U
    U = x*V;
       
    
    count = 1;
    
    XU=U;
    SU=U;
    XV = V;
    SV = V;
    
    while(count<1000)
        
        if optimizeProjections<2
            %calculate the matrix of indices of intervals
            for k=1:m
                inds = splitPointIntervalsFast1(x(:,k)-V(k)*U(:),potentialFunction.intervals(k,:));
                RS(:,k) = inds(:);
            end

            %Initialize arrays XU and SU
            XU(:)=0;
            SU(:)=0;
            for k=1:m
                XU = XU + potentialFunction.A(k,RS(:,k))'.*x(:,k)*V(k);
                SU = SU + potentialFunction.A(k,RS(:,k))'*V(k)*V(k);
            end
            ind = SU==0;
            U(ind) = 0;
            U(~ind) = XU(~ind)./SU(~ind);
        else
            U = optimizePQSQprojections(x,V,potentialFunction);
        end
        
        V1=V;

        %calculate the matrix of indices of intervals
        for k=1:m
            inds = splitPointIntervalsFast1(x(:,k)-V(k)*U(:),potentialFunction.intervals(k,:));
            XV(k) = sum(potentialFunction.A(k,inds)'.*x(:,k).*U(:));
            SV(k) = sum(potentialFunction.A(k,inds)'.*U(:).*U(:));
        end
        ind = SV==0;
        V(ind) = 0;
        V(~ind) = XV(~ind)./SV(~ind);
        
        %renormalize V1
        V = V/norm(V1);
        
        delta = (norm(V-V1));
        
        if(delta<eps)
            break;
        end
        
        count = count+1;
    end
    
    %Standartize direction
    if(V(1)<0)
        V=-V;
        U=-U;
    end
    
    %Final projection normalization
    if optimizeProjections>0
        U = optimizePQSQprojections(x,V,potentialFunction);
    end
end