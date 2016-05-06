function U = optimizePQSQprojections(x, V, potentialFunction, eps)
%optimizePQSQprojections calculates the vector U which provides minimum of
%PQSQ_NORM(x-U*V')
%   x is n-by-m matrix of data to approximate
%   V is m-by-1 vector to project data
%   potentialFunction is structure which define PQSQ potential function
%       (see definePotentialFunction.m for details)
%   eps is tolerance level for iterations convergence

    [n, m] = size(x);
    XU = zeros(n,1);
    SU = zeros(n,1);
    UU = zeros(n,m);
    
    for i=1:m %Loop of several initiations
        %Initiations
        uk = x(:,i)/V(i);
        %Optimizations
        count = 0;
        while (count<100)
            u1 = uk;
            XU(:)=0;
            SU(:)=0;
            for j=1:m
                inds = identifyIntervals(x(:,j)-uk.*V(j),potentialFunction.intervals(j,:));
                XU = XU + potentialFunction.A(j,inds)'.*(x(:,j))*V(j);
                SU = SU + potentialFunction.A(j,inds)'.*(V(j)*V(j));
            end;
            ind = SU==0;
            uk(ind) = 0;
            uk(~ind) = XU(~ind)./SU(~ind);
            if norm(uk-u1)<eps
                break;
            end
            count = count+1;
        end
        UU(:,i) = uk;
    end

    %Calculate potential for each projection
    err = zeros(n,m);
    for i=1:m
        err(:,i) = PQSQ_Norm( x-UU(:,i)*V', potentialFunction, 2);
    end
    %Search minimal positions
    [~, ind] = min(err,[],2);
    %Get projections
    U = UU(sub2ind([n,m],(1:n)',ind));
end

