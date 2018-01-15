function [idxbest, Cbest, sumDbest, Dbest] = kmeansPQSQ(X, k, potential_function_handle, varargin)
%kmeansPQSQ K-means clustering with PQSQ potential instead of distance.
%   IDX = kmeansPQSQ(X, K, @fun) partitions the points in the N-by-P data
%   matrix X into K clusters. This partition minimizes the sum, over all
%   clusters, of the within-cluster sums of point-to-cluster-centroid PQSQ
%   potential. Rows of X correspond to points, columns correspond to
%   variables.  Note: when X is a vector, kmeansPQSQ treats it as an N-by-1
%   data matrix, regardless of its orientation.  kmeansPQSQ returns an
%   N-by-1 vector IDX containing the cluster indices of each point.
%
%   kmeansPQSQ do not work with missed data. All NaNs are treated as missed
%   data and caused error.
%
%   [IDX, C] = kmeansPQSQ(X, K, @fun) also returns the K cluster centroid
%   locations in the K-by-P matrix C.
%
%   [IDX, C, SUMD] = kmeansPQSQ(X, K, @fun) also returns the within-cluster
%   sums of point-to-centroid distances in the K-by-1 vector sumD.
%
%   [IDX, C, SUMD, D] = kmeansPQSQ(X, K, @fun) also returns distances from
%   each point to every centroid in the N-by-K matrix D.
%
%   [ ... ] = kmeansPQSQ(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control the iterative algorithm
%   used by KMEANS.  Parameters are:
%
%   'intervals', intervals pair serves to specify user defined intervals.
%       intervals is P-by-R matrix each row of which corresponds to one
%       dimension. In each row the first element must be zero. 
%   'number_of_intervals', number_of_intervals pair specifies the number of
%       intervals to automatic interval calculation.  
%   'intshrinkage', delta serves to specify delta which is coefficient for
%       intervals shrinkage (see argument delta in defineIntervals.m). 
%
%   'Start' - Method used to choose initial cluster centroid positions,
%       sometimes known as "seeds". Choices are:
%          'plus'    - The Default. Select K observations from X according
%                      to the k-means++ algorithm: the first cluster center
%                      is chosen uniformly at random from X, after which
%                      each subsequent cluster center is chosen randomly
%                      from the remaining data points with probability
%                      proportional to its distance from the point's
%                      closest existing cluster center.
%          'sample'  - Select K observations from X at random.
%          'uniform' - Select K points uniformly at random from the range
%                      of X.  Not valid for Hamming distance.
%          'cluster' - Perform preliminary clustering phase on random 10%
%                      subsample of X.  This preliminary phase is itself
%                      initialized using 'sample'.
%           matrix   - A K-by-P matrix of starting locations.  In this case,
%                      you can pass in [] for K, and KMEANS infers K from
%                      the first dimension of the matrix.  You can also
%                      supply a 3D array, implying a value for 'Replicates'
%                      from the array's third dimension.
%
%   'Replicates', replicates is Number of times to repeat the clustering,
%       each with a new set of initial centroids.  A positive integer,
%       default is 1. 
%
%   'MaxIter', MaxIter  - Maximum number of iterations allowed.  Default is
%       100. 
%
%   'eps', positive number is tolerance level for PQSQ mean calculation.
%
%   Example:
%
%       X = [randn(20,2)+ones(20,2); randn(20,2)-ones(20,2)];
%       [cidx, ctrs] = kmeansPQSQ(X, 2, @L1, 'Replicates',5);
%
%       plot(X(cidx==1,1),X(cidx==1,2),'r.', ...
%            X(cidx==2,1),X(cidx==2,2),'b.', ctrs(:,1),ctrs(:,2),'kx');
%
%
%   kmeansPQSQ uses a iterative algorithm to minimize the sum of
%   point-to-centroid PQSQ potentials, summed over all K clusters.  Method
%   uses what the literature often describes as "batch" updates, where each
%   iteration consists of reassigning points to their nearest cluster
%   centroid, all at once, followed by recalculation of cluster centroids.
%   If during the optimization one or several of centroids do not have
%   associated points then this centroid is dropped, its coordinates are
%   the same as before lost of the last point.

    if nargin < 3
        error('Too few inputs');
    end

    if any(any(isnan(X)))
        error('Missed data are unacceptable');
    end

    % n points in p dimensional space
    [n, p] = size(X);

    if n==1 || p==1
        X = X(:);
        [n, p] = size(X);
    end

    %Parse vargin
    intervals = 0;
    number_of_intervals = 5;
    delta = 1;
    start = 'plus';
    reps = 1;
    maxIter = 100;
    eps = 0.001;

    %Process Name, Value pairs
    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'intervals')
            intervals = varargin{i+1};
        elseif strcmpi(varargin{i},'number_of_intervals')
            number_of_intervals = varargin{i+1};
        elseif strcmpi(varargin{i},'intshrinkage')
            delta = varargin{i+1};
        elseif strcmpi(varargin{i},'Start')
            start = varargin{i+1};
        elseif strcmpi(varargin{i},'Replicates')
            reps = varargin{i+1};
        elseif strcmpi(varargin{i},'MaxIter')
            maxIter = varargin{i+1};
        elseif strcmpi(varargin{i},'eps')
            eps = varargin{i+1};
        end
    end

    %Form potential function structure
    if isscalar(intervals)
        potentialFunction = definePotentialFunction(X, number_of_intervals,...
            potential_function_handle, delta);
    else
        potentialFunction.intervals = [intervals, Inf(size(X,2),1)];
        [potentialFunction.A,potentialFunction.B] = ...
            computeABcoefficients(intervals, potential_function_handle);
    end

    %Start
    if ischar(start)
        startNames = {'uniform','sample','cluster','plus','kmeans++'};
        j = find(strncmpi(start,startNames,length(start)));
        if length(j) > 1
            error(['Ambiguous start:', start]);
        elseif isempty(j)
            error(['Unknown start', start]);
        elseif isempty(k)
            error('Undefined number of clusters');
        end
        start = startNames{j};
        if strcmp(start, 'uniform')
            Xmins = min(X,[],1);
            Xmaxs = max(X,[],1);
        end
    elseif isnumeric(start)
        CC = start;
        start = 'numeric';
        if isempty(k)
            k = size(CC,1);
        elseif k ~= size(CC,1);
            error('Start has bad number of rows');
        elseif size(CC,2) ~= p
            error('Start has bad number of columns');
        end
        if reps==1
            reps = size(CC,3);
        elseif reps ~= size(CC,3);
            error('Incorrect number of sheets in 3D array start');
        end
    else
        error(message('stats:kmeans:InvalidStart'));
    end

    %Check of types of some input atributes
    if ~(isscalar(maxIter) && isnumeric(maxIter) && isreal(maxIter) && ...
            maxIter > 0 && (round(maxIter)==maxIter))
        error('MaxIter must be positive integer');
    end

    if ~(isscalar(reps) && isnumeric(reps) && isreal(reps) && ...
            reps > 0 && (round(reps)==reps))
        error('Replicates must be positive integer');
    end

    if ~(isscalar(k) && isnumeric(k) && isreal(k) && k > 1 && (round(k)==k))
        error('Number of clusters k must be positive integer');
    elseif n < k
        error(message('Too many clusters'));
    end

    %Prepare to loop
    idxbest = []; 
    Cbest = []; 
    sumDbest = [];
    Dbest = [];
    BestTotal = Inf;
    S = RandStream.getGlobalStream;
    for rep = 1:reps
        %initiate
        switch start
            case 'uniform'  
                C = Xmins(ones(k,1),:) + rand(S,[k,p]).*(Xmaxs(ones(k,1),:)-Xmins(ones(k,1),:));
            case 'sample'
                C = X(randsample(S,n,k),:);
            case 'cluster'
                Xsubset = X(randsample(S,n,floor(.1*n)),:);
                [~, C] = kmeansPQSQ(Xsubset, k, potential_function_handle, varargin{:}, 'start','sample', 'replicates',1);
            case 'numeric'
                C = CC(:,:,rep);
            case {'plus','kmeans++'}
                % Select the first seed by sampling uniformly at random
                index = zeros(k,1);
                [C(1,:), index(1)] = datasample(S,X,1);

                % Select the rest of the seeds by a probabilistic model
                for ii = 2:k                    
                    sampleProbability = min(distfun(X,C(1:ii-1,:),potentialFunction),[],2);
                    denominator = sum(sampleProbability);
                    if denominator==0 || isinf(denominator) || isnan(denominator)
                        C(ii:k,:) = datasample(S,X,k-ii+1,1,'Replace',false);
                        break;
                    end
                    sampleProbability = sampleProbability/denominator;

                    [C(ii,:), index(ii)] = datasample(S,X,1,1,'Replace',false,...
                        'Weights',sampleProbability);
                end
        end
        if isa(X,'single')
            C = single(C);
        end
        if ~isfloat(C)      % X may be logical
            C = double(C);
        end
        
        %Centroids are initialized. Start optimization.
        %create list of associated centroids
        indx = repmat(-1,n,1);
        while (true)
            %Associate points with centroids 
            D = distfun(X, C, potentialFunction);
            [d, idx] = min(D, [], 2);
            %Calculate the number of points associated with each centroid.
            m = accumarray(idx,1,[k,1]);
            %check the convergence
            if ~any(indx-idx)
                break;
            end
            indx = idx;
            %Recalculate centroids
            for kk = 1:k
                %If centroid does not have associated points we drop it
                %without renewing.
                if m(kk)==0
                    continue;
                end
                ind = indx==kk;
                C(kk,:) = PQSQ_Mean(X(ind,:), potentialFunction, eps);
            end
        end
        
        %Select the best result. Now we have D, C, idx. 
        %Total sum is
        dd = sum(d);
        if BestTotal>dd
            idxbest = indx;
            Cbest = C;
            Dbest = D;
            sumDbest = accumarray(indx,d,[k,1]); 
            BestTotal = dd;
        end
    end

end

function D = distfun(X, C, potentialFunction)
%DISTFUN Calculate point to cluster centroid PQSQ potentials.
    k = size(C,1);
    D = zeros(size(X,1),k);

    for i = 1:k
        D(:,i) = PQSQ_Norm(bsxfun(@minus,X,C(i,:)), potentialFunction, 2);
    end
end % function


