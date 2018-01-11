function [V,U,C] = pcaPQSQ_fast(x, ncomp, varargin)
%pcaPQSQ calculates the first ncomp principal components for dataset in
%matrix x.
%   Input attributes
%       x is n-by-m matrix with n observations and m attributes.
%       ncomp i a number of principal components to calculate
%       Other input features are optional and can be specified by Name,
%       Value pairs or by Value only:
%       'verbose', 1 serves to switch on detailed output form all functions.
%               Mostly used for debuging.
%       'intervals', intervals serves to specify user defined intervals.
%               intervals is matrix each row of which corresponds to one
%               dimension. In each row the first element must be zero.
%   Output arguments
%       V is m-by-ncomp matrix of ncomp columns each of which is the
%           principal component.
%       U is n by ncomp matrix of ncomp columns each of which contains
%           lenths of projection of data point (row) onto corresponding
%           principal component.
%       C is 1-by-m the vector of PQSQ mean (coordinates of central point).
%       D is 1-by-ncomp vctor of values of fraction of explained energy.

    verbose=0;
    intervals = 0;
    number_of_intervals = 5;
    alpha =1;
    useJavaImplementation = 0;
    typeOfProjection = 'PQSQ'; 

    for i=1:length(varargin)
        if strcmp(varargin{i},'verbose')
            verbose = varargin{i+1};
        end
        if strcmp(varargin{i},'intervals')
            intervals = varargin{i+1};
        end
        if strcmp(varargin{i},'javacode')
            useJavaImplementation = 1;
        end
        if strcmp(varargin{i},'potential')
            potential_function_handle = varargin{i+1};
        end
        if strcmp(varargin{i},'number_of_intervals')
            number_of_intervals = varargin{i+1};
        end
        if strcmp(varargin{i},'alpha')
            alpha = varargin{i+1};
        end
        if strcmp(varargin{i},'subtract')
            typeOfProjection = varargin{i+1};
        end
    end

    if ~useJavaImplementation

        if isscalar(intervals)
            potentialFunction = definePotentialFunction(x, number_of_intervals, alpha, potential_function_handle);
        else
            potentialFunction.intervals = [intervals, Inf(size(x,2),1)];
            [potentialFunction.A,potentialFunction.B] = computeABcoefficients_fast(intervals, potential_function_handle);
        end
    

        %calculate central point
        %C = PQSQ_Mean(x,intervals,potential_function_handle);
        C = PQSQ_Mean_fast(x, potentialFunction, varargin);
        %C = mean(x);

        % initiate xwork like x-C
        xwork = bsxfun(@minus,x,C);
        V = zeros(size(x,2),ncomp);
        U = zeros(size(x,1),ncomp);
        
        for i=1:ncomp
            if verbose
                display(sprintf('Component %i',i));
            end

            %Calculate one component
            [Vi, Ui] = firstPrincipalComponentPQSQ_fast(xwork, potentialFunction, varargin{:});

            V(:,i) = Vi;
            U(:,i) = Ui;

            %%%%
            % You might want to use PQSQ-based projections: this should increase precision?
            %%%%%
            if strcmp(typeOfProjection,'L2')
                xwork = xwork - (xwork*Vi)*Vi';
            end
            if strcmp(typeOfProjection,'PQSQ')
                xwork = xwork - Ui*Vi';
            end
        end

    else % java implementation
        if isscalar(intervals)
            potentialFunction = definePotentialFunction(x, 5, @L2);
        end
        intervals = potentialFunction.intervals;
        javaclasspath({'VDAOEngine.jar'});
        pca = vdaoengine.analysis.PCAMethodPQSQ;
        pca.setData(x);
        if strcmp(potential_function_handle,'L1')
            pca.PQSQpotential = vdaoengine.analysis.PQSQPotential.getTrimmedLinearPQSQPotential(x);
        end
        if strcmp(potential_function_handle,'L2')
            pca.PQSQpotential = vdaoengine.analysis.PQSQPotential.getTrimmedQuadraticPQSQPotential(x);
        end
        if strcmp(potential_function_handle,'LSQRT')
            pca.PQSQpotential = vdaoengine.analysis.PQSQPotential.getTrimmedSqrtPQSQPotential(x);
        end

        pca.PQSQpotential.numberOfIntervals = size(intervals,2);
        pca.PQSQpotential.intervals = intervals;
        pca.PQSQpotential.computeABCoefficients();

        pca.verboseMode = 0;
        pca.calcBasis(ncomp);
        C = pca.getBasis().a0;
        C = C';
        U = pca.pointProjections;
        U = U';
        U = U(:,1:ncomp);
        V = pca.getBasis().basis;
    end

end

