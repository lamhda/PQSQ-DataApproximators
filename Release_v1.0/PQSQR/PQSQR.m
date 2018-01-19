function y = PQSQR(x, potentialFunction)
%PQSQR calculates the PQSQ energy potentials u(x) in accordance with
%formula (2) in paper for kth coordinate
%   y is vector of PQSQ potentials for all points.
%   x is vector of squared length of row vectors of data or the matrix of
%       data. 
%   potentialFunction is structure which define PQSQ potential function
%       (see PQSQR_definePotentialFunction.m for details)

    if size(x,2)>1
        %x contains data matrix
        x = sqrt(sum(x.^2,2));
    end
    inds = identifyIntervals(x,potentialFunction.sqint);
    y = potentialFunction.A(inds)'.*x + potentialFunction.B(inds)';
end