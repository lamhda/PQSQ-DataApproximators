function [A,B] = computeABcoefficients_fast(intervals, potential_function_handle)
%computeABcoefficients_fast calculates the coefficients a abd b for
%quadratic fragments of potential function.
%   intervals is the m-by-K matrix of intervals' boudaries without final
%       infinit boundary.
%   potential_function_handle is a handle of majorant function.

    %Get dimensions of intervals
    [m, K] = size(intervals);

    %Preallocate memory
    A = zeros(m,K);
    B = zeros(m,K);

    %Calculate value of function all boundaries
    pxk = potential_function_handle(intervals);
    sxk = intervals.^2;

    A(:,1:K-1) = (pxk(:,1:K-1)-pxk(:,2:K))./(sxk(:,1:K-1)-sxk(:,2:K));
    B(:,1:K-1) = (pxk(:,2:K).*sxk(:,1:K-1)-pxk(:,1:K-1).*sxk(:,2:K))./...
        (sxk(:,1:K-1)-sxk(:,2:K));
    B(:,K) = pxk(:,K);
end