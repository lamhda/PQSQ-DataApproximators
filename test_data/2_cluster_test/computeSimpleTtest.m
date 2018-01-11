function [ttest] = computeSimpleTtest(p,npoints)

mean1 = mean(p(1:npoints,:));
mean2 = mean(p(npoints+1:2*npoints,:));
st = std(p(1:npoints,:));
dist = norm(mean1-mean2);

ttest = dist/st(1)/sqrt(1/npoints);
