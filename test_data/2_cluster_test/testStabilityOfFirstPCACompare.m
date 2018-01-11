function [corrsL2,corrsL1,corrsInter] = testStabilityOfFirstPCACompare(x,num_components)

% by removing one point

[pcref] = pca(x,'NumComponents',num_components);
[pcrefL1] = pcaPQSQ_fast(x,num_components,'potential',@L1,'optimize',0);

n = size(x,1);

for i=1:n
    %display(sprintf('Point %i',i));
    xprime = x([1:i-1,(i+1):end],:);
    [pc] = pca(xprime,'NumComponents',num_components);
    [pcL1] = pcaPQSQ_fast(xprime,num_components,'potential',@L1,'optimize',0);
    for k=1:num_components
        crL2 = corr(pcref(:,k),pc(:,k));
        if(crL2<0) crL2=-crL2; end
        corrsL2(i,k) = 1-crL2;
    end
    for k=1:num_components
        crL1 = corr(pcrefL1(:,k),pcL1(:,k));
        if(crL1<0) crL1=-crL1; end
        corrsL1(i,k) = 1-crL1;
    end
    for k=1:num_components
        crInter = corr(pc(:,k),pcL1(:,k));
        if(crInter<0) crInter=-crInter; end
        corrsInter(i,k) = 1-crInter;
    end
    
end
