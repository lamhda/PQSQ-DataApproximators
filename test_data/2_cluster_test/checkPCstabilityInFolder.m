function [meansL2,meansL1,mediansL2,mediansL1,varsL2,varsL1,madsL2,madsL1,maxDiff,ns] = checkPCstabilityInFolder(folder)

files = dir(folder);

k=1;

for i=1:size(files)
    file = strcat(folder,files(i).name);
    display(sprintf('%i:%s',i,file));
    if(strendswith(file,'.txt'))
    dat = importdata(file);
    ns(k)=size(dat.data,1);
    [corrL2,corrL1,corrInter] = testStabilityOfFirstPCACompare(dat.data,1);
    meansL2(k) = mean(corrL2);
    meansL1(k) = mean(corrL1);
    mediansL2(k) = median(corrL2);
    mediansL1(k) = median(corrL1);
    varsL2(k) = var(corrL2);
    varsL1(k) = var(corrL1);
    madsL2(k) = mad(corrL2);
    madsL1(k) = mad(corrL1);
    maxDiff(k) = max(corrInter);
    k=k+1;
    end
end

    
    
function b = strendswith(s, pat)
%STRENDSWITH Determines whether a string ends with a specified pattern
%
%   b = strstartswith(s, pat);
%       returns whether the string s ends with a sub-string pat.
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 9, 2008
%

%% main

sl = length(s);
pl = length(pat);

b = (sl >= pl && strcmp(s(sl-pl+1:sl), pat)) || isempty(pat);


    