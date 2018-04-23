function [ corCnt, DI ] = statCalc( data, idx, centroids, nClust, fName, lims )
% data is data matrix
% idx is index of clasters for each point
% Centroids is matrix of centroids (each row is one centroid)
% fName is file name to write figures. Empty fle name prevent figure
%   creation.
% nClust is number of elements in original clusters.
% lims is array of limits for figure

    %Calculate Dunn index
    DI = indexDN(data,idx);

    %Threshold
    th = 0.95;
    lth = (th * 1 + (1 - th) * 2) * nClust;
    gth = (th * 2 + (1 - th) * 1) * nClust;
    
    %Is it acceptable?
    a = sum(idx(1:nClust));
    b = sum(idx(nClust + 1:2 * nClust));
    corCnt = (a < lth && b > gth) || (b < lth && a > gth);
    
    if isempty(fName)
        return;
    end
    
    %Draw and save figures
    plot(data(idx == 1, 1),data(idx == 1, 2), 'ro');
    hold on;
    plot(data(idx == 2, 1),data(idx == 2, 2), 'b+');
    
    plot(centroids(1, 1), centroids(1, 2), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
    plot(centroids(2, 1), centroids(2, 2), 'bx', 'MarkerSize', 10, 'LineWidth', 3);
    if nargin > 5
        axis(lims);
    end
    set(gca,'FontSize',15);
    hold off;
    
    %Save png and eps figures
    print(gcf,'-depsc', '-noui', '-loose', fName);
    print(gcf,'-dpng', '-noui', '-loose', fName);
end

