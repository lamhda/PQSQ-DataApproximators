%Calculate statistics
corr = sum(corrs);
corrPQ = sum(corrsPQ);
DI = mean(DunnInd);
DIPQ = mean(DunnIndPQ);
%Restore used number of noise points
nNoise = length(DI) * 20;
noises = 20:20:nNoise;

%Draw graph number of corrected cases
figure;
plot(noises,corr,'b-');
hold on;
plot(noises,corrPQ,'r--');
xlabel('Number of noise points');
set(gca,'FontSize',15);
ylabel('Percent of correct clustering');
legend('kmeans','kmeansPQSQ');
hold off

%Draw graph of Dunn Index
figure;
plot(noises,DI,'b-');
hold on;
plot(noises,DIPQ,'r--');
set(gca,'FontSize',15);
xlabel('Number of noise points');
ylabel('Dunn index');
legend('kmeans','kmeansPQSQ');
hold off;