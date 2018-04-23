addpath ../../Release_v1.0/

% This script 
% 1. Generate/load data set for k-means test,
% 2. Test standard Matlab kmeans and kmeansPQSQ
% 3. Draw figures for each test if required
% 4. Calculate stratistics for applied tests
% 5. Draw graph of fraction of correct clusterings and Dunn index

% Define parameters of script
generate = 1;   % 0 to load file data.mat and non zero to genrate randomly
nClust = 100;   % Number of points in each cluster
nNoise = 200;   % Number of noise points in set with maximal noise.
minNoise = 20;  % Minimal number of noise points
stepNoise = 20; % Step to increase number of noise points in test sets

% Prefix of file name. If prefix is not empty then the best clustering of
% each test will be written into two files with name 
% [fileNamePrefix, '-ML-n-', mat2str(n), '-t-', mat2str(t)]
% [fileNamePrefix, '-PQ-n-', mat2str(n), '-t-', mat2str(t)]
% and with two etensions eps and png. n is number of noise points and t is
% number of test.
fileNamePrefix = [];

% Number of tests is number of attempts for each test set
nTests = 100;

% Number of repetitions of procedure to select the best in each kmeans call
nReps = 5;

% Tolerance for correct clustering identification: Clustering is correct if
% at least tol * nClust of r1 points are located in one cluster and at least
% tol * nClust of r2 points are located in another cluster.
tol = 0.95;

% Type of PQSQ potential
pot = @L1;

%% Part 1. Load or generate data
if generate == 0
    %Load
    load('data.mat');
    nNoise = size(noise, 1);
else
    % Data generation
    % Check data
    nNoise = floor(nNoise);
    if nNoise < 0
        error('Number of noise points must be at least non negative');
    end
    
    % Means of distributions for two clusters
    mu1 = zeros(1, 2);
    mu2 = zeros(1, 2);
    mu1(1) = -1; 
    mu2(1) = 1;

    % List of standard deviasions for each dimension for clusters 
    s1 = repmat(0.1, 1, 2);
    s2 = repmat(0.1, 1, 2); 

    % Generate normally distributed points for each cluster
    r1 = mvnrnd(mu1, s1, nClust); 
    r2 = mvnrnd(mu2, s2, nClust);
    
    % Generate noise points with Laplassian distribution with zero mean and
    % unit standard deviation
    noise = laprnd(nNoise, 2);
    %Change standard deviation by 2 in x direction and by 4 in y direction
    noise(:, 1) = noise(:, 1) * 2;
    noise(:, 2) = noise(:, 2) * 4;
    
    % Save generated data to file
    save('data.mat', 'noise', 'r1', 'r2');
    clear mu1 mu2 s1 s2     
end    

%Chack correctness of minNoise and stepNoise
if minNoise < 0 || minNoise > nNoise
    minNoise = floor(nNoise / 10);
end
if stepNoise < 0 || stepNoise > nNoise
    stepNoise = minNoise;
end
 
%% Part 2. Test
% Calculate sizes of expected ouputs
% Number of test sets is 
nSets = floor((nNoise - minNoise) / stepNoise) + 1;

% Create arrays for test results
DunnInd = zeros(nTests, nSets);     % Dunn index for standard kmeans
DunnIndPQ = zeros(nTests, nSets);   % Dunn index for kmeansPQSQ
correct = zeros(nTests, nSets);     % Number of correct clusterings for standard kmeans
correctPQ = zeros(nTests, nSets);   % Number of correct clusterings for kmeansPQSQ

% Loops of testing
for kSet = 1:nSets
    display(sprintf('Iteration %i/%i',kSet,nSets));
    nN = minNoise + (kSet - 1) * stepNoise;
    data = [r1; r2; noise(1:nN, :)];
    for kTest = 1: nTests
        % Call standard kmeans
        [idx, cent] = kmeans(data, 2, 'Replicates', nReps);
        % Call PQSQ kmeans
        [idxPQ, centPQ] = kmeansPQSQ(data, 2, pot, 'Replicates', nReps);
        if isempty(fileNamePrefix)
            [correct(kTest, kSet), DunnInd(kTest, kSet)] =...
                statCalc( data, idx, cent, nClust, []);
            [correctPQ(kTest, kSet), DunnIndPQ(kTest, kSet)] =...
                statCalc( data, idxPQ, centPQ, nClust, []);
        else
            fName = [fileNamePrefix, '-ML-n-', mat2str(nM), '-t-', mat2str(kTest)];
            [correct(kTest, kSet), DunnInd(kTest, kSet)] =...
                statCalc( data, idx, cent, nClust, fName);
            fName = [fileNamePrefix, '-PQ-n-', mat2str(nM), '-t-', mat2str(kTest)];
            [correctPQ(kTest, kSet), DunnIndPQ(kTest, kSet)] =...
                statCalc( data, idxPQ, centPQ, nClust, fName);
        end
    end
end

%% Part 3. Draw graphs
%Calculate statistics
corr = sum(correct);
corrPQ = sum(correctPQ);
DI = mean(DunnInd);
DIPQ = mean(DunnIndPQ);
%Restore used number of noise points
nNoise = length(DI) * 20;
noises = minNoise:stepNoise:nNoise;

%Draw graph number of correct cases
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