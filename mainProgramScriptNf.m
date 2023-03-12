clc; clear;

%% Read the original dataset
% Define all paths
currentDir = 'C:\Users\61946\Desktop\CollaborativeFiltering';
rawDataDir = strcat(currentDir, '\RawDataset');
functionDir = strcat(currentDir, '\Functions');
targetDir = strcat(currentDir, '\AllResults\EffectofNf')

% Load the IDA results
cd (rawDataDir)
maximum_story_drift = importdata('SMRF-3.csv');
[NGM, NIL] = size(maximum_story_drift.data);


%% Define all input parameters
% Nf = 20;
% lambda = 0.1;
N = 2;
DL = 0.1;

% ----------------- These are for OLD buildings --------------------
% SMRF-1
% SaMCE = 1.22;
% SMRF-2
% SaMCE = 0.81;
% SMRF-3
% SaMCE = 0.51;
% -----------------------------------------------------------------


% ----------------- These are for NEW buildings -------------------
% SMRF-1
% SaMCE = 1.623;
% SMRF-2
% SaMCE = 1.078;
% SMRF-3
SaMCE = 0.673;
% ----------------------------------------------------------------

IDA_scales = [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300];

intensityLevels = SaMCE * IDA_scales / 100;

% TList = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70];
% New TList
TList = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90];
NfList = [10, 20, 30, 40, 50];
lambdaList = [0.01, 0.05, 0.1, 0.15, 0.2];
% NList = [2,3,4,5,6,7,8,9,10];

%% Collaborative filtering based collapse fragility assessment
% Randomness is introduced, so multiple simulations are run here to see the
% overall performance
numberSimulations = 100;

for Nf = NfList
for lambda = lambdaList
thetaList = zeros([length(TList), numberSimulations]);
betaList = zeros([length(TList), numberSimulations]);

for Tindex = 1:length(TList)
    
    T = TList(Tindex)

    for eachSimulation = 1:numberSimulations

    % 1st module: prepare the matrix P and R
    P = maximum_story_drift.data';
    R = zeros([NIL, NGM]);
    for i = 1:NIL
        for j = 1:NGM
            n = randi(100);
            if n > T
                R(i, j) = 1;
            else
                P(i, j) = -1;
            end
        end
    end
    
    % 2nd module: formulate the collaborative filtering algorithm
    cd(functionDir)
    completePerformanceRating = collaborativeFilteringFun(N, P, R, DL, Nf, lambda);
    
    % 3rd module: fit the empirical collapse probability using lognormal distribution
    % function
    [theta, beta, collapseFractions] = fitLognormalDistribution(completePerformanceRating, intensityLevels);
    
    thetaList(Tindex, eachSimulation) = theta;
    betaList(Tindex, eachSimulation) = beta;

    end
end

%% Save results
cd(targetDir)
finalTheta = [TList', thetaList];
finalBeta = [TList', betaList];
writematrix(finalTheta, strcat('Nf', num2str(Nf), 'lambda', num2str(lambda),'thetaList.csv'));
writematrix(finalBeta, strcat('Nf', num2str(Nf), 'lambda', num2str(lambda), 'betaList.csv'));
end
end
cd(currentDir)