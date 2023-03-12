clc; clear;

% Define all paths
currentDir = 'C:\Users\61946\Desktop\CollaborativeFiltering';
rawDataDir = strcat(currentDir, '\RawDataset');
functionDir = strcat(currentDir, '\Functions');
targetDir = strcat(currentDir, '\AllResults\EffectofFittingApproach')

% Load the IDA results
cd (rawDataDir)
maximum_story_drift = importdata('SMRF-3.csv');
[NGM, NIL] = size(maximum_story_drift.data);


Nf = 20;
lambda = 0.1;
N = 5;
DL = 0.1;

% 1 - MLE;
% 2 - Probit
% 3 - SSE
% fitting = 3;


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

% OLD TList
% TList = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50];
% NEW TList
TList = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70];
% NfList = [10, 20, 30, 40, 50];
% lambdaList = [0.01, 0.05, 0.1, 0.15, 0.2];
% NList = [2,3,4,5,6,7,8];


for fitting = [1,2,3]

% Randomness is introduced, so multiple simulations are run here to see the
% overall performance
numberSimulations = 30;

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
    if fitting == 1
        [theta, beta, collapseFractions] = fitLognormalDistribution(completePerformanceRating, intensityLevels);
    end
    if fitting == 2
        [theta, beta, collapseFractions] = fitLognormalDistribution_probit(completePerformanceRating, intensityLevels);
    end
    if fitting == 3
        [theta, beta, collapseFractions] = fitLognormalDistribution_sse(completePerformanceRating, intensityLevels);
    end
    thetaList(Tindex, eachSimulation) = theta;
    betaList(Tindex, eachSimulation) = beta;

    end
end

cd(functionDir)
[NIL, NGM] = size(maximum_story_drift.data')
num_collapse = zeros([1, NIL])
story_drift = maximum_story_drift.data';
for i = 1:NIL
    num_collapse(1,i) = sum(story_drift(i,:) >= 0.1);
end
collapseFractions = num_collapse / NGM;
[theta, beta] = fn_mle_pc(intensityLevels, NGM, num_collapse);

% figure(1)
% font_size = 20;
% scatter(intensityLevels, collapseFractions, 'MarkerEdgeColor', 'red')'
% hold on;
% IM_levels = 0.01:0.01:6.0;
% probability = logncdf(IM_levels, log(theta), beta);
% plot(IM_levels, probability, 'linewidth', 2, 'color', 'blue');
% hold off;
% xlabel('Sa(T1) (g)', 'fontname', 'times new roman', 'fontsize', font_size);
% ylabel('Lognormal fragility', 'fontname', 'times new roman', 'fontsize', font_size);
% % legend('Empirical collapse probability', 'Lognormal fragility');
% set(gca, 'fontname', 'times new roman', 'fontsize', font_size);
% box on;

cd (currentDir)


% figure(2)
% font_size = 20;
% boxplot((thetaList(1:7,:))', 'whisker', inf);
% xlabel('Threshold (T)', 'fontname', 'times new roman', 'fontsize', font_size);
% ylabel('Median collapse capacity (g)', 'fontname', 'times new roman', 'fontsize', font_size);
% % set(gca,'XTickLabel',{'0','5', '10', '15', '20', '25', '30', '35', '40', '45', '50'});
% set(gca,'XTickLabel',{'0','5', '10', '15', '20', '25', '30'});
% set(gca, 'fontname', 'times new roman', 'fontsize', font_size);


% figure(3)
% font_size = 20;
% boxplot((betaList(1:7,:))', 'whisker', inf);
% xlabel('Threshold (T)', 'fontname', 'times new roman', 'fontsize', font_size);
% ylabel('Standard deviation', 'fontname', 'times new roman', 'fontsize', font_size);
% % set(gca,'XTickLabel',{'0','5', '10', '15', '20', '25', '30', '35', '40', '45', '50'});
% set(gca,'XTickLabel',{'0','5', '10', '15', '20', '25', '30'});
% set(gca, 'fontname', 'times new roman', 'fontsize', font_size);

cd(targetDir)
finalTheta = [TList', thetaList];
finalBeta = [TList', betaList];
writematrix(finalTheta, strcat('fitting', num2str(fitting), 'thetaList.csv'));
writematrix(finalBeta, strcat('fitting', num2str(fitting), 'betaList.csv'));

cd (currentDir)
end
