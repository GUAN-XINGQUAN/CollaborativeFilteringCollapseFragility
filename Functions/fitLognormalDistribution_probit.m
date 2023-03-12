function [theta, beta, collapseFraction] = ...
    fitLognormalDistribution_probit(completePerformanceRating, intensityLevels)

% This function is used to count the number of collapse cases at each
% intensity.

[NIL, NGM] = size(completePerformanceRating);
collapseCase = zeros([1, NIL]);

for i = 1:NIL
    collapseCase(1,i) = sum(completePerformanceRating(i,:) == 1);
end

collapseFraction = collapseCase / NGM;

[theta, beta] = fn_mle_pc_probit(intensityLevels, NGM, collapseCase);

end

