function [theta, beta, collapseFraction] = ...
    fitLognormalDistribution_sse(completePerformanceRating, intensityLevels)

% This function is used to count the number of collapse cases at each
% intensity.

[NIL, NGM] = size(completePerformanceRating);
collapseCase = zeros([1, NIL]);

for i = 1:NIL
    collapseCase(1,i) = sum(completePerformanceRating(i,:) == 1);
end

collapseFraction = collapseCase / NGM;

[theta, beta] = fn_sse_pc(intensityLevels, NGM, collapseCase);

end

