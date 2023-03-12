function final_prediction = collaborativeFilteringFun(N, P, R, DL, Nf, ...
    lambda)

% This function is used to perform the collaborative filtering based
% collapse fragility assessment

% Developed by GUAN, XINGQUAN @ UCLA in Oct. 2021

% Reference:
% Guan, X. and Burton, H. Collaborative filtering-based collapse fragility
% assessment. (Under review)

% Function arguments:
% N: an integer used for the performance rating
% P: a NIL*NGM matrix that saves the maximum peak story drift ratio where 
% NGM is the number of ground motions and NIL is the number of total intensity levels
% R: a **nIL matrix serves as indicator
% DL: a float number denotes the story drift limit for the collapse
% Nf: number of efatures
% lambda: regularization parameters

% Convert the matrix P into a rating table Y
[NIL, NGM] = size(P);
Y = zeros([NIL, NGM]);
for i = 1:NIL
    for j = 1:NGM
        if R(i,j) == 0
            Y(i,j) = 0;
        else
            temp = P(i,j) / (DL/(N-1));
            Y(i,j) = N+1 - min(floor(temp)+1, N);
        end
    end
end

% Perform the normalization to the matrices
[Ynorm, Ymean] = normalizeRatings(Y, R);

% Set Initial Parameters (Theta, X)
X = randn(NIL, Nf);
Theta = randn(NGM, Nf);
initial_parameters = [X(:); Theta(:)];

% Set options for fmincg
options = optimset('GradObj','on','MaxIter',100);
theta = fmincg(@(t)(cofiCostFunc(t, Ynorm, R, NGM, NIL, Nf, lambda)), ...
    initial_parameters, options);

% Unfold the returned theta back into U and W
X = reshape(theta(1:NIL*Nf), NIL, Nf);
Theta = reshape(theta(NIL*Nf+1:end), NGM, Nf);

% Generate the performance prediction
p = X * Theta';
performance_prediction = round(p + Ymean);

% Assemble the NRHA-based and collaborative filtering based performance
final_prediction = zeros(NIL, NGM);
for i = 1:NIL
    for j = 1:NGM
        if R(i,j) == 0
            % Collaborative filtering generated performance rating
            final_prediction(i,j) = performance_prediction(i, j);
        else
            % NRHA-based performance rating
            final_prediction(i,j) = Y(i,j);
        end
    end      
end

end