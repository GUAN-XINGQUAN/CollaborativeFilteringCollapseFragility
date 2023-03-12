% Load data
load('IDAResults.mat');

% From the matrix, we can compute statistics like average rating.
fprintf('Average rating for scale 1: %f / 5\n\n', mean(Y(1, R(1, :))));

%  We can "visualize" the ratings matrix by plotting it with imagesc
imagesc(Y);
ylabel('Intensity scales');
xlabel('Ground motion');

% Initialize my ratings
% my_ratings = zeros(20, 1);

%  Y is a 20x44 matrix, containing ratings (1-5) of 20 movies by 44 users
%  R is a 20x44 matrix, where R(i,j) = 1 if and only if user j gave a rating to movie i
%  Add our own ratings to the data matrix
% Y = [my_ratings Y];
% R = [(my_ratings ~= 0) R];

%  Normalize Ratings
[Ynorm, Ymean] = normalizeRatings(Y, R);

%  Useful Values
num_users = size(Y, 2);
num_movies = size(Y, 1);
num_features = 10;

% Set Initial Parameters (Theta, X)
X = randn(num_movies, num_features);
Theta = randn(num_users, num_features);
initial_parameters = [X(:); Theta(:)];

% Set options for fmincg
options = optimset('GradObj','on','MaxIter',100);

% Set Regularization
lambda = 0.05;
theta = fmincg(@(t)(cofiCostFunc(t, Ynorm, R, num_users, num_movies, num_features,lambda)), initial_parameters, options);

% Unfold the returned theta back into U and W
X = reshape(theta(1:num_movies*num_features), num_movies, num_features);
Theta = reshape(theta(num_movies*num_features+1:end), num_users, num_features);

% Predict
p = X * Theta';
my_predictions = p + Ymean;

% Round the prediction results
predictions = zeros(size(my_predictions));
for i = 1:num_movies
    for j = 1:num_users
        predictions(i, j) = round(my_predictions(i,j));
    end
end

% Check the prediction results
diff = 0;
for i = 1:num_movies
    for j = 1:num_users
        if R(i,j) == 0
           diff = diff + abs(predictions(i,j) - real_res(i, j)); 
           fprintf('Real = %f, Prediction = %f\n', real_res(i, j), predictions(i, j));
        end
    end
end