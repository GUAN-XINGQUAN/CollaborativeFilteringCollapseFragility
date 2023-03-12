% Load data
% load('IDAResults.mat');

num_users = size(Y, 2);
num_movies = size(Y, 1);

count = 0
predictions = zeros(size(Y));
for i = 1:num_movies
    for j = 1:num_users
        if R(i, j) == 1
            predictions(i, j) = Y(i, j);
        else
            count = count + 1;
            predictions(i, j) = mean(Y(:, j));
        end
    end
end