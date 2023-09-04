function [R, t, RMSE, haus_dist, mean_ae, RMSLE] = ...
    rigidtransformation(Point_Cloud_1, Point_Cloud_2, dim, method)
% Method 1 = SVD, 2 = SVD - Umeyama, 3 = CM F/G Matrices, 4 = MATLAB's
% in-built fucntion

% Make sure point clouds are in single precision as MATLAB uses single
% precision by default when we extract the coordinates. 
% So I guess it is better to use single precision
Point_Cloud_1 = single(Point_Cloud_1);
Point_Cloud_2 = single(Point_Cloud_2);

s = size(Point_Cloud_1);

% Check dimension, transpose inputs if needed
if nargin < 3
    dim = 1;
else
    if ~isscalar(dim)
        error('rigidtform:NonscalarDimension',...
            'Dimension must be a scalar.');
    end
    if dim ~= 1 && dim ~= 2
        error('rigidtform:InvalidDimension','Dimension must be 1 or 2.');
    end
    if dim == 2
        Point_Cloud_1 = Point_Cloud_1';
        Point_Cloud_2 = Point_Cloud_2';
        s = s([2 1]);
    end
end

if s(2) > s(1) && s(2) > 2^15
    warning('rigidtform:LargeDimension',...
        ['Number of dimensions (%d) is very large. Datasets may be '...
        'transposed.'], s(2));
end

%% Ensure Data are Nx3
s2 = size(Point_Cloud_2);

if s(2) ~=3
    Point_Cloud_1 = Point_Cloud_1';
end
if s2(2) ~=3
    Point_Cloud_2 = Point_Cloud_2';
end

%% Calculate centroids
c1 = mean(Point_Cloud_1,1);
c2 = mean(Point_Cloud_2,1);

%% Center data and form the Cross-Covariance Matrix
% 1st Way - FASTER
H = bsxfun(@minus,Point_Cloud_1,c1)'*bsxfun(@minus,Point_Cloud_2,c2);

if method == 1
    % Calculate optimal rotation matrix
    %% SVD
    [U,~,V] = svd(H);
    R = U*V';

    % Correct R if needed to ensure right-handed coordinate system
    if s(2) == 3 && det(R) < 0
        V(:,end) = -V(:,end);
        R = U*V';
    end

elseif method == 2
    % Umeyama
    n = size(Point_Cloud_1,1);
    m = size(Point_Cloud_1,2);
    sigma = 1/n*H;
    [U,~,V] = svd(sigma);

    % Define S
    S = eye(m);
    if det(sigma) < 0 || (rank(sigma) == m-1 && det(U)*det(V) < 0)
        S(m,m) = -1;
    end

    R = (U*S*V');

elseif method ==3
    %% F/G Matrices
    F_matrix = matrix_F_optimized(Point_Cloud_1, Point_Cloud_2);
    G_matrix = matrix_G_optimized(Point_Cloud_2);

    % We can REPLACE F_matrix with H for faster computations - THEY ARE THE SAME
    [~, ~, R] = CM_FGmatrices(F_matrix,G_matrix);
    % Convert Rotation to single for fair result comparison as everything
    % is calculated using single precision and Clifford needs double
    R = single(R);

else
    [~,transformed_data_MATLAB, transformation_MATLAB] = ...
        procrustes(Point_Cloud_2,Point_Cloud_1);
    R = transformation_MATLAB.T;
end

% Calculate optimal translation vector
t = c2-c1*R;

% Registration Result
registration_result = Point_Cloud_1 * R + t;

%% Calculate root mean squared distance error for each point, if requested
if method == 1 || method ==2 || method ==3
    % Calculate squared error between actual and prediction
    score = (Point_Cloud_2(:) - (registration_result(:))).^2;
    % Calculate Mean Squared Error
    MSE = mean(score);
    % Calculate RMSE
    RMSE = sqrt(MSE);
else
    % MATLAB
    score = (Point_Cloud_2(:) - (transformed_data_MATLAB(:))).^2;
    % Calculate Mean Squared Error
    MSE = mean(score);
    % Calculate RMSE
    RMSE = sqrt(MSE);
end

%% Hausdorff Distance
if method == 1 || method ==2 || method ==3
    haus_dist = hausdorff_distance2(Point_Cloud_2,registration_result);
else
    haus_dist = hausdorff_distance2(Point_Cloud_2,transformed_data_MATLAB);
end

%% Mean Absolute Error (MAE)
if method == 1 || method ==2 || method ==3
    mean_ae = mean_absolute_error(Point_Cloud_2,registration_result);
else
    mean_ae = mean_absolute_error(Point_Cloud_2,transformed_data_MATLAB);
end

%% Root Mean Squared Log Error
if method == 1 || method ==2 || method ==3
    % First Calculate the squared log error between actual and prediction
    score2 = (log(1+Point_Cloud_2(:)) - log(1+registration_result(:))).^2;
    % Calculate Mean squared log error
    MSLE = mean(score2);
    % Calculate RMSLE
    RMSLE = sqrt(MSLE);
else
    score2 = (log(1+Point_Cloud_2(:)) - log(1+transformed_data_MATLAB(:))).^2;
    % Calculate Mean squared log error
    MSLE = mean(score2);
    % Calculate RMSLE
    RMSLE = sqrt(MSLE);
end

R = R';
end
