%% NORMAL VERSION - 2 PATCHES
%
rng(1);
clear;
close all;
clc;
clifford_signature(4,0);
addpath C:\Users\Harismats\Desktop\Updated_New_ICP\ICP\Paper_TAROS\Functions\

%% Options
show = 1; % for live plotting
method = 1; % 1 - SVD, 2 - CM,

%% Read Point Clouds
ptCloud = pcread('C:\Users\Harismats\Desktop\Point_Clouds\plant.ply');

% Downsample
ptCloud = pcdownsample(ptCloud,"random",0.01);
pcshow(ptCloud)
title('Model Point Cloud')

%% Extract the coordinates
% Initial Point Cloud
xyz_coordinates = ptCloud.Location; % This is Nx3

%% Rescaling
%
% Move source_Data to origin
mean_ptCloud = mean(ptCloud.Location,1)'; % 3x1 vector
ptCloud_origin = ptCloud.Location' - mean_ptCloud; % 3xN - 3x1 vector

% Calculate the distances of every point of the POINT CLOUD from the origin
Distance = sqrt(sum((ptCloud_origin - [0,0,0]').^2, 2));

% Find the max Distance and divide every coordinate by this max Distance
max_Distance = max(Distance);
source_Data_rescaled = (ptCloud_origin ./ max_Distance); % 3xN

%% Show Rescaled Point Cloud
ptCloud_rescaled = pointCloud(source_Data_rescaled');
ptCloud_rescaled.Color = ptCloud.Color; % add the colour information

figure
pcshow(ptCloud_rescaled)
title('Rescaled Model Point Cloud')

%% SURFACE PATCH
%
% Specify a query point and the number of nearest neighbors to be identified.
point = ptCloud_rescaled.Location(1,1:3); % 1st point
K = 100; % Number of neighbouring points

% Get the indices and the distances of K nearest neighboring points.
[index,dists] = findNearestNeighbors(ptCloud_rescaled,point,K);

point2 = ptCloud_rescaled.Location(1690,1:3); % point 690
K = 1150; % Number of neighbouring points

% Get the indices and the distances of K nearest neighboring points.
[index2,dists2] = findNearestNeighbors(ptCloud_rescaled,point2,K);

% Plot
figure
pcshow(ptCloud_rescaled)
hold on
plot3(point(1),point(2),point(3),'*r')
plot3(vertcat(ptCloud_rescaled.Location(index,1),ptCloud_rescaled.Location(index2,1)),...
    vertcat(ptCloud_rescaled.Location(index,2),ptCloud_rescaled.Location(index2,2)),...
    vertcat(ptCloud_rescaled.Location(index,3),ptCloud_rescaled.Location(index2,3)),'*')
legend('Rescaled Model Point Cloud','Query Point','Nearest Neighbours','Location','southoutside','Color',[1 1 1])
hold off
title('Rescaled Model Point Cloud with Nearest Neigbours')

ptCloud_rescaled_segm = [vertcat(ptCloud_rescaled.Location(index,1),ptCloud_rescaled.Location(index2,1)),...
    vertcat(ptCloud_rescaled.Location(index,2),ptCloud_rescaled.Location(index2,2)),...
    vertcat(ptCloud_rescaled.Location(index,3),ptCloud_rescaled.Location(index2,3))];

%% Show Segmented Point Cloud
indices = vertcat(index,index2);
ptCloud_rescaled_segm = pointCloud(ptCloud_rescaled_segm);
ptCloud_rescaled_segm.Color = ptCloud_rescaled.Color(indices,:); % add the colour information

figure
pcshow(ptCloud_rescaled_segm)
title('Rescaled Segmented Model Point Cloud')

%% Transform the Data - Rotation and Translation
% Rotation
x = 90;
y = 60;
z = 120;

% Translation vector
t=[65, 75, 45] ./max_Distance; % Need to divide to be correct

% Transform source Data and calculate initial Transform matrix
[moved_ptCloud_rescaled_segm,T0, Rotation_GT, Translation_GT]= ...
    transformation(ptCloud_rescaled_segm.Location', x,y,z,t); % 3xN

%% Show Segmented Measured Point Cloud
moved_ptCloud_rescaled_segm = moved_ptCloud_rescaled_segm'; % Transpose the data for the in-built function pointCloud

moved_ptCloud_rescaled_segm = pointCloud(moved_ptCloud_rescaled_segm); % Always Nx3
moved_ptCloud_rescaled_segm.Color = ptCloud_rescaled.Color(indices,:); % Same color information

figure
pcshow(moved_ptCloud_rescaled_segm)
title('Segmented Transformed Point Cloud')

%% Show both Point Clouds
figure
pcshowpair(ptCloud_rescaled, moved_ptCloud_rescaled_segm)
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('The Point Clouds for Registration'));
legend('Rescaled Model Point Cloud', 'Rescaled Segmented Measured Point Cloud','TextColor','w')

%% Initializations
% Assign variables - Both must be 3xN
source_Data_new = ptCloud_rescaled.Location';
moved_Data_new = moved_ptCloud_rescaled_segm.Location';

% Transformation Matrix
T_final = eye(4,4);

% Split the Transformation Matrix to Rotation and Translation components
Rotation = T_final(1:3,1:3);
Translation = T_final(1:3,4);

T_t = zeros(4,4);

% Iteration
no_iteration = 0;

%% Main Program
tic;
while 1
    no_iteration = no_iteration + 1; 
    disp(['Ieration =',num2str(no_iteration)]);

    % Convert Colour form unit8 to single/double
    ptCloud_Color = im2double(ptCloud.Color)';
    ptCloud3_Color = im2double(moved_ptCloud_rescaled_segm.Color)';

    k = size(moved_Data_new,2);
    for i = 1:k
        % x coordinate, Channel R
        data_difference_x(1,:) = source_Data_new(1,:) - moved_Data_new(1,i);
        data_difference_R(1,:) = ptCloud_Color(1,:) - ptCloud3_Color(1,i);

        % y coordinate, Channel G
        data_difference_y(2,:) = source_Data_new(2,:) - moved_Data_new(2,i);
        data_difference_G(1,:) = ptCloud_Color(2,:) - ptCloud3_Color(2,i);

        % z coordinate, Channel B
        data_difference_z(3,:) = source_Data_new(3,:) - moved_Data_new(3,i);
        data_difference_B(1,:) = ptCloud_Color(3,:) - ptCloud3_Color(3,i);

        % Distances
        distance_xyz = (data_difference_x(1,:).^2 + data_difference_y(2,:).^2 +...
            data_difference_z(3,:).^2);
        distance_RGB = (data_difference_R(1,:).^2 + data_difference_G(1,:).^2 +...
            data_difference_B(1,:).^2); % (1,:) to all as they are 1 dimensional

        % Distance
        distance = sqrt(distance_xyz + distance_RGB);
        [minimum_dist, minimum_index] = min(distance);   % Find the point with the smallest distance
        data_corresp(:,i) = source_Data_new(:,minimum_index);
    end

    %% Centroids and Decentralization
    center_moved_Data = mean(moved_Data_new,2); % 3x1
    center_corresp_Data = mean(data_corresp,2); % 3x1

    center_moved_data_dec = moved_Data_new - center_moved_Data * ones(1,size(moved_Data_new,2));
    center_corresp_Data_dec = data_corresp - center_corresp_Data * ones(1,size(data_corresp,2));

    %% SVD - CM

    % Covariance Matrix
    H2 = bsxfun(@minus,data_corresp, center_corresp_Data) * ...
            bsxfun(@minus,moved_Data_new, center_moved_Data)';

    if method == 1
        % SVD
        [U,S,V]=svd(H);
        Rotation = U*V';
        if det(Rotation) < 0
            fprintf("det(R) < R, reflection detected!, correcting for it ...\n");
            V(:,3) = V(:,3) * -1;
            Rotation = U*V';
        end
    else
        %% CM Method
        F_matrix = matrix_F_optimized(data_corresp,moved_Data_new);
        G_matrix = matrix_G_optimized(moved_Data_new);
        [~, ~, Rotation] = CM_FGmatrices(F_matrix,G_matrix);

    end

    % Calculate Translation Vector
    Translation = center_corresp_Data - Rotation * center_moved_Data;

    T_t= [Rotation,Translation]; % 3x3
    T_t= [T_t;0,0,0,1]; % 4x4

    % Update transformation matrix
    T_final = T_t * T_final;
    % disp(['Error =',num2str(err)]);
    disp('Transform matrix T =');
    disp(inv(T_final));

    % Update points
    moved_Data_new = Rotation * moved_Data_new + Translation * ones(1,size(moved_Data_new,2));

    %% Live Plotting
    if show == 1
        h=figure(7);
        scatter3(source_Data_new(1,:),source_Data_new(2,:),source_Data_new(3,:),'b.');
        hold on;
        scatter3(moved_Data_new(1,:),moved_Data_new(2,:),moved_Data_new(3,:),'r.');
        hold off;
        title(sprintf(['Measured Point Set to be aligned with Model Point Set\n' ...
        'Characteristic Multivector - Registration after %d iterations.\n'],no_iteration),'fontweight','bold','fontsize',12);

        legend('Model Point Set','Measured Point Set');
        xlabel('X'); ylabel('Y'); zlabel('Z');
        daspect([1 1 1]);
        pause(0.1);
        drawnow
    end

    %% Specify Number of Iterations
    if no_iteration == 120
        break
    end
end

% Print on the screen
disp('True value');
disp(T0);  % Rotating matrix and translation vector true values