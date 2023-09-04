%% Clear Everything
clear;
close all;
clc;
clifford_signature(4,0);
addpath C:\Users\Harismats\Desktop\Updated_New_ICP\ICP\Functions\
addpath C:\Users\Harismats\Desktop\Updated_New_ICP\ICP\One_One_Correspondences\Paper\Functions_Paper

%% Replication of Results
rng(100)

%% Read Point Clouds
ptCloud = pcread('C:\Users\Harismats\Desktop\Updated_New_ICP\helix.ply');
nsample, dragon 0.1, armadillo 0.1

f = figure;
colormap(f,'copper')
pcshow(ptCloud, 'BackgroundColor',[0.9 0.9 0.9]);
xlabel('X'); ylabel('Y'); zlabel('Z');
set(gca,'FontSize',16); % axes font size

%% Extract the coordinates
% Initial Point Cloud
xyz_coordinates = ptCloud.Location; % This is Nx3

model_Data = xyz_coordinates;
model_Data=model_Data'; % 3xN

size_model_Data = size(model_Data);

%% Number of Iterations
no_test_cases = 200;

%% Preallocation
measured_Data = zeros(3,size_model_Data(2),no_test_cases);
x = zeros(1,1,no_test_cases);
y = zeros(1,1,no_test_cases);
z = zeros(1,1,no_test_cases);
t = zeros(1,3,no_test_cases);
Rotation_GT = zeros(3,3,no_test_cases);
Translation_GT = zeros(1,3,no_test_cases);
P2 = zeros(3,size_model_Data(2),no_test_cases);
R_SVD = zeros(3,3,no_test_cases);
R_CM = zeros(3,3,no_test_cases);
R_MATLAB = zeros(3,3,no_test_cases);
t_SVD = zeros(1,3,no_test_cases);
t_CM = zeros(1,3,no_test_cases);
t_MATLAB = zeros(1,3,no_test_cases);

RMSE_SVD = zeros(1,no_test_cases);
RMSE_CM = zeros(1,no_test_cases);
RMSE_MATLAB = zeros(1,no_test_cases);
haus_dist_SVD = zeros(1,no_test_cases);
haus_dist_CM = zeros(1,no_test_cases);
haus_dist_MATLAB = zeros(1,no_test_cases);
mean_ae_SVD = zeros(1,no_test_cases);
mean_ae_CM = zeros(1,no_test_cases);
mean_ae_MATLAB = zeros(1,no_test_cases);
RMSLE_SVD = zeros(1,no_test_cases);
RMSLE_CM = zeros(1,no_test_cases);
RMSLE_MATLAB = zeros(1,no_test_cases);

convergence_error_SVD = zeros(1,no_test_cases);
convergence_error_CM = zeros(1,no_test_cases);
convergence_error_MATLAB = zeros(1,no_test_cases);
convergence_error_SVD_chordal = zeros(1,no_test_cases);
convergence_error_CM_chordal = zeros(1,no_test_cases);
convergence_error_MATLAB_chordal = zeros(1,no_test_cases);
translation_error_SVD = zeros(1,no_test_cases);
translation_error_CM = zeros(1,no_test_cases);
translation_error_MATLAB = zeros(1,no_test_cases);

time_finish_SVD = zeros(1,no_test_cases);
time_finish_CM = zeros(1,no_test_cases);
time_finish_MATLAB = zeros(1,no_test_cases);

%% Transform the Data - Rotation and Translation
for k = 1:no_test_cases
    fprintf('Number of Test: %.d \n', k)
    % Rotation angle
    x(:,:,k) = randi([0 180]);
    y(:,:,k) = randi([0 180]);
    z(:,:,k) = randi([0 180]);

    % Translation vector
    % Generate values from the uniform distribution on the interval [a, b]:
    % r = a + (b-a).*rand(100,1);
    t(:,:,k) = [0 + (500).*rand(1,1),0 + (500-0).*rand(1,1), 0 + (500-0).*rand(1,1)];

    %% Transform source Data and calculate initial Transform matrix
    [measured_Data(:,:,k), ~, Rotation_GT(:,:,k), Translation_GT(:,:,k)]=...
        transformation(model_Data, x(:,:,k),y(:,:,k),z(:,:,k),t(:,:,k));

    %% Show both Point Clouds

    ptCloud2 = pointCloud(permute(measured_Data(:,:,k),[2 1 3])); % Always Nx3
    figure
    pcshowpair(ptCloud, ptCloud2)
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('Test Case %d'), k);
    legend({'Model Point Cloud', 'Measured Point Cloud'},'TextColor','w')


    %% Each row represents the coordinates of a point
    P1 = model_Data';% Nx3 - Uncomment when using wrong correspondences
    P2(:,:,k) = measured_Data(:,:,k);

    m = size(model_Data,1);  % Number of dimensions (integer > 0)
    n = size(model_Data,2); % Number of points (integer > 0)

    %% Find optimal transform
    %% SVD
    time_start_SVD = tic;
    [R_SVD(:,:,k), t_SVD(:,:,k), RMSE_SVD(:,k), haus_dist_SVD(:,k),...
        mean_ae_SVD(:,k), RMSLE_SVD(:,k)] = ...
        rigidtransformation(P1,P2(:,:,k),1,1); % 3rd argument (dim) = 1, 4th argument (method) = 1 - SVD, 2 - CM, 3 - F/G Matrices
    time_finish_SVD(k) = toc(time_start_SVD);

    %% CM
    time_start_CM = tic;
    [R_CM(:,:,k), t_CM(:,:,k), RMSE_CM(:,k), haus_dist_CM(:,k), ...
        mean_ae_CM(:,k), RMSLE_CM(:,k)] = ...
        rigidtransformation(P1,P2(:,:,k),1,3);
    time_finish_CM(k) = toc(time_start_CM);

    %% MATLAB
    time_start_MATLAB = tic;
    [R_MATLAB(:,:,k), ~, RMSE_MATLAB(:,k),  haus_dist_MATLAB(:,k),...
        mean_ae_MATLAB(:,k), RMSLE_MATLAB(:,k)] =...
        rigidtransformation(P1,P2(:,:,k),1,4);
    time_finish_MATLAB(k) = toc(time_start_MATLAB);

    %% 2 Ways to achive Registration -  Apply optimal transform to points

    % FIRST WAY: Rotate and translate P1(transformed) to P2(model)
    %     P3_SVD{k} = bsxfun(@plus,P1{k} *R_SVD{k},t_SVD{k});
    %     P3_CM{k} = bsxfun(@plus,P1{k} *R_CM{k},t_CM{k});

    % SECOND WAY: Rotate and translate P2 to P1
    % P3 = bsxfun(@minus,P2,t) * R';

    %% Convergences
    %% Angular Distance
    convergence_error_SVD(:,k) = 1/sqrt(2)*norm(logm(Rotation_GT(:,:,k) * R_SVD(:,:,k)'), 'fro');
    convergence_error_CM(:,k) = 1/sqrt(2)*norm(logm(Rotation_GT(:,:,k) * R_CM(:,:,k)'), 'fro');
    convergence_error_MATLAB(:,k) = 1/sqrt(2)*norm(logm(Rotation_GT(:,:,k) * R_MATLAB(:,:,k)'), 'fro');

    %% Chordal Distance
    convergence_error_SVD_chordal(:,k) = norm(Rotation_GT(:,:,k) - R_SVD(:,:,k) ,'fro');
    convergence_error_CM_chordal(:,k) = norm(Rotation_GT(:,:,k) - R_CM(:,:,k) ,'fro'); %
    convergence_error_MATLAB_chordal(:,k) = norm(Rotation_GT(:,:,k) - R_MATLAB(:,:,k) ,'fro'); 
end

%% 
%% FOR Plots for paper - LOAD workspace for each point cloud

%% RMSE SVD - CM
fig = figure;

subplot(3,2,1)
plot(RMSE_SVD_apple,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(RMSE_CM_apple,'LineWidth',1,'Color','r','LineStyle','none','Marker','+','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('SVD with μ: %.3e', mean(RMSE_SVD_apple));
legend_cm = sprintf('  CM with μ: %.3e', mean(RMSE_CM_apple));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.290629577453454 0.918330777747412 0.174963392677279 0.0806451590929163]);
grid on


subplot(3,2,2)
plot(RMSE_SVD_helix,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(RMSE_CM_helix,'LineWidth',1,'Color','#EDB120','LineStyle','none','Marker','+','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('SVD with μ: %.3e', mean(RMSE_SVD_helix));
legend_cm = sprintf('  CM with μ: %.3e', mean(RMSE_CM_helix));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.730600294876587 0.918330777747412 0.174963392677279 0.0806451590929163]);
grid on


subplot(3,2,3)
plot(RMSE_SVD_tree,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(RMSE_CM_tree,'LineWidth',1,'Color', '#00FF00','LineStyle','none','Marker','+','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('SVD with μ: %.3e', mean(RMSE_SVD_tree));
legend_cm = sprintf('  CM with μ: %.3e', mean(RMSE_CM_tree));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.289897513031785 0.6187916072405 0.174963392677279 0.0806451590929163]);
grid on


subplot(3,2,4)
plot(RMSE_SVD_horse ,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(RMSE_CM_horse ,'LineWidth',1,'Color', '#7E2F8E','LineStyle','none','Marker','+','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('SVD with μ: %.3e', mean(RMSE_SVD_horse));
legend_cm = sprintf('  CM with μ: %.3e', mean(RMSE_CM_horse));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.729136166033247 0.618791607240501 0.174963392677279 0.0806451590929163]);
grid on

subplot(3,2,5)
plot(RMSE_SVD_armadillo ,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(RMSE_CM_armadillo ,'LineWidth',1,'Color', '[1 0.5 0.8]','LineStyle','none','Marker','+','MarkerSize',9);
set(gca,'FontSize',15); % axes font size
legend_svd = sprintf('SVD with μ: %.3e', mean(RMSE_SVD_armadillo));
legend_cm = sprintf('  CM with μ: %.3e', mean(RMSE_CM_armadillo));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.289897513031784 0.319252436733588 0.174963392677279 0.0806451590929164]);
grid on

subplot(3,2,6)
plot(RMSE_SVD_buddha ,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(RMSE_CM_buddha ,'LineWidth',1,'Color', '#D95319','LineStyle','none','Marker','+','MarkerSize',9);
set(gca,'FontSize',15); % axes font size
legend_svd = sprintf('SVD with μ: %.3e', mean(RMSE_SVD_buddha));
legend_cm = sprintf('  CM with μ: %.3e', mean(RMSE_CM_buddha));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.729136166033247 0.319252436733588 0.174963392677279 0.0806451590929164]);
grid on


% Give common xlabel, ylabel and title to the figure
han=axes(fig,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'RMSE','FontSize', 18);
xlabel(han,'Test Case','FontSize', 18);


%% RMSE MATLAB - CM

figb = figure;

subplot(3,2,1)
plot(RMSE_MATLAB_apple,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(RMSE_CM_apple,'LineWidth',1,'Color','r','LineStyle','none','Marker','+','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('MATLAB with μ: %.3e', mean(RMSE_MATLAB_apple));
legend_cm = sprintf('         CM with μ: %.3e', mean(RMSE_CM_apple));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.263543194244407 0.918330777747412 0.201317711071947 0.0806451590929163]);
grid on


subplot(3,2,2)
plot(RMSE_MATLAB_helix,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(RMSE_CM_helix,'LineWidth',1,'Color','#EDB120','LineStyle','none','Marker','+','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('MATLAB with μ: %.3e', mean(RMSE_MATLAB_helix));
legend_cm = sprintf('         CM with μ: %.3e', mean(RMSE_CM_helix));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.70351391166754 0.918330777747412 0.201317711071947 0.0806451590929163]);
grid on


subplot(3,2,3)
plot(RMSE_MATLAB_tree,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(RMSE_CM_tree,'LineWidth',1,'Color', '#00FF00','LineStyle','none','Marker','+','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('MATLAB with μ: %.3e', mean(RMSE_MATLAB_tree));
legend_cm = sprintf('         CM with μ: %.3e', mean(RMSE_CM_tree));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.24926793825094 0.6187916072405 0.216691063468837 0.0806451590929163]);
grid on


subplot(3,2,4)
plot(RMSE_MATLAB_horse ,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(RMSE_CM_horse ,'LineWidth',1,'Color', '#7E2F8E','LineStyle','none','Marker','+','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('MATLAB with μ: %.3e', mean(RMSE_MATLAB_horse));
legend_cm = sprintf('         CM with μ: %.3e', mean(RMSE_CM_horse));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.703513911667538 0.620327705550793 0.201317711071947 0.0806451590929163]);
grid on

subplot(3,2,5)
plot(RMSE_MATLAB_armadillo ,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(RMSE_CM_armadillo ,'LineWidth',1,'Color', '[1 0.5 0.8]','LineStyle','none','Marker','+','MarkerSize',9);
set(gca,'FontSize',15); % axes font size
legend_svd = sprintf('MATLAB with μ: %.3e', mean(RMSE_MATLAB_armadillo));
legend_cm = sprintf('         CM with μ: %.3e', mean(RMSE_CM_armadillo));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.262811129822737 0.319252436733588 0.201317711071947 0.0806451590929164]);
grid on

subplot(3,2,6)
plot(RMSE_MATLAB_buddha ,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(RMSE_CM_buddha ,'LineWidth',1,'Color', '#D95319','LineStyle','none','Marker','+','MarkerSize',9);
set(gca,'FontSize',15); % axes font size
legend_svd = sprintf('MATLAB with μ: %.3e', mean(RMSE_MATLAB_buddha));
legend_cm = sprintf('         CM with μ: %.3e', mean(RMSE_CM_buddha));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.703513911667538 0.319252436733588 0.201317711071947 0.0806451590929164]);
grid on


% Give common xlabel, ylabel and title to the figure
han=axes(figb,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'RMSE','FontSize', 18);
xlabel(han,'Test Case','FontSize', 18);

%% HAUSDORFF SVD - CM

fig2 = figure;

subplot(3,2,1)
plot(haus_dist_SVD_apple,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(haus_dist_CM_apple,'LineWidth',1,'Color','r','LineStyle','none','Marker','x','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('SVD with μ: %.3e', mean(haus_dist_SVD_apple));
legend_cm = sprintf('  CM with μ: %.3e', mean(haus_dist_CM_apple));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.290629577453454 0.918330777747412 0.174963392677279 0.0806451590929163]);
grid on


subplot(3,2,2)
plot(haus_dist_SVD_helix,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(haus_dist_CM_helix,'LineWidth',1,'Color','#EDB120','LineStyle','none','Marker','x','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('SVD with μ: %.3e', mean(haus_dist_SVD_helix));
legend_cm = sprintf('  CM with μ: %.3e', mean(haus_dist_CM_helix));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.730600294876587 0.918330777747412 0.174963392677279 0.0806451590929163]);
grid on


subplot(3,2,3)
plot(haus_dist_SVD_tree,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(haus_dist_CM_tree,'LineWidth',1,'Color', '#00FF00','LineStyle','none','Marker','x','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('SVD with μ: %.3e', mean(haus_dist_SVD_tree));
legend_cm = sprintf('  CM with μ: %.3e', mean(haus_dist_CM_tree));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.289897513031785 0.6187916072405 0.174963392677279 0.0806451590929163]);
grid on


subplot(3,2,4)
plot(haus_dist_SVD_horse ,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(haus_dist_CM_horse ,'LineWidth',1,'Color', '#7E2F8E','LineStyle','none','Marker','x','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('SVD with μ: %.3e', mean(haus_dist_SVD_horse));
legend_cm = sprintf('  CM with μ: %.3e', mean(haus_dist_CM_horse));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.729136166033247 0.618791607240501 0.174963392677279 0.0806451590929163]);
grid on

subplot(3,2,5)
plot(haus_dist_SVD_armadillo ,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(haus_dist_CM_armadillo ,'LineWidth',1,'Color', '[1 0.5 0.8]','LineStyle','none','Marker','x','MarkerSize',9);
set(gca,'FontSize',15); % axes font size
legend_svd = sprintf('SVD with μ: %.3e', mean(haus_dist_SVD_armadillo));
legend_cm = sprintf('  CM with μ: %.3e', mean(haus_dist_CM_armadillo));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.289897513031784 0.319252436733588 0.174963392677279 0.0806451590929164]);
grid on

subplot(3,2,6)
plot(haus_dist_SVD_buddha ,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(haus_dist_CM_buddha,'LineWidth',1,'Color', '#D95319','LineStyle','none','Marker','x','MarkerSize',9);
set(gca,'FontSize',15); % axes font size
legend_svd = sprintf('SVD with μ: %.3e', mean(haus_dist_SVD_buddha));
legend_cm = sprintf('  CM with μ: %.3e', mean(haus_dist_CM_buddha));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.729136166033247 0.319252436733588 0.174963392677279 0.0806451590929164]);
grid on


% Give common xlabel, ylabel and title to the figure
han=axes(fig2,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Hausdorff Distance','FontSize', 18);
xlabel(han,'Test Case','FontSize', 18);


%% HAUSDORFF MATLAB - CM

fig2b = figure;

subplot(3,2,1)
plot(haus_dist_MATLAB_apple,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(haus_dist_CM_apple,'LineWidth',1,'Color','r','LineStyle','none','Marker','x','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('MATLAB with μ: %.3e', mean(haus_dist_MATLAB_apple));
legend_cm = sprintf('         CM with μ: %.3e', mean(haus_dist_CM_apple));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.263543194244407 0.918330777747412 0.201317711071947 0.0806451590929163]);
grid on


subplot(3,2,2)
plot(haus_dist_MATLAB_helix,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(haus_dist_CM_helix,'LineWidth',1,'Color','#EDB120','LineStyle','none','Marker','x','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('MATLAB with μ: %.3e', mean(haus_dist_MATLAB_helix));
legend_cm = sprintf('         CM with μ: %.3e', mean(haus_dist_CM_helix));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.70351391166754 0.918330777747412 0.201317711071947 0.0806451590929163]);
grid on


subplot(3,2,3)
plot(haus_dist_MATLAB_tree,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(haus_dist_CM_tree,'LineWidth',1,'Color', '#00FF00','LineStyle','none','Marker','x','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('MATLAB with μ: %.3e', mean(haus_dist_MATLAB_tree));
legend_cm = sprintf('         CM with μ: %.3e', mean(haus_dist_CM_tree));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.24926793825094 0.6187916072405 0.216691063468837 0.0806451590929163]);
grid on


subplot(3,2,4)
plot(haus_dist_MATLAB_horse ,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(haus_dist_CM_horse ,'LineWidth',1,'Color', '#7E2F8E','LineStyle','none','Marker','x','MarkerSize',9);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('MATLAB with μ: %.3e', mean(haus_dist_MATLAB_horse));
legend_cm = sprintf('         CM with μ: %.3e', mean(haus_dist_CM_horse));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.703513911667538 0.620327705550793 0.201317711071947 0.0806451590929163]);
grid on

subplot(3,2,5)
plot(haus_dist_MATLAB_armadillo ,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(haus_dist_CM_armadillo ,'LineWidth',1,'Color', '[1 0.5 0.8]','LineStyle','none','Marker','x','MarkerSize',9);
set(gca,'FontSize',15); % axes font size
legend_svd = sprintf('MATLAB with μ: %.3e', mean(haus_dist_MATLAB_armadillo));
legend_cm = sprintf('         CM with μ: %.3e', mean(haus_dist_CM_armadillo));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.262811129822737 0.319252436733588 0.201317711071947 0.0806451590929164]);
grid on

subplot(3,2,6)
plot(haus_dist_MATLAB_buddha ,'LineWidth',1,'Color','k','LineStyle','none','Marker','o','MarkerSize',9);
hold on
plot(haus_dist_CM_buddha ,'LineWidth',1,'Color', '#D95319','LineStyle','none','Marker','x','MarkerSize',9);
set(gca,'FontSize',15); % axes font size
legend_svd = sprintf('MATLAB with μ: %.3e', mean(haus_dist_MATLAB_buddha));
legend_cm = sprintf('         CM with μ: %.3e', mean(haus_dist_CM_buddha));
legend(legend_svd,legend_cm, 'Fontsize',14, 'Position',[0.703513911667538 0.319252436733588 0.201317711071947 0.0806451590929164]);
grid on


% Give common xlabel, ylabel and title to the figure
han=axes(fig2b,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Hausdorff Distance','FontSize', 18);
xlabel(han,'Test Case','FontSize', 18);


%% Angular Distance PLOTS
svd_ang_apple = convergence_error_SVD_apple;
cm_ang_apple = convergence_error_CM_apple;
MATLAB_ang_apple = convergence_error_MATLAB_apple;

svd_ang_helix = convergence_error_SVD_helix;
cm_ang_helix = convergence_error_CM_helix;
MATLAB_ang_helix = convergence_error_MATLAB_helix;

svd_ang_tree = convergence_error_SVD_tree;
cm_ang_tree = convergence_error_CM_tree;
MATLAB_ang_tree = convergence_error_MATLAB_tree;

svd_ang_horse = convergence_error_SVD_horse;
cm_ang_horse = convergence_error_CM_horse;
MATLAB_ang_horse = convergence_error_MATLAB_horse;

svd_ang_armadillo = convergence_error_SVD_armadillo;
cm_ang_armadillo = convergence_error_CM_armadillo;
MATLAB_ang_armadillo = convergence_error_MATLAB_armadillo;

svd_ang_buddha = convergence_error_SVD_buddha;
cm_ang_buddha = convergence_error_CM_buddha;
MATLAB_ang_buddha = convergence_error_MATLAB_buddha;

%% ALL
fig3 = figure;

subplot(3,2,1)
plot(svd_ang_apple,'LineWidth',1.8,'Color','k','LineStyle','-');
hold on
plot(cm_ang_apple,'LineWidth',2.5,'Color','r','LineStyle','--');
hold on
plot(MATLAB_ang_apple,'LineWidth',2.5,'Color','#4DBEEE','LineStyle',':');
% xlabel('No Test','FontSize', 23);
% ylabel('RMSE','FontSize', 23);
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('μ(SVD): %.1e \n', mean(svd_ang_apple));
legend_cm = sprintf('μ(CM): %.1e', mean(cm_ang_apple));
legend_MATLAB = sprintf('μ(MATLAB): %.1e', mean(MATLAB_ang_apple));
legend(legend_svd, legend_cm, legend_MATLAB, 'Fontsize',14, 'Position',[0.173865303506264 0.906041991265078 0.30014640616463 0.0806451590929163],...
    'Orientation','horizontal','NumColumns',2);
grid on
legend boxoff


subplot(3,2,2)
plot(svd_ang_helix,'LineWidth',1.8,'Color','k','LineStyle','-');
hold on
plot(cm_ang_helix,'LineWidth',2.5,'Color','#EDB120','LineStyle','--')
hold on
plot(MATLAB_ang_helix,'LineWidth',2.5,'Color','#4DBEEE','LineStyle',':');
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('μ(SVD): %.1e', mean(svd_ang_helix));
legend_cm = sprintf('μ(CM): %.1e', mean(cm_ang_helix));
legend_MATLAB = sprintf('μ(MATLAB): %.1e', mean(MATLAB_ang_helix));
legend(legend_svd, legend_cm, legend_MATLAB, 'Fontsize',14, 'Position',[0.615300149772736 0.90757808957537 0.30014640616463 0.0806451590929163],...
    'Orientation','horizontal','NumColumns',2);
grid on
legend boxoff


subplot(3,2,3)
plot(svd_ang_tree,'LineWidth',1.8,'Color','k','LineStyle','-');
hold on
plot(cm_ang_tree,'LineWidth',2.5,'Color', '#00FF00','LineStyle','--')
hold on
plot(MATLAB_ang_tree,'LineWidth',2.5,'Color','#4DBEEE','LineStyle',':');
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('μ(SVD): %.1e', mean(svd_ang_tree));
legend_cm = sprintf('μ(CM): %.1e', mean(cm_ang_tree));
legend_MATLAB = sprintf('μ(MATLAB): %.1e', mean(MATLAB_ang_tree));
legend(legend_svd, legend_cm, legend_MATLAB, 'Fontsize',14, 'Position',[0.173865303506264 0.606502820758166 0.30014640616463 0.0806451590929163],...
    'Orientation','horizontal','NumColumns',2);
grid on
legend boxoff

subplot(3,2,4)
plot(svd_ang_horse ,'LineWidth',1.8,'Color','k','LineStyle','-');
hold on
plot(cm_ang_horse ,'LineWidth',2.5,'Color', '#7E2F8E','LineStyle','--')
hold on
plot(MATLAB_ang_horse,'LineWidth',2.5,'Color','#4DBEEE','LineStyle',':');
set(gca,'FontSize',15,'xtick',[]); % axes font size
legend_svd = sprintf('μ(SVD): %.1e', mean(svd_ang_horse));
legend_cm = sprintf('μ(CM): %.1e', mean(cm_ang_horse));
legend_MATLAB = sprintf('μ(MATLAB): %.1e', mean(MATLAB_ang_horse));
legend(legend_svd, legend_cm, legend_MATLAB, 'Fontsize',14, 'Position',[0.614568085351067 0.606502820758166 0.30014640616463 0.0806451590929163],...
    'Orientation','horizontal','NumColumns',2);
grid on
legend boxoff

subplot(3,2,5)
plot(svd_ang_armadillo ,'LineWidth',1.8,'Color','k','LineStyle','-');
hold on
plot(cm_ang_armadillo ,'LineWidth',2.5,'Color', '[1 0.5 0.8]','LineStyle','--')
hold on
plot(MATLAB_ang_armadillo,'LineWidth',2.5,'Color','#4DBEEE','LineStyle',':');
set(gca,'FontSize',15); % axes font size
legend_svd = sprintf('μ(SVD): %.1e', mean(svd_ang_armadillo));
legend_cm = sprintf('μ(CM): %.1e', mean(cm_ang_armadillo));
legend_MATLAB = sprintf('μ(MATLAB): %.1e', mean(MATLAB_ang_armadillo));
legend(legend_svd, legend_cm, legend_MATLAB, 'Fontsize',14, 'Position',[0.173865303506264 0.303891453630669 0.30014640616463 0.0806451590929162],...
    'Orientation','horizontal','NumColumns',2);
grid on
legend boxoff

subplot(3,2,6)
plot(svd_ang_buddha ,'LineWidth',1.8,'Color','k','LineStyle','-');
hold on
plot(cm_ang_buddha ,'LineWidth',2.5,'Color', '#D95319','LineStyle','--')
hold on
plot(MATLAB_ang_buddha,'LineWidth',2.5,'Color','#4DBEEE','LineStyle',':');
set(gca,'FontSize',15); % axes font size
legend_svd = sprintf('μ(SVD): %.1e', mean(svd_ang_buddha));
legend_cm = sprintf('μ(CM): %.1e', mean(cm_ang_buddha));
legend_MATLAB = sprintf('μ(MATLAB): %.1e', mean(MATLAB_ang_buddha));
legend(legend_svd, legend_cm, legend_MATLAB, 'Fontsize',14, 'Position',[0.614568085351067 0.303891453630669 0.30014640616463 0.0806451590929163],...
    'Orientation','horizontal','NumColumns',2);
grid on
legend boxoff

% Give common xlabel, ylabel and title to the figure
han=axes(fig3,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Radians','FontSize', 18);
xlabel(han,'Test Case','FontSize', 18);


%% TIMES
fig4 = figure;

subplot(3,2,1)
plot(time_finish_SVD_apple,'LineWidth',2.5,'Color','k','LineStyle','-');
hold on
plot(time_finish_CM_apple,'LineWidth',2.5,'Color','r','LineStyle',':');
hold on
plot(time_finish_MATLAB_apple,'LineWidth',2.5,'Color','#A2142F','LineStyle','-.');
set(gca,'FontSize',14,'xtick',[]); % axes font size
grid on

subplot(3,2,2)
plot(time_finish_SVD_helix,'LineWidth',2.5,'Color','k','LineStyle','-');
hold on
plot(time_finish_CM_helix,'LineWidth',2.5,'Color','#EDB120','LineStyle',':');
hold on
plot(time_finish_MATLAB_helix,'LineWidth',2.5,'Color','#A2142F','LineStyle','-.');
set(gca,'FontSize',14,'xtick',[]); % axes font size
grid on


subplot(3,2,3)
plot(time_finish_SVD_tree,'LineWidth',2.5,'Color','k','LineStyle','-');
hold on
plot(time_finish_CM_tree,'LineWidth',2.5,'Color','#00FF00','LineStyle',':');
hold on
plot(time_finish_MATLAB_tree,'LineWidth',2.5,'Color','#A2142F','LineStyle','-.');
set(gca,'FontSize',14,'xtick',[]); % axes font size
grid on


subplot(3,2,4)
plot(time_finish_SVD_horse,'LineWidth',2.5,'Color','k','LineStyle','-');
hold on
plot(time_finish_CM_horse,'LineWidth',2.5,'Color','#7E2F8E','LineStyle',':');
hold on
plot(time_finish_MATLAB_horse,'LineWidth',2.5,'Color','#A2142F','LineStyle','-.');
set(gca,'FontSize',14,'xtick',[]); % axes font size
grid on

subplot(3,2,5)
plot(time_finish_SVD_armadillo,'LineWidth',2.5,'Color','k','LineStyle','-');
hold on
plot(time_finish_CM_armadillo,'LineWidth',2.5,'Color','[1 0.5 0.8]','LineStyle',':');
hold on
plot(time_finish_MATLAB_armadillo,'LineWidth',2.5,'Color','#A2142F','LineStyle','-.');
set(gca,'FontSize',14); % axes font size
grid on


subplot(3,2,6)
plot(time_finish_SVD_buddha,'LineWidth',2.5,'Color','k','LineStyle','-');
hold on
plot(time_finish_CM_buddha,'LineWidth',2.5,'Color','#D95319','LineStyle',':');
hold on
plot(time_finish_MATLAB_buddha,'LineWidth',2.5,'Color','#A2142F','LineStyle','-.');
set(gca,'FontSize',14); % axes font size
grid on

% Modify legend (Create our legend)
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'-k',LineWidth= 3);
h(2) = plot(NaN,NaN,':k', LineWidth= 3);
h(3) = plot(NaN,NaN,'-.', LineWidth= 3, Color ='#A2142F');
legend(h, 'SVD', 'CM', 'MATLAB', 'Fontsize',14, 'Orientation','horizontal',...
    'Position',[0.399219134425232 0.953158321532331 0.217423130094045 0.0430107515894689]);

% Give common xlabel, ylabel and title to the figure
han=axes(fig4,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylh = ylabel(han,'Time (s)','FontSize', 20);
xlabel(han,'Test Case','FontSize', 20);

% Move Ylabel as it coincides with the axis
ylh.Position(2) = ylh.Position(2) - 0.0000000008;

