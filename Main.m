%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, clc, close all
set(0,'DefaultLineLineWidth',0.2);

%% scenario settings (4000x4000m)
parameters.xmin = -2000; parameters.ymin = -2000;
parameters.xmax =  2000; parameters.ymax =  2000;

%% TASK - 1a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UE_a = [0, 0];

load('Task1a_rhoUEAP')

parameters.NumOfAP = size(rhoUEAP,1);
parameters.PosOfAP = getPositionOfAP(parameters.NumOfAP, rhoUEAP, UE_a);

figure (1)
plot(parameters.PosOfAP(:,1), parameters.PosOfAP(:,2), '^','MarkerSize', 10); hold on;
plot(UE_a(:,1), UE_a(:,2), 'o');
xlabel('[m]'), ylabel('[m]');
xlim([parameters.xmin parameters.xmax])
ylim([parameters.ymin parameters.ymax])
grid on;
axis equal 

%% TASK - 1b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UE_b = [500, -800];

load('Task1b_rhoUEAP')

CovMatrix = computeCovMat(parameters.NumOfAP, rhoUEAP, UE_b, parameters.PosOfAP);

%% TASK - 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("Task2_trajectory_GR12");

[wrongT, parameters.StdAcc] = computeStdAcceleration(UEtrajectory);
parameters.MinMaxUV = computeExtrValues(UEtrajectory);

%% TASK - 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("Task3_rhoUEAP_GR12");
%parameters
TotalSimulationTime = 200; %s
Ts = 1; %s
sigma_a   = parameters.StdAcc;
sigma_upP = 0.5;

%model matrices
F = [eye(2)     , Ts*eye(2);
     zeros(2,2) ,    eye(2)];
L = [0.5*Ts^2*eye(2); Ts*eye(2)];
Q = sigma_a^2 .* (L * transpose(L));
R = CovMatrix;

% initialization
x_hat_3    = cell(1, length(UEtrajectory)); 

for Traj = 1:100
    
    x_hat_3{Traj}(:, 1) = UEtrajectory{Traj}(1, :)';
    
    x_hat_3{Traj} = computeKFTraj(x_hat_3{Traj}, TotalSimulationTime,  rhoUEEAP{Traj},...
                                     sigma_upP, parameters.NumOfAP, parameters.PosOfAP, F, Q, R);
    
%     figure(2)
%     title('Trajectory')
%     plot(UEtrajectory{Traj}(:, 1), UEtrajectory{Traj}(:, 2), '-o');hold on;
%     plot(parameters.PosOfAP(:,1), parameters.PosOfAP(:,2), '^','MarkerSize', 10);
%     plot(x_hat_3{Traj}(1,:), x_hat_3{Traj}(2,:), '-^');
%     xlabel('[m]'), ylabel('[m]');
%     xlim([parameters.xmin parameters.xmax])
%     ylim([parameters.ymin parameters.ymax])
%     txt_Start = " START"; text(UEtrajectory{Traj}(1, 1),UEtrajectory{Traj}(1, 2), txt_Start)
%     txt_End   = " END"; text(UEtrajectory{Traj}(end, 1),UEtrajectory{Traj}(end, 2), txt_End)
 
end    

% evaluate performances
MSE = zeros(1,length(UEtrajectory));

for a = 1:length(UEtrajectory)
     if a~= wrongT
        MSE(a) = mean(sqrt(sum((x_hat_3{a}([1,2],:) - UEtrajectory{a}(:,[1,2]).').^2, 1)), 2);
     else
         MSE(a) = 0;
     end
     
%      MSE(a) = mean(sqrt(sum((x_hat{a}([1,2],:) - UEtrajectory{a}(:,[1,2]).').^2, 1)), 2);

end
MSE = nonzeros(MSE);
RMSE = sqrt(mean(MSE));

%% TASK - 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Task4_rhoUEAP_GR12.mat')

TotalSimulationTime = length(rhoUEEAP);
x_hat_4    = zeros(4, TotalSimulationTime);

for a =1:4
    
    x_hat_4(1, 1) = -1000 + 100*randn;
    x_hat_4(2, 1) = 1000 + 100*randn;
    x_hat_4(3, 1) = 13.9;
    x_hat_4(4, 1) = 0;
    
    x_hat_4 = computeKFTraj(x_hat_4, TotalSimulationTime,  rhoUEEAP,...
        sigma_upP, parameters.NumOfAP, parameters.PosOfAP, F, Q, R);
    
    figure(3), hold on;
    plot(x_hat_4(1,:), x_hat_4(2,:), '-o');
    plot(parameters.PosOfAP(:,1), parameters.PosOfAP(:,2), '^','MarkerSize', 10);
    txt_Start = " START"; text(x_hat_4(1, 1),x_hat_4(2, 1), txt_Start)
    txt_End   = " END"; text(x_hat_4(1, end),x_hat_4(2, end), txt_End)
    
end

%% TASK - 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Task5_rhoUEAP_GR12.mat')

TotalSimulationTime = length(rhoUEEAP);
x_hat_5    = zeros(4, TotalSimulationTime);

for a =1:4
    
    x_hat_5(1:2, 1) = 1000*randn(2, 1);
    x_hat_5(3:4, 1) = 10*randn(2, 1);
    
    x_hat_5 = computeKFTraj(x_hat_5, TotalSimulationTime,  rhoUEEAP,...
        sigma_upP, parameters.NumOfAP, parameters.PosOfAP, F, Q, R);
    
    figure(4), hold on;
    plot(x_hat_5(1,:), x_hat_5(2,:), '-o');
    plot(parameters.PosOfAP(:,1), parameters.PosOfAP(:,2), '^','MarkerSize', 10);
    txt_Start = " START"; text(x_hat_5(1, 1),x_hat_5(2, 1), txt_Start)
    txt_End   = " END"; text(x_hat_5(1, end),x_hat_5(2, end), txt_End)
    
end

%% TASK - 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Task6_rhoUEAP_GR12.mat')

TotalSimulationTime = length(rhoUEEAP);
for a =1:4
    x_traj = computeKFMissing(TotalSimulationTime,  rhoUEEAP,...
        parameters.PosOfAP, F, Q, R);
    
    figure(5), hold on;
    plot(x_traj(1,:), x_traj(2,:), '-o');
    plot(parameters.PosOfAP(:,1), parameters.PosOfAP(:,2), '^','MarkerSize', 10);
    txt_Start = " START"; text(x_traj(1, 1),x_traj(2, 1), txt_Start)
    txt_End   = " END"; text(x_traj(1, end),x_traj(2, end), txt_End)
    xlim([parameters.xmin parameters.xmax])
    ylim([parameters.ymin parameters.ymax])

end
