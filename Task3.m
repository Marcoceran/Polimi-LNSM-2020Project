%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK - 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

parameters.StdAcc = computeStdAcceleration(UEtrajectory);
parameters.MinMaxUV = computeExtrValues(UEtrajectory);

%% TASK - 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("Task3_rhoUEAP_GR12");
%parameters
TotalSimulationTime = 200; %s
Ts = 1; %s
sigma_a   = parameters.StdAcc;
sigma_upP = 10000;
%model matrices
F = [eye(2)     , Ts*eye(2);
     zeros(2,2) ,    eye(2)];
L = [0.5*Ts^2*eye(2); Ts*eye(2)];
Q = sigma_a^2 .* (L * transpose(L));
R = CovMatrix;

% initialization
x_hatU    = zeros(4, TotalSimulationTime);
update_P  = cell (1, TotalSimulationTime);
x_hatP    = zeros(4, TotalSimulationTime);
predict_P = cell (1, TotalSimulationTime);

update_P{1}  = sigma_upP^2 * eye(4);

for Traj = 1:100
    
    x_hatU(:, 1) = UEtrajectory{Traj}(1, :)';
    for t = 2:TotalSimulationTime
   
        %prediction
        x_hatP(:, t) = F * x_hatU(:, t - 1);
        predict_P{t} = F .* update_P{t - 1} .* transpose(F) + Q;
    
        %update
        h = sqrt( sum( (parameters.PosOfAP - x_hatP(1:2, t).').^2, 2) );
        H = calculateH(x_hatP(:, t), parameters.PosOfAP, parameters.NumOfAP);
        G = predict_P{t} * H.' * inv(H * predict_P{t} * H.' + R);
    
        x_hatU(:, t) = x_hatP(:,t) + G * (rhoUEEAP{Traj}(t,:).' - h);
        update_P{t}  = predict_P{t} - G * H * predict_P{t};
    
    end

    figure(2)
    title('Trajectory')
    plot(UEtrajectory{Traj}(:, 1), UEtrajectory{Traj}(:, 2), '-o');hold on;
    plot(x_hatU(1,:), x_hatU(2,:), '-^');
    xlabel('[m]'), ylabel('[m]');
    xlim([parameters.xmin parameters.xmax])
    ylim([parameters.ymin parameters.ymax])
    txt_Start = " START"; text(UEtrajectory{Traj}(1, 1),UEtrajectory{Traj}(1, 2), txt_Start)
    txt_End   = " END"; text(UEtrajectory{Traj}(end, 1),UEtrajectory{Traj}(end, 2), txt_End)
    
end    