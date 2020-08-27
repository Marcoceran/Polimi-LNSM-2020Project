%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK - 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, clc, close all

%% scenario settings (4000x4000m)
parameters.xmin = -2000; parameters.ymin = -2000;
parameters.xmax =  2000; parameters.ymax =  2000;
TotalSimulationTime = 200; %s
Ts = 1; %s
load("Task2_trajectory_GR12");

%% plot trajectory and velocity 

for a = 1:size(UEtrajectory, 2)
    figure(1)
    title('Trajectory')
    plot(UEtrajectory{a}(:, 1), UEtrajectory{a}(:, 2), '.');
    xlabel('[m]'), ylabel('[m]');
    xlim([parameters.xmin parameters.xmax])
    ylim([parameters.ymin parameters.ymax])
    txt_Start = " START"; text(UEtrajectory{a}(1, 1),UEtrajectory{a}(1, 2), txt_Start)
    txt_End = " END"; text(UEtrajectory{a}(end, 1),UEtrajectory{a}(end, 2), txt_End)
    
    
%     figure(2); hold on
%     title('Velocity')
%     plot(UEtrajectory{a}(:, 3));plot(UEtrajectory{a}(:, 4));
%     xlabel('[s]'), ylabel('[m/s]');
%     xlim([0 TotalSimulationTime])
    
end

%% calculate model parameters and statistics
Acceleration = cell(1, size(UEtrajectory, 2));

startU = zeros(2, length(UEtrajectory));
startV = zeros(2, length(UEtrajectory));

stdAcc = zeros(2, length(UEtrajectory));
meanAcc = zeros(2, length(UEtrajectory));
%wrongTraj = zeros (1, length(UEtrajectory));

% Extract acceleration, starting position(u0) and starting 
% velocity (v0) from each trajectory
for j = 1:size(UEtrajectory, 2)
    
    Acceleration{j} = diff(UEtrajectory{j}(:, [3, 4]), 1, 1);
    
    startV(:,j) = UEtrajectory{j}(1, 3:4);
    startU(:,j) = UEtrajectory{j}(1, 1:2);
    
end
%startV = sort(startV, 2);
%startU = sort(startU, 2);

% Select trajectories that bounce back and eliminate them
for a = 1:size(Acceleration, 2)
    for i = 1:size(Acceleration{1, 1}, 1)
        for j = 1:size(Acceleration{1, 1}, 2)
        
            if(abs(Acceleration{a}(i,j)) >1)
                wrongTraj(a) = a;
            end
        
        end
    end
end
wrongTraj = nonzeros(wrongTraj)';

for a = 1:length(wrongTraj) 
    
    Acceleration{wrongTraj(a)} = 0;
   
end

% Compute the standard deviation and the mean of the 
% acceleration of each trajectory
for j = 1:size(UEtrajectory, 2)
    
        stdAcc(:,j)  = std(Acceleration{j}, 0, 1);   
        meanAcc(:,j) = mean(Acceleration{j}, 1);
end
stdAcc ( :, ~any(stdAcc, 1) ) = [];
meanAcc ( :, ~any(meanAcc, 1) ) = [];


% Parameter & statistic values
meanStdA = round(mean(stdAcc, 2),2);
meanMeA = round(mean(meanAcc, 2),3);

maxStartingU = round(max(startU, [], 2)', -3);
minStartingU = round(min(startU, [], 2)', -3);

maxStartingV = round(max(startV, [], 2)');
minStartingV = round(min(startV, [], 2)');


% for j = 1:size(Acceleration, 2)
%     figure(2); hold on
%     title('Velocity')
%     plot(diff(UEtrajectory{j}(:, 1) + Ts*UEtrajectory{j}(:, 3), 1, 1));
%     plot(diff(UEtrajectory{j}(:, 2) + Ts*UEtrajectory{j}(:, 4), 1, 1));
%     xlabel('[s]'), ylabel('[m/s]');
%     xlim([0 TotalSimulationTime])
%     
% %     figure(3); hold on; grid on
% %     title('Acceleration')
% %     plot((Acceleration{j}(:,1)));
% %     plot((Acceleration{j}(:,2)));
% %     xlabel('[s]'), ylabel('[m/s^2]');
% %     xlim([0 TotalSimulationTime])
%     
% end

