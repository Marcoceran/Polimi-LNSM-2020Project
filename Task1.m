%% TASK - 1a %%%%%%%%%%
clear all, clc, close all
set(0,'DefaultTextFontSize',18)
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',16)

%% scenario settings (4000x4000m)
parameters.xmin = -2000; parameters.ymin = -2000;
parameters.xmax =  2000; parameters.ymax =  2000;

%% user position
UE_a = [0, 0];

%% AP position from measurement 
load('Task1a_rhoUEAP.mat')
rhoA_TOA = rhoUEAP(:,1);
rhoA_AOA = rhoUEAP(:,2);

positionOfAP = zeros(size(rhoUEAP, 1),2);
parameters.NumOfAP = size(rhoUEAP,1);

for a = 1:parameters.NumOfAP
    z = rhoA_TOA(a) * exp(1i * rhoA_AOA(a));
   
    positionOfAP(a,1) = round(UE_a(1) - real(z));
    positionOfAP(a,2) = round(UE_a(2) - imag(z));
end

%% plot 2D
figure (1)
plot(positionOfAP(:,1), positionOfAP(:,2), '^','MarkerSize', 10); hold on;
plot(UE_a(:,1), UE_a(:,2), 'o');
xlabel('[m]'), ylabel('[m]');
xlim([parameters.xmin parameters.xmax])
ylim([parameters.ymin parameters.ymax])
grid on;
axis equal   


%% TASK - 1b %%%%%%%%%%
%% user position
UE_b = [500, -800];

%% initialization
load('Task1b_rhoUEAP');
TOA_meas = rhoUEAP;

errTOA = zeros(size(TOA_meas, 1), parameters.NumOfAP);
StdMatrix = zeros(parameters.NumOfAP, parameters.NumOfAP);

%AP are in the same position of the task-1a
distanceUEAP = sqrt ( sum([UE_b - positionOfAP].^2, 2));

%% Computing the accuracy matrix
for b = 1:parameters.NumOfAP
    
    for a = 1:size(rhoUEAP,1)
        errTOA(a,b) = abs(distanceUEAP(b) - TOA_meas(a,b)); 
    end
    
    StdMatrix(b,b) = round(mean(errTOA(:,b)));
    
end

CovarianceMatrix = StdMatrix.^2;
