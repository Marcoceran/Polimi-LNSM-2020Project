function CovMat = computeCovMat(NumOfAP, rhoUEAP, UE, PosOfAP)

errTOA = zeros(size(rhoUEAP, 1), NumOfAP);
StdMatrix = zeros(NumOfAP, NumOfAP);

%AP are in the same position of the task-1a
distanceUEAP = sqrt ( sum((PosOfAP - UE).^2, 2));

for b = 1:NumOfAP
    
    for a = 1:size(rhoUEAP,1)
        errTOA(a,b) = abs(distanceUEAP(b) - rhoUEAP(a,b)); 
    end
    
    StdMatrix(b,b) = round(mean(errTOA(:,b)));
    
end

CovMat = StdMatrix.^2;
%CovMat = StdMatrix;
end