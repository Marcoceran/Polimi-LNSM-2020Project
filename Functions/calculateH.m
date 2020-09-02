function H = calculateH(UE, PosOfAP, NumOfAP)

distanceUEAP = sqrt ( sum((PosOfAP - UE([1,2])').^2, 2));
H = zeros(NumOfAP, 4);

for a = 1:NumOfAP
    
    H(a,1) = -(PosOfAP(a, 1)- UE(1,1))/distanceUEAP(a);
    H(a,2) = -(PosOfAP(a, 2)- UE(2,1))/distanceUEAP(a);
    
end

end