function posOfAP = getPositionOfAP(NumOfAP, rhoUEAP, UE)

rhoTOA = rhoUEAP(:,1);
rhoAOA = rhoUEAP(:,2);

posOfAP = zeros(size(rhoUEAP, 1),2);

for a = 1:NumOfAP
    z = rhoTOA(a) * exp(1i * rhoAOA(a));
   
    posOfAP(a,1) = round(UE(1) - real(z));
    posOfAP(a,2) = round(UE(2) - imag(z));
end

end