function posOfAP = getPositionOfAP(NumOfAP, rhoUEAP, UE)
% NumOfAP is the number of access points, rhoUEAP is the vector that contains the TOA and AOA
% measurements of the user equipment whose coordinates are stored in the vector UE, which is [0 0] in this case
rhoTOA = rhoUEAP(:,1);
rhoAOA = rhoUEAP(:,2);

posOfAP = zeros(size(rhoUEAP, 1),2);

for a = 1:NumOfAP
    z = rhoTOA(a) * exp(1i * rhoAOA(a));
    % we combine the angle and the magnitude measurements in the vector to form complex numbers
    % that we subsequently subtract from the UE position
    posOfAP(a,1) = round(UE(1) - real(z));
    posOfAP(a,2) = round(UE(2) - imag(z));
    % and we get the x and y coordinates as the real and imaginary part of those complex numbers
    % the positions of the APs have to be integer numbers
end

end
