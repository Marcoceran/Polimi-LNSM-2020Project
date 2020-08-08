user = [0 0];               %those are the coordinates of our reference point,
                            %since they are both zero we can ignore them in
                            %the expressions used to compute the positions
                            %of the access points, since they're additive
                            %in that expressions
positions = zeros(8, 2);    %this is the vector of positions, expressed as integers coordinates on an xy plane
for i = 1:8
    hTOA = rhoUEAP(i, 1);   %this is the result of the TOA position computation
    hAOA = rhoUEAP(i, 2)/2; %this is the result of the AOA position computation
    positions(i, 2) = ((hTOA^2*(tan(hAOA))^2)/(1+(tan(hAOA))^2))^0.5;   %this expression computes the y-coordinate of each AP
    positions(i, 1) = round((hTOA^2-positions(i, 2).^2)^0.5);           %this expression computes the x-coordinate of each AP
    positions(i, 2) = round(positions(i, 2));                           %both coordinates have to be rounded to nearest integer
end