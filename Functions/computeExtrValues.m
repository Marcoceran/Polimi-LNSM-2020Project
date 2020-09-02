function MinMaxUV = computeExtrValues(UEtraj)

startU = zeros(2, length(UEtraj));
startV = zeros(2, length(UEtraj));

for j = 1:size(UEtraj, 2)
    
    startV(:,j) = UEtraj{j}(1, 3:4);
    startU(:,j) = UEtraj{j}(1, 1:2);
    
end

MinU = round(min(startU, [], 2), -3);
MaxU = round(max(startU, [], 2), -3);

MinV = round(min(startV, [], 2));
MaxV = round(max(startV, [], 2));

MinMaxUV = [MinU(1), MaxU(1), MinV(1), MaxV(1)];

end