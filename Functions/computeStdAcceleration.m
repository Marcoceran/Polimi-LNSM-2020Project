function meanStdAcc = computeStdAcceleration(UEtraj)

Acceleration = cell(1, length(UEtraj));
stdAcc = zeros(2, length(UEtraj));

for j = 1:size(UEtraj, 2)
    
    Acceleration{j} = diff(UEtraj{j}(:, [3, 4]), 1, 1);
      
end

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

for j = 1:size(Acceleration, 2)
    
        stdAcc(:,j)  = std(Acceleration{j}, 0, 1);   
        
end
stdAcc ( :, ~any(stdAcc, 1) ) = [];

meanSigma = round(mean(stdAcc, 2),2);
meanStdAcc = meanSigma(1);
end