function x_hatU = computeKFTraj(x_hatU, TST, rhoTraj, sigma_upP, NumOfAP, PosOfAP, F, Q, R)

update_P  = cell (1, TST);
x_hatP    = zeros(4, TST);
predict_P = cell (1, TST);

update_P{1}  = sigma_upP^2 * eye(4);

for t = 2:TST
    
    %prediction
    x_hatP(:, t) = F * x_hatU(:, t - 1);
    predict_P{t} = F .* update_P{t - 1} .* transpose(F) + Q;
    
    %update
    h = sqrt( sum( (PosOfAP - x_hatP(1:2, t).').^2, 2) );
    H = calculateH(x_hatP(:, t), PosOfAP, NumOfAP);
    G = predict_P{t} * H.' * inv(H * predict_P{t} * H.' + R);
    
    x_hatU(:, t) = x_hatP(:,t) + G * (rhoTraj(t,:).' - h);
    update_P{t}  = predict_P{t} - G * H * predict_P{t};
    
end

end