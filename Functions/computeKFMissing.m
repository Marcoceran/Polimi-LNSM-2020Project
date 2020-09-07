function x_hatU = computeKFMissing(TST, rhoTraj, PosOfAP, F, Q, R)

sigma_upP = 1;

x_hatU    = zeros(4, TST);
update_P  = cell (1, TST);
x_hatP    = zeros(4, TST);
predict_P = cell (1, TST);

x_hatU(1:2, 1) = 1000*randn(2, 1);
x_hatU(3:4, 1) = 10*randn(2, 1);

update_P{1}  = sigma_upP^2 * eye(4);

for t = 2:TST
    
    %prediction
    x_hatP(:, t) = F * x_hatU(:, t - 1);
    predict_P{t} = F .* update_P{t - 1} .* transpose(F) + Q;
    
    %magic for AP
    N = ~isnan(rhoTraj(t, :));
    rhoT = rhoTraj(t, :);
    rhoT(isnan(rhoT)) = 0;
    N_AP = nnz(N);
    AP_pos = nonzeros(N' .* PosOfAP)';
    AP_pos = [AP_pos(1:N_AP);
              AP_pos(N_AP + 1:end)];  
    R_AP = nonzeros(N.' .* R) .* eye(N_AP);
    rhoT = nonzeros(N.* rhoT);
    
    %update
    h = sqrt( sum( (AP_pos.' - x_hatP(1:2, t).').^2, 2) );
    H = calculateH(x_hatP(:, t), AP_pos.', N_AP);
    G = predict_P{t} * H.' * inv(H * predict_P{t} * H.' + R_AP);
    
    x_hatU(:, t) = x_hatP(:,t) + G * (rhoT(:) - h);
    update_P{t}  = predict_P{t} - G * H * predict_P{t};
    
end

end