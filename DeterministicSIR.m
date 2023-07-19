function [rmse, Y] = DeterministicSIR(S0, r, a, actual_data)
    S = []; I = []; R = [];
    S(1) = S0; I(1) = actual_data(1, 1);
    R(1) = actual_data(1, 2);
    
    for n = 1:size(actual_data, 1)-1
        S(n+1) = abs(S(n) - r*S(n)*I(n));
        I(n+1) = abs(I(n) + (r*S(n)*I(n) - a*I(n)));
        R(n+1) = abs(R(n) + a*I(n));
    end
    
    Y = [S' I' R'];
    
    % Compute error
    rmse = mean(mean(sqrt((Y(:,2:3)-actual_data).^2)));
%     mape_r = mean(abs(((actual_data(:,1)-Y(:,2))./actual_data(:,1))));
%     mape_a = mean(abs(((actual_data(:,2)-Y(:,3))./actual_data(:,2))));
%     mape = mean(mean(mape_r+mape_a))*100;
    
end