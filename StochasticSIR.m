function [rmse, Y] = StochasticSIR(S0, r, a, sig1, sig2, sig3, dt, actual_data)
    % SIR Model Initialization
    S = []; I = []; R = [];
    S(1) = S0; I(1) = actual_data(1, 1);
    R(1) = actual_data(1, 2);
    
    % Brownian Motion Initialization
    S1 = []; S2 = []; S3 = [];
    S1(1) = 0; S2(1) = 0; S3(1) = 0;
%     N=length(actual_data);
    
    for n = 1:size(actual_data, 1)-1
        S(n+1) = abs(S(n) - r*S(n)*I(n));
        I(n+1) = abs(I(n) + (r*S(n)*I(n) - a*I(n)));
        R(n+1) = abs(R(n) + a*I(n));
        
        S1(n+1) = abs(S1(n) + (sig1*randn()*sqrt(dt)));
        S2(n+1) = abs(S2(n) + (sig2*randn()*sqrt(dt)));
        S3(n+1) = abs(S3(n) + (sig3*randn()*sqrt(dt)));
    end
    
    S = S + S1;
    I = I + S2;
    R = R + S3;
    
    Y = [S' I' R'];
    
    % Compute error
    rmse = mean(mean(sqrt((Y(:,2:3)-actual_data).^2)));
%     mape_r = mean(abs(((actual_data(:,1)-Y(:,2))./actual_data(:,1))));
%     mape_a = mean(abs(((actual_data(:,2)-Y(:,3))./actual_data(:,2))));
%     mape = mean(mean(mape_r+mape_a))*100;
    
end