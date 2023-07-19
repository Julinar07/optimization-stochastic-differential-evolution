% Global Project Based Learning (Group 6)
% Modelling COVID-19 in Japan using SIR model (DETERMINISTIC)
clc; clear; close all; format longg;

data = readtable('JAPAN_data.xlsx');                % Calling data into MATLAB
infect_real = table2array(data(:,8));
recover_real = table2array(data(:,4));
actual_data = [infect_real recover_real];

% Problem settings
lb = [0 0 0];                                       % Lower bound
ub = [2e-2*125000000 6e-8 2e-6];                    % Upper bound
T = 15000;                                           % Maximum iterations
Np = 150;                                           % Number of generate
PCr = 0.9;                                          % Crossover probability
F = 0.85;                                           % Mutation rate
d = length(lb);

% Inisialization
In = infect_real(1,1);                              % Get the initial value of infected
Re = recover_real(1,1);                             % Get the initial value of recovered
rmse = zeros(1,Np);                                 % Just make some vector of error measure
save = zeros(1,T);                                  % Just make some vector of save best rmse
ValueReal = [];
ValueEst = [];
dt = 1/length(infect_real);

% Generate parameter
ParamS = repmat(lb(1),Np,1) + repmat((ub(1)-lb(1)),Np,1).*rand(Np,1);
Paramr = repmat(lb(2),Np,1) + repmat((ub(2)-lb(2)),Np,1).*rand(Np,1);
Sigma3 = repmat(lb(3),Np,1) + repmat((ub(3)-lb(3)),Np,1).*rand(Np,1);
vecX = [ParamS Paramr Sigma3];

for n=1:T
    for i = 1:Np
        % Mutation stage
        id = randi(Np,3,1);
        xp = vecX(id(1),:); xq = vecX(id(2),:); xr = vecX(id(3),:);
        v = abs(xp + F*(xq - xr));
        Jr = randi(d);
        r = rand();
        % Crossover stage
        for j = 1:d
            if r <= PCr || j == Jr
                u = v;
            else
                u = vecX(i,:);
            end
        end
        S = vecX(i,1);
        a = vecX(i,2);
        rPar = vecX(i,3);
        [FitReal, ValueReal] = DeterministicSIR(S, a, rPar, actual_data);
        Su = u(1);
        au = u(2);
        ru = u(3);
        [FitEsti, ValueEst] = DeterministicSIR(Su, au, ru, actual_data);
        % Selection Stage
        if FitEsti <= FitReal
            vecX(i,:) = u;
            FitReal=FitEsti;
        end
        rmse(i)=FitReal;
    end
    ParamBest = vecX(i,:);
    save(n)=min(rmse);
    fprintf('Iteration %3d completed, RMSE %5d\n',n,min(rmse));
end

S0=ParamBest(1);
rBest=ParamBest(2);
aBest=ParamBest(3);
[FitBest, y] = DeterministicSIR(S0, rBest,aBest,actual_data);

% Display the parameter
disp('The Parameter of SIR Model');
disp('-------------------------------------------------');
fprintf('Parameter S0 = %d\n',S0);
fprintf('Parameter r  = %d\n',rBest);
fprintf('Parameter a  = %d\n',aBest);
fprintf('RMSE         = %d\n',save(T));

N=length(infect_real);
DATE=datetime(2020,01,16)+caldays(0:N-1);

% Plot RMSE
figure(1)
plot(save,'-ok','MarkerFaceColor','r','MarkerSize',4);
title('Root Mean Square Error (RMSE)');
xlabel('Iteration');
ylabel('RMSE');

% Plot Estimation SIR Model
figure(2)
plot(DATE,y(:,1),'-b','linewidth',1.5);
hold on;
plot(DATE,y(:,2),'-r','linewidth',1.5);
hold on;
plot(DATE,y(:,3),'-g','linewidth',1.5);
plot(DATE,actual_data(:,1),'*r',DATE,actual_data(:,2),'*g');
xlim([DATE(1) DATE(end)]);
title('Estimation COVID-19 Japan');
xlabel('Date');
grid on;
ylabel('Population');
legend('Susceptible','DE Infected','DE Recovered','Actual Infected','Actual Recovered');
grid on;