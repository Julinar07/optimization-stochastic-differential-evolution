% Global Project Based Learning (Group 6)
% Modelling COVID-19 in Japan using SIR model (DETERMINISTIC+STOCHASTIC)
clc; clear; close all; format longg;

%% DETERMINISTIC SECTION
data = readtable('JAPAN_data.xlsx');                % Calling data into MATLAB
infect_real = table2array(data(:,8));
recover_real = table2array(data(:,4));
actual_data = [infect_real recover_real];

% Problem settings
lb = [2e-3*125000000 0 0];                          % Lower bound
ub = [4e-2*125000000 5e-8 2e-4];                    % Upper bound
T = 1500;                                           % Maximum iterations
Np = 150;                                           % Number of generate
PCr = 0.9;                                          % Crossover probability
F = 0.85;                                           % Mutation rate
d = length(lb);

% Inisialization
In = infect_real(1,1);                              % Get the initial value of infected
Re = recover_real(1,1);                             % Get the initial value of recovered
rmse = zeros(1,Np);                                 % Just make some vector of error measure
save = zeros(1,T);                                  % Just make some vector of save best mape
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
        
        for f = 1:3
            if v(f) < lb(f)
                v(f) = (1+rand())*lb(f);
            end
            if v(f) > ub(f)
               v(f) = rand()*ub(f);
            end
        end
        
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
disp('');
disp('The Parameter of SIR Model');
disp('-------------------------------------------------');
fprintf('Parameter S0 = %d\n',S0);
fprintf('Parameter r  = %d\n',rBest);
fprintf('Parameter a  = %d\n',aBest);
fprintf('RMSE         = %d\n',save(T));

N=length(infect_real);
DATE=datetime(2020,01,29)+caldays(0:N-1);

% Plot Estimation SIR Model
figure(1)
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

%% STOCHASTIC SECTION
VarMinStoch = [0.2 1200 700];                      % Lower bound sigma
VarMaxStoch = [1e4 9e4 1e3];                        % Upper bound sigma

saveS = zeros(length(infect_real),Np);
saveInf = zeros(length(infect_real),Np);
saveRec = zeros(length(infect_real),Np);
saveStoch = zeros(1,T);

% Generate parameter sigma
Sigma1 = repmat(VarMinStoch(1),Np,1) + repmat((VarMaxStoch(1)-VarMinStoch(1)),Np,1).*rand(Np,1);
Sigma2 = repmat(VarMinStoch(2),Np,1) + repmat((VarMaxStoch(2)-VarMinStoch(2)),Np,1).*rand(Np,1);
Sigma3 = repmat(VarMinStoch(3),Np,1) + repmat((VarMaxStoch(3)-VarMinStoch(3)),Np,1).*rand(Np,1);
vecX = [Sigma1 Sigma2 Sigma3];

for n=1:T
    for i = 1:Np
        % Mutation stage
        id = randi(Np,3,1);
        xp = vecX(id(1),:); xq = vecX(id(2),:); xr = vecX(id(3),:);
        v = abs(xp + F*(xq - xr));
        
        for f = 1:3
            if v(f) < VarMinStoch(f)
               v(f) = (1+rand())*VarMinStoch(f);
            end
            if v(f) > VarMaxStoch(f)
               v(f) = rand()*VarMaxStoch(f);
            end
        end
        
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
        if n == 1
            sig1 = vecX(i,1);
            sig2 = vecX(i,2);
            sig3 = vecX(i,3);
            [FitReal, ValueReal] = StochasticSIR(S0, rBest, aBest, sig1, sig2, sig3, dt, actual_data);
        else
            ValueReal = [saveS(:,i) saveInf(:,i) saveRec(:,i)];
            FitReal = rmse(i);
        end
        sig1u = u(1);
        sig2u = u(2);
        sig3u = u(3);
        [FitEsti, ValueEst] = StochasticSIR(S0, rBest, aBest, sig1u, sig2u, sig3u, dt, actual_data);
        % Selection Stage
        if FitEsti <= FitReal
            vecX(i,:) = u;
            FitReal=FitEsti;
            ValueReal = ValueEst;
        end
        saveS(:,i) = ValueReal(:,1);
        saveInf(:,i) = ValueReal(:,2);
        saveRec(:,i) = ValueReal(:,3);
        rmse(i)=FitReal;
    end
    
    sigBest = vecX(i,:);
    saveStoch(n)=min(rmse);
    fprintf('Iteration %3d completed\n',n);
end

[m,t] = min(rmse);
BestSol = [saveInf(:,t) saveRec(:,t)];
sig1=sigBest(1);
sig2=sigBest(2);
sig3=sigBest(3);

% Display the parameter
disp('');
disp('The Parameter of SIR Model');
disp('-------------------------------------------------');
fprintf('Parameter Sigma S = %d\n',sig1);
fprintf('Parameter Sigma r = %d\n',sig2);
fprintf('Parameter Sigma a = %d\n',sig3);
fprintf('RMSE Stochastic   = %d\n',saveStoch(T));

N=length(infect_real);
DATE=datetime(2020,01,16)+caldays(0:N-1);

% Plot RMSE Stochastic
figure(2)
plot(saveStoch,'-ok','MarkerFaceColor','r','MarkerSize',4);
title('Root Mean Square Error (RMSE) Stochastic');
xlabel('Iteration');
ylabel('RMSE');

% Plot Estimation SIR Model
figure(3)
plot(DATE,BestSol(:,1),'-r','linewidth',1.5);
hold on;
plot(DATE,BestSol(:,2),'-g','linewidth',1.5);
hold on;
plot(DATE,actual_data(:,1),'*r',DATE,actual_data(:,2),'*g');
xlim([DATE(1) DATE(end)]);
title('Estimation COVID-19 Japan Using SIR Stochastic');
xlabel('Date');
ylabel('Number of Cases');
legend('DE Infected','DE Recovered','Actual Infected','Actual Recovered');
grid on;