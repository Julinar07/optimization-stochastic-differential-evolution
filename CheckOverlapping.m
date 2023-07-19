% Global Project Based Learning (GPBL) Group 6
% Check Estimate SIR Model (only 7 Days Data)
clc; clear; close all; format longg;

excel = readtable('JAPAN_data.xlsx');
number = length(table2array(excel(:,8)));
segmen = 6;

Simpan.ParameterDeterm = [];
Simpan.ParameterStok = [];
Simpan.RMSE = [];
Simpan.Time = [];
KumpulanInfect = zeros(number-segmen, number);
KumpulanRecover = zeros(number-segmen, number);

q=294;                                              % Change this value into days you want to check

%% DETERMINISTIC SECTION
data = readtable('JAPAN_data.xlsx');                % Calling data into MATLAB
infect_real = table2array(data(q:q+segmen,8));
recover_real = table2array(data(q:q+segmen,4));
actual_data = [infect_real recover_real];

% Problem settings
S0 = 125000000;
lb = 0.8*[186129.951489941 3.33324204888468e-07 0.036669480506023];
% lb = 0.8*[1e-3*S0 1e-09 1e-05];
ub = [1e-2*S0 0.99 0.99];                           % Upper bound
T = 300;                                            % Maximum iterations
Np = 250;                                           % Number of generate
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
        sig1 = rand();
        % Crossover stage
        for j = 1:d
            if sig1 <= PCr || j == Jr
                u = v;
            else
                u = vecX(i,:);
            end
        end
        S = vecX(i,1);
        a = vecX(i,2);
        r = vecX(i,3);
        [FitReal, ValueReal] = DeterministicSIR(S, a, r, actual_data);
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
%         fprintf('Iteration %3d completed\n',n);
end

S0=ParamBest(1);
rBest=ParamBest(2);
aBest=ParamBest(3);
%     [FitBest, y] = ModelDetermSIR(S0, rBest,aBest,actual_data);
ParamDeterm = [S0 rBest aBest];


%% STOCHASTIC SECTION
VarMinStoch = [3 0.1 5000];                      % Lower bound sigma
VarMaxStoch = [1e4 9e3 16e3];                    % Upper bound sigma

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
        v = xp + F*(xq - xr);

        for f = 1:3
            if v(f) < VarMinStoch(f)
                v(f) = (1+rand())*VarMinStoch(f);
            end
            if v(f) > VarMaxStoch(f)
               v(f) = rand()*VarMaxStoch(f);
            end
        end

        Jr = randi(d);
        sig1 = rand();
        % Crossover stage
        for j = 1:d
            if sig1 <= PCr || j == Jr
                u = v;
            else
                u = vecX(i,:);
            end
        end
        sig1 = vecX(i,1);
        sig2 = vecX(i,2);
        sig3 = vecX(i,3);
        [FitReal, ValueReal] = StochasticSIR(S0, rBest, aBest, sig1, sig2, sig3, dt, actual_data);
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

    Sigbest = vecX(i,:);
    saveStoch(n)=min(rmse);
%         fprintf('Iteration %3d completed\n',n);
end

[RMSE,t] = min(rmse);
sig1=Sigbest(1);
sig2=Sigbest(2);
sig3=Sigbest(3);
ParamStoch = [sig1 sig2 sig3];

KumpulanInfect(q,q:q+segmen) = saveInf(:,t);
KumpulanRecover(q,q:q+segmen) = saveRec(:,t);

Simpan(q).ParameterDeterm = array2table(ParamDeterm);
Simpan(q).ParameterStok = array2table(ParamStoch);
Simpan(q).RMSE = RMSE;
St = S0 - saveInf(segmen+1) - saveRec(segmen+1);
Np = saveS(1) + saveInf(1) + saveRec(1);

R0 = S0*(rBest/aBest);
Rt = (St/Np)*R0;

disp(['Segmen ' num2str(q)]);
disp(['Error: ' num2str(RMSE)]);
disp(['Rt: ' num2str(Rt)]);
disp('------------------------------------------------------------------------');

NilaiPakai = zeros(1,number);
tes=0;

for j=1:number
for i=1:number-segmen
    if KumpulanInfect(i,j) ~= 0
        tes = tes+KumpulanInfect(i,j);
        r=r+1;
    end
end
NilaiPakai(j)=tes/r;
tes=0; r=0;
end

excel = readtable('JAPAN_data.xlsx');
Penyakit = table2array(excel(:,8));

N=number;
[peak_infected, B] = max(NilaiPakai(1:N));
DATE = datetime(2020,01,29)+caldays(0:N-1);
DATE = DATE';

figure(1)
bar(DATE,Penyakit','y');
hold on;
plot(DATE,NilaiPakai,'-r','linewidth',2);
xlim([DATE(1) DATE(end)]);
legend('Actual Infected','DE Infected','Location', 'NorthWest');
xlabel('Date');
ylabel('Number of Cases');
title('Estimation Infected Case Japan');
grid on;

NilaiPakai = zeros(1,number);
tes=0;

for j=1:number
    for i=1:number-segmen
        if KumpulanRecover(i,j) ~= 0
            tes = tes+KumpulanRecover(i,j);
            r=r+1;
        end
    end
    NilaiPakai(j)=tes/r;
    tes=0; r=0;
end

excel = readtable('JAPAN_data.xlsx');
Sehat = table2array(excel(:,4));

N=number;
[peak_infected, B] = max(NilaiPakai(1:N));
DATE = datetime(2020,01,16)+caldays(0:N-1);
DATE = DATE';

figure(2)
bar(DATE,Sehat','g');
hold on;
plot(DATE,NilaiPakai,'-b','linewidth',1.5);
xlim([DATE(1) DATE(end)]);
legend('Actual Recovered','DE Recovered','Location', 'NorthWest');
xlabel('Date');
ylabel('Number of Cases');
title('Estimation Recovered Case Japan');
grid on;
set(gca, 'fontsize', 15);
set(gca, 'fontsize', 15);