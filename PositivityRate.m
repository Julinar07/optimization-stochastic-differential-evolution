excel = readtable('JAPAN_data.xlsx');
N = length(table2array(excel(:,8)));

infected_case = table2array(excel(:,8));
total_positive = table2array(excel(:,3));

DATE = datetime(2020,01,16)+caldays(0:N-1);
DATE = DATE';

% Based on actual data
PosRate = (infected_case./total_positive)*100;

% Based on estimation data
segmen = 6;
r=1;
KumpulanInfect = xlsread('Infected.xlsx');
KumpulanRecover = xlsread('Recovered.xlsx');

NilaiInfected = zeros(1,N);
tes=0;

for j=1:N
    for i=1:N-segmen
        if KumpulanInfect(i,j) ~= 0
            tes = tes+KumpulanInfect(i,j);
            r=r+1;
        end
    end
    NilaiInfected(j)=tes/r;
    tes=0; r=0;
end
NilaiInfected(1) = KumpulanInfect(1,1);

NilaiRecover = zeros(1,N);
tes=0;

for j=1:N
    for i=1:N-segmen
        if KumpulanRecover(i,j) ~= 0
            tes = tes+KumpulanRecover(i,j);
            r=r+1;
        end
    end
    NilaiRecover(j)=tes/r;
    tes=0; r=0;
end
NilaiRecover(1) = KumpulanRecover(1,1);

TotPos = NilaiInfected + NilaiRecover;
Positif = (NilaiInfected./TotPos)*100;


% Plot
figure(1);
plot(DATE,PosRate,'y','linewidth',1.5);
title('Postivity Rate COVID-19 in Japan');
ylabel('Percentage');
xlabel('Date');
xlim([DATE(1) DATE(end)]);
grid on;

figure(2);
plot(DATE,Positif,'y','linewidth',1.5);
title('Postivity Rate COVID-19 in Japan');
ylabel('Percentage');
xlabel('Date');
xlim([DATE(1) DATE(end)]);
grid on;