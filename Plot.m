%% Plot Infected
excel = readtable('JAPAN_data.xlsx');
confirmed = table2array(excel(:,8));

% plot(confirmed,'r','linewidth',2);

%% Plot Rt beserta CI nya
cd 'D:\GPBL SIT\Project\Final Project';
file = xlsread('RtClean.xlsx');

for i = 1:size(file,1)
    file(i,:) = smooth(file(i,:),10);
end
plot(file');

mu = mean(file);
sigma = std(file);
z = 1.96;

X1 = mu + z*sigma;
X2 = mu - z*sigma;

atas = X1';
bawah = X2';
tengah = mu';
% xData=1:1:300;
% yyaxis left;
% x=7:size(file,2)+6;
% plot(x,mu,'b','linewidth',1);
% hold on;
% X=[x,fliplr(x)];                %#create continuous x value array for plotting
% Y=[X1,fliplr(X2)];              %#create y values for out and then back
% fill( X, Y ,1,'facecolor','blue','edgecolor','none','facealpha', 0.1);
% ylim([0 7]);
% xlim([0 300]);
% yyaxis right;
% yyaxis right;
% plot(xData,confirmed,'r','linewidth',1.5);
% ylabel('Number of Cases');
% set(gca,'fontsize',10);

N = length(confirmed);
DATE = datetime(2020,01,16)+caldays(0:N-1);
DATE = DATE';

startDate = datenum('1-16-2020');
endDate = datenum('11-10-2020');
xData = linspace(startDate,endDate,300);
xData7 = xData(7:300);
x=1:1:300;

yyaxis left;
plot(xData7,tengah,'b','linewidth',1.5);
hold on;
X=[xData7,fliplr(xData7)];
Y=[atas',fliplr(bawah')];
fill(X, Y ,1,'facecolor','blue','edgecolor','none','facealpha', 0.2);
ylim([0 6.3]);
set(gca,'fontsize',15);
ylabel('Effective Reproduction Number (Rt)');
xlabel('Month');
title('Effective Reproduction Number COVID-19 in Japan');
yline(1,'g','linewidth',1.25);
yticks(1:6);
xticks(0:30:300);
ax = gca;
ax.XTick = xData;
datetick('x','mmm','keeplimits');
xlim([xData(1) xData(end)]);
yyaxis right;
plot(xData,confirmed,'r','linewidth',1.5);
ylabel('Number of Cases');
xticks(0:30:300);
ax = gca;
ax.XTick = xData;
legend('Rt','95% Confidence Interval','Rt Boundary','Infected Case','Location', 'NorthEast');
datetick('x','mmm','keeplimits');
xlim([xData(1) xData(end)]);
grid on;