close all;
clear all;

% 
% load("G:\FiberTest\ch1_25_fA\Data.mat");
% ch1_1_fa = Data;
% load("G:\FiberTest\ch4_25_fa\Data.mat");
% ch4_1_fa = Data;

load("G:\FiberTest\ch4_10_fA\Data.mat");
ch1_1_fa = Data;
load("G:\FiberTest\ch4_15_fA\Data.mat");
ch1_15_fa = Data;
load("G:\FiberTest\ch4_20_fA\Data.mat");
ch1_20_fa = Data;
load("G:\FiberTest\ch4_25_fA\Data.mat");
ch1_25_fa = Data;

%% Checking the changes in the g2 curves over time
for  i=1:6000

    semilogx(Data_tau, squeeze(Data(i,3,:)));
    ylim([1 1.5])

    pause(0.01)
    

end

%% Analyzing g2 curves for the fiber test.

semilogx(Data_tau, mean(squeeze(ch1_1_fa(:,4,:))));
hold on; 
semilogx(Data_tau, mean(squeeze(ch1_15_fa(:,4,:))));
semilogx(Data_tau, mean(squeeze(ch1_20_fa(:,4,:))));
semilogx(Data_tau, mean(squeeze(ch1_25_fa(:,4,:))));
ylabel("g_2(\tau)");
xlabel("\tau(s)");
title("ch1 vs ch4 1cm")
legend("1 cm","1.5 cm","2 cm", "2.5 cm");



%% Averaging the g2 curves and then fitting for the response and to remove the pulsatility and noise due to it.


t_res=0.05; % seconds
t_avg=1; % window used for averaging in seconds

t_avg_pt=t_avg/t_res; % window used for averaging in points
time=t_res*(1:1:size(Data,1));
% Data_time(1,i)=i*t_res;
%     aDb1 = adb_avg

j=1;
for i=2:t_avg_pt:size(Data,1)
    Data_avg(j,1,:)=mean(Data(i-1:i+t_avg_pt-2,1,:));
    Data_avg(j,2,:)=mean(Data(i-1:i+t_avg_pt-2,2,:));
    Data_avg(j,3,:)=mean(Data(i-1:i+t_avg_pt-2,3,:));
    Data_avg(j,4,:)=mean(Data(i-1:i+t_avg_pt-2,4,:));

%     time_avg(1,j)=(Data_time(1,i+t_avg_pt-1));
    j=j+1;
end

