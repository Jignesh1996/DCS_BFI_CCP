close all;

% Read files;
dcs_dir = "C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\MPCM Study\DCS\20220218 - 8\Data.mat" ;     % Directory for DCS data
param_dir = nan;   % Directory for ECG,BP data

load(dcs_dir);
% param_data = load(param_dir);

%% Deviding the data in individual parameter
% Change it according to the available parameters in dataset
ecg_a = data(datastart(1):dataend(1));
bp_a = data(datastart(2):dataend(2));
% tcd_a = data(datastart(3):dataend(3));

ecg1 = ecg_a;
ecg1 = normalize(ecg1);
    
% tcd = tcd_a(1:length(ecg1));
% % tcd = normalize(tcd);
% tcd = lpf_ffilts(tcd,30,1000);

bp_a = bp_a(1:length(ecg1));
bp_a = lpf(bp_a,3,1000);

%% Processing of the DCS signal

aDb1 = hybrid_dcs(Data,Data_tau);

%% Data plotting

% time resultion - aqusition time used to aquire data
figure();
t_res=1; % seconds
time=t_res*(1:1:size(aDb1,2));

subplot(2,2,1)

plot(time,aDb1(1,:))
title('{\itr}_{SD}=1 cm')
% set(gca,'xticklabel',{})

subplot(2,2,2)

plot(time,aDb1(2,:))
title('{\itr}_{SD}=1.5 cm')
xlabel('Time (s)')

subplot(2,2,3)

plot(time,aDb1(3,:))
title('{\itr}_{SD}=2 cm')
% set(gca,'xticklabel',{})

subplot(2,2,4)

plot(time,aDb1(4,:))
title('{\itr}_{SD}=2.5 cm')
xlabel('Time (s)')

%% Assigning the channels
dcs_1cm = aDb1(1,:).*10^9;
dcs_1lp = lpf_ffilts(dcs_1cm,10,20);
dcs_15 = aDb1(2,:).*10^9;
dcs_15lp = lpf_ffilts(dcs_15,10,20);
dcs_2 = aDb1(3,:).*10^9;
dcs_2lp = lpf_ffilts(dcs_2,10,20);
dcs_25 = aDb1(4,:).*10^9;
dcs_25lp = lpf_ffilts(dcs_25,7,20);

%% Recombine the DCS signal
adb_lp = [dcs_1lp;dcs_15lp;dcs_2lp;dcs_25lp];
breakpoints = 1:50:length(dcs_1lp);
for k=1:size(adb_lp,1)
    adb_lp(k,:) = detrend(adb_lp(k,:),1,breakpoints);
end
    
%% Detrending the signal
str_var = ["DCS 1cm","DCS 1.5cm","DCS 2cm","DSC 2.5cm"];
per_change = zeros(size(adb_lp));
for i=1:size(aDb1,1)
    s = adb_lp(i,:);
    s_str = str_var(i);
    breakpoints = 1:50:length(s);
    y = detrend(s,1,breakpoints);
%     plot(y);
    
    % Plotting the mean average smoothing curve with the envelope
    x = (1:1:length(y))*0.05;
    [envHigh, envLow] = envelope(y,20,'peak');
    envMean = (envHigh+envLow)/2;
    env_diff = envHigh - envLow;
    figure();
    plot(x,y,x,envLow,x,envHigh,x,envMean);
    [mx,my] = max(y);
    snr_val = zeros(4,5);
    snr_val(i,1) = db2mag(snr(s(1:800),20));
    snr_val(i,2) = db2mag(snr(s(800:1200),20));
    snr_val(i,3) = db2mag(snr(s(1200:4900),20));
    snr_val(i,4) = db2mag(snr(s(4900:5400),20));
    snr_val(i,5) = db2mag(snr(s(5400:6000),20));
    str = [strcat("SNR=",string(snr_val(1))),strcat("SNR=",string(snr_val(2))),strcat("SNR=",string(snr_val(3))),strcat("SNR=",string(snr_val(4))),strcat("SNR=",string(snr_val(5)))];
    
    xt = [10,40,135,240,265];
    yt = [mx,mx-2,mx-4,mx-2,mx];
    text(xt,yt,str,'Color','red','FontSize',14);
    text((max(x)/2)-10, max(envHigh)-1,s_str, 'FontSize',16)
    
    
    title("Pulsatility envelope");
    xlabel("Time (s)");
    ylabel("DCS 1cm aDb");
    
    figure();
    env_diff_smooth = smooth(env_diff,200);
    plot(x,env_diff_smooth);
    title("Pulse amplitude normalized "+s_str);
    xlabel("Time (s)");
    ylabel("Pulse Amplitude");
    % % Change in the pulsatility
    avg_sig = smooth(s,150);
    per_change(i,:) = (100*(env_diff_smooth)./mean(env_diff_smooth(1:700))) ;
    plot(x,per_change(i,:));
    title("Percent change in pulse amplitude with pressure cuff " +s_str);
    xlabel("Time (s)");
    ylabel("% change in Pulse Amplitude");
end
%% Plotting the changes in each channel in single plot
plot(x,per_change(1,:),'DisplayName','per_change_dcs_25');
hold on;
plot(x,per_change(2,:),'DisplayName','per_change_dcs_15');
plot(x,per_change(3,:),'DisplayName','per_change_dcs_1');
plot(x,per_change(4,:),'DisplayName','per_change_dcs_2');
hold off;
legend("DCS 1cm","DCS 1.5cm","DCS 2cm","DCS 2.5cm"); 
title("Percent change in pulse amplitude with pressure cuff " +s_str);
xlabel("Time (s)");
ylabel("% change in Pulse Amplitude");

%% Calculating the piecewise SNR in signal


% Steps to finish in this processing
% Plot the changes in the pulsitility of the signal and SNR 
% Plot the SNR graph with frequency specreum of SNR for each phase
