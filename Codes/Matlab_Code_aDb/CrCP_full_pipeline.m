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
    

