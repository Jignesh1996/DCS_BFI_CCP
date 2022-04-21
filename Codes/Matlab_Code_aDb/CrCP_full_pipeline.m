close all;

% Read files;
dcs_dir = "C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\MPCM Study\DCS\20220218 - 8\Data.mat" ;     % Directory for DCS data
param_dir = nan;   % Directory for ECG,BP data

load(dcs_dir);
% param_data = load(param_dir);

%% Deviding the data in individual parameter
% Change it according to the available parameters in dataset
l = 60; % length of signal as ROI in seconds
s = 1; % Signal ROI starting point in seconds;

ecg_a = data(datastart(1):dataend(1));
bp_a = data(datastart(2):dataend(2));
tcd_a = data(datastart(3):dataend(3));

ecg1 = ecg_a(1:300000);
ecg1 = normalize(ecg1);
ecg1 = lpf(ecg1,5,1000);
    
tcd = tcd_a(1:length(ecg1));
% tcd = normalize(tcd);
tcd = lpf_ffilts(tcd,12,1000);

bp = bp_a(1:length(ecg1));
bp = lpf(bp,3,1000);

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


%%
close all;

ecg_ad = circshift(ecg1,-950); % Advancing the ECG signal to match the DCS signal
% Finding the R-R peaks of ECG signal

if exist("g2")
    clear g2;
end
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10; 

g2 = Data;
g2_avg = zeros(size(g2))*NaN;



[h_pks,l_pks] = findpeaks(normalize(ecg_ad),"MinPeakHeight",2.5);

fig1=figure('units','centimeters', 'Position',[2 2 35 13]); %18 width 15 heigh
hold on
plot(normalize(ecg_ad),'r')
plot(l_pks, h_pks,'*k')
hd_pks = floor(h_pks./50);
ld_pks = floor(l_pks./50);
% g2_avg = zeros(size(g2,1),size(g2,2),size(g2,3));
avg_window_width = 50;
for i=1:size(ld_pks,2)
    if size(ld_pks,2)-i < avg_window_width
        avg_window_width = avg_window_width-1
    end
    min_length = min(diff(ld_pks(i:i+avg_window_width)));
    base_sig = g2(ld_pks(i):ld_pks(i)+min_length,:,:);
   
    for j=1:avg_window_width-1
        base_sig = base_sig + g2(ld_pks(i+j):ld_pks(i+j)+min_length,:,:);
%         adb_1 = hybrid_dcs(base_sig,Data_tau);
       
    end
    g2(ld_pks(i):ld_pks(i)+min_length,:,:) = base_sig./avg_window_width;
end

adb_avg = hybrid_dcs(g2,Data_tau);
adb_avg = adb_avg; % This is cut out specific portion of the signal to plot
figure();
aDb1 = hybrid_dcs(Data,Data_tau);
aDb1 = aDb1;
figure();
snr(adb_avg(4,1200:5000),20);
figure();
plot(adb_avg(4,:),'b'); hold on; plot(aDb1(4,:),'r');

%% Filtering the signal

dcs_1cm = adb_avg(1,:).*10^9;
dcs_1lp = lpf_ffilts(dcs_1cm,10,20);
dcs_15 = adb_avg(2,:).*10^9;
dcs_15lp = lpf_ffilts(dcs_15,10,20);
dcs_2 = adb_avg(3,:).*10^9;
dcs_2lp = lpf_ffilts(dcs_2,10,20);
dcs_25 = adb_avg(4,:).*10^9;
dcs_25lp = lpf_ffilts(dcs_25,10,20);

%% Recombine the DCS signal
adb_lp = [dcs_1lp;dcs_15lp;dcs_2lp;dcs_25lp]; %Recombining the filtered signal
breakpoints = 1:50:length(dcs_1lp);
for k=1:size(adb_lp,1)
    adb_lp(k,:) = detrend(adb_lp(k,:),1,breakpoints);
end

%%
% Write a code for checking the quality of the signal

%% Temporal alignment the TCD, ABP, and DCS signals
break_pt = 1:100:size(tcd,2);
bp_shift = circshift(bp,0)
[a,l_bp] = findpeaks(normalize(bp_shift),"MinPeakHeight",1.5,'MinPeakDistance',500) ;
[a,l_tcd] = findpeaks(normalize(detrend(tcd,1,break_pt)),"MinPeakHeight",1);
shift = l_bp(1)-l_tcd(1);
tcd_shift = circshift(tcd,-shift);

%% CrCP calculation 

% Using full cycle method
close all;
l = 2400;
avg_pt = 60;
ccp = ccp_measure(ecg1(1:l*50),tcd_shift(1:l*50),bp_shift(1:l*50),avg_pt);
CrCP = zeros(5,size(ccp,2));
CrCP(1,:) = ccp;
ccp_tcd = mean(CrCP(1,1:end-1));
CrCP(2,:) = ccp_measure(ecg1(1:l*50),dcs_1lp(1:l),bp_shift(1:l*50),avg_pt);
ccp_dcs1 = mean(CrCP(2,1:end-1));
CrCP(3,:) = ccp_measure(ecg1(1:l*50),dcs_15lp(1:l),bp_shift(1:l*50),avg_pt);
ccp_dcs15 = mean(CrCP(3,1:end-1));
CrCP(4,:) = ccp_measure(ecg1(1:l*50),dcs_2lp(1:l),bp_shift(1:l*50),avg_pt);
ccp_dcs2 = mean(CrCP(4,1:end-1));
CrCP(5,:) = ccp_measure(ecg1(1:l*50),dcs_25lp(1:l),bp_shift(1:l*50),avg_pt);
ccp_dcs25 = mean(CrCP(4,1:end-1));

% close all
figure();
plot(CrCP(1,1:end-1));
hold on;
plot(CrCP(2,1:end-1));
plot(CrCP(3,1:end-1));
plot(CrCP(4,1:end-1));
plot(CrCP(5,1:end-1));
ylabel("CrCP (mmHg)");
title("CrCP using Tail, TCD averaged over 25 cycles")
legend("TCD","DCS 1cm","DCS 1.5cm","DCS 2cm","DCS 2.5cm");


%% Using tail of the cycle
ccp = ccp_measure_tail(ecg1(1:l*50),tcd_shift(1:l*50),bp_shift(1:l*50),avg_pt);
CrCP_t = zeros(5,size(ccp,2));
CrCP_t(1,:) = ccp;
ccp_tcd_t = mean(CrCP_t(1,1:end-1));
CrCP_t(2,:) = ccp_measure_tail(ecg1(1:l*50),dcs_1lp(1:l),bp_shift(1:l*50),avg_pt);
ccp_dcs1_t = mean(CrCP_t(2,1:end-1));
CrCP_t(3,:) = ccp_measure_tail(ecg1(1:l*50),dcs_15lp(1:l),bp_shift(1:l*50),avg_pt);
ccp_dcs15_t = mean(CrCP_t(3,1:end-1));
CrCP_t(4,:) = ccp_measure_tail(ecg1(1:l*50),dcs_2lp(1:l),bp_shift(1:l*50),avg_pt);
ccp_dcs2_t = mean(CrCP_t(4,1:end-1));
CrCP_t(5,:) = ccp_measure_tail(ecg1(1:l*50),dcs_25lp(1:l),bp_shift(1:l*50),avg_pt);
ccp_dcs25_t = mean(CrCP_t(4,1:end-1));

close all;

figure();
plot(CrCP_t(1,1:end-1));
hold on;
plot(CrCP_t(2,1:end-1));
plot(CrCP_t(3,1:end-1));
plot(CrCP_t(4,1:end-1));
plot(CrCP_t(5,1:end-1));
ylabel("CrCP_t (mmHg)");
title("CrCP_t using Tail, TCD averaged over 25 cycles")
legend("TCD","DCS 1cm","DCS 1.5cm","DCS 2cm","DCS 2.5cm");

%% Using the Harmonics (Frequency) method
close all;
CrCP_harmonics = zeros(4,2)*NaN;
CrCP_harmonics(1,:) = crcp_freq_method(ecg1,tcd,bp_a,dcs_1lp);
CrCP_harmonics(2,:) = crcp_freq_method(ecg1,tcd,bp_a,dcs_15lp);
CrCP_harmonics(3,:) = crcp_freq_method(ecg1,tcd,bp_a,dcs_2lp);
CrCP_harmonics(4,:) = crcp_freq_method(ecg1,tcd,bp_a,dcs_25lp);

%% Statistical analysis of the Critical closing pressure


