close all;


% This is a list of files corresponding to the LBNP pressure cuff experiment
% file_list = ["20220218 - 8\","20220217\","20220217 - 6\","20220218 - 4\","20220223 - 2\","20220223 - 6\","20220304\","20220309 - 5\","20220309 - 7\"];
% file_name = ["MPCM002","MPCM003","MPCM004","MPCM005","MPCM006","MPCM007","MPCM008","MPCM009","MPCM010"];



% Read files;
dcs_dir = "C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\MPCM Study\DCS\20220223 - 2\Data.mat" ;     % Directory for DCS data
param_dir = nan;   % Directory for ECG,BP data

load(dcs_dir);
% param_data = load(param_dir);

%% Deviding the data in individual parameter
% Change it according to the available parameters in dataset

l = 60; % length of signal as ROI in seconds
s = 1; % Signal ROI starting point in seconds;

ecg_a = data(datastart(1):dataend(1));
bp_a = data(datastart(2):dataend(2));
% tcd_a = data(datastart(3):dataend(3));

ecg1 = ecg_a(1:300000);
ecg1 = normalize(ecg1);
ecg1 = lpf(ecg1,5,1000);
    
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
%% Averaging thr g2 curve
ecg_ad = circshift(ecg1,-700); % Advancing the ECG signal to match the DCS signal
% Finding the R-R peaks of ECG signal

if exist("g2")
    clear g2
end
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10; 

g2 = Data;
g2_avg = zeros(size(g2))*NaN;



[h_pks,l_pks] = findpeaks(normalize(ecg1),"MinPeakHeight",2.5);

fig1=figure('units','centimeters', 'Position',[2 2 35 13]); %18 width 15 heigh
hold on
plot(normalize(ecg1),'r')
plot(l_pks, h_pks,'*k')
hd_pks = floor(h_pks./50);
ld_pks = floor(l_pks./50);
% g2_avg = zeros(size(g2,1),size(g2,2),size(g2,3));

for i=1:size(ld_pks,2)-5
    min_length = min(diff(ld_pks(i:i+5)));
    base_sig = g2(ld_pks(i):ld_pks(i)+min_length,:,:);
    for j=1:4
        base_sig = base_sig + g2(ld_pks(i+j):ld_pks(i+j)+min_length,:,:);
%         adb_1 = hybrid_dcs(base_sig,Data_tau);
       
    end
    g2(ld_pks(i):ld_pks(i)+min_length,:,:) = base_sig./5;
end

adb_avg = hybrid_dcs(g2,Data_tau);
figure();
aDb1 = hybrid_dcs(Data,Data_tau);
