close all;


% This is a list of files corresponding to the LBNP pressure cuff experiment
% file_list = ["20220218 - 8\","20220217\","20220217 - 6\","20220218 - 4\","20220223 - 2\","20220223 - 6\","20220304\","20220309 - 5\","20220309 - 7\"];
% file_name = ["MPCM002","MPCM003","MPCM004","MPCM005","MPCM006","MPCM007","MPCM008","MPCM009","MPCM010"];



% Read files;
dcs_dir = "C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\MPCM Study\DCS\20220217\Data.mat" ;     % Directory for DCS data
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

ecg1 = ecg_a(1:120000);
ecg1 = normalize(ecg1);
ecg1 = lpf(ecg1,5,1000);
    
% tcd = tcd_a(1:length(ecg1));
% % tcd = normalize(tcd);
% tcd = lpf_ffilts(tcd,30,1000);

bp_a = bp_a(1:length(ecg1));
bp_a = lpf(bp_a,3,1000);
%% Processing of the DCS signal

aDb1 = standalone_dcs(Data,Data_tau);

%% Data plotting

% time resultion - aqusition time used to aquire data
figure();
t_res=0.05; % seconds
time=t_res*(1:1:size(aDb1,2));

subplot(2,1,1)

plot(time,aDb1(1,:))
title('{\itr}_{SD}=1 cm')
% set(gca,'xticklabel',{})

subplot(2,1,2)

plot(time,aDb1(2,:))
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


%% Calculating BFI
aDb1 = standalone_dcs(Data(9600:14400,:,:),Data_tau);

%% Finding the shift

close all;
ecg_ad = circshift(ecg1,-400);
x_e = (1:1:length(ecg_ad))/1000;
plot(x_e,normalize(ecg_ad),'b');
hold on;
x_d = (1:1:size(aDb1,2))/20;
plot(x_d,normalize(aDb1(1,:)),'r');
%% Averaging thr g2 curve
close all;
shift = 1:45:800;
SNR = zeros(1,length(shift));
bp_shift = circshift(bp_a,0);
ecg_ad = circshift(ecg1,0); % Advancing the ECG signal to match the DCS signal
%% Finding the R-R peaks of ECG signal

if exist("g2")
    clear g2;
end
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10; 

% g2 = Data;
% g2_avg = zeros(size(g2))*NaN;

% Defining the protocol ltime line


% ecg_1 = ecg_ad(1:300000);   % 5 min baseline data
% Data_1 = Data(1:6000,:,:);
% bp_1 = bp_shift(1:300000);
% ecg_2 = ecg_ad(300001:480000); % 3 min TNQT @80 mmHg
% Data_2 = Data(6001:9600,:,:);
% bp_2 = bp_shift(300001:480000);
% ecg_3 = ecg_ad(480001:720000); % 4 min TNQT @160 mmHg
% Data_3 = Data(9601:14400,:,:);
% bp_3 = bp_shift(480001:720000);
% ecg_4 = ecg_ad(720001:780000); % 1 min baseline
% Data_4 = Data(14401:15600,:,:);
% bp_4 = bp_shift(720001:780000);
% 
g2 = Data(9600:14400,:,:);
ecg = ecg_ad(480000:720000);

% start_time =  421;   %time in seconds
% stop_time =  540;    %time in seconds
% 
% g2 = Data_all(start_time*20:stop_time*20,:,:);
% ecg = ecg_ad(start_time*1000:stop_time*1000);

% D = Data(9600:14400,:,:);
% g2(:,1,:)=squeeze(D(:,1,:)-1); %g2-1 curve generation
% g2_2_temp=squeeze(D(:,2,:)-1); %g2-1 curve generation
% g2_3_temp=squeeze(D(:,3,:)-1); %g2-1 curve generation
% g2_4_temp=squeeze(D(:,4,:)-1); %g2-1 curve generation
% 
% average g2 curve for large source detector separation
% for i=1:size(g2,1)
%     g2(i,2,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
% end

[h_pks,l_pks] = findpeaks(normalize(ecg),"MinPeakHeight",2.5);

fig1=figure('units','centimeters', 'Position',[2 2 35 13]); %18 width 15 heigh
hold on
plot(normalize(ecg),'r')
plot(l_pks, h_pks,'*k')
hd_pks = floor(h_pks./50);
ld_pks = floor(l_pks./50);
% g2_avg = zeros(size(g2,1),size(g2,2),size(g2,3));
avg_window_width = 50;
for i=1:size(ld_pks,2)-1
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


% for i=1:size(g2,1)
%     g2(i,2,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
% end



% adb_avg = standalone_tr_dcs(g2,Data_tau);
adb_avg = standalone_dcs(g2,Data_tau);
% SNR(1,k) = snr(adb_avg(2,:),20);
figure
snr(adb_avg(2,:),20)

% adb_avg = adb_avg; % This is cut out specific portion of the signal to plot
% figure();
% adb = standalone_tr_dcs(Data_all(start_time*20:stop_time*20,:,:),Data_tau);
dcs_1lp = lpf(adb_avg(1,:),4,20);
dcs_25lp = lpf(adb_avg(2,:),4,20);
% dcs_25lp_tr = lpf(adb_avg(3,:),7,20);

% aDb1 = aDb1;
% figure();
% snr(adb_avg(2,:),20);
% figure();
plot(adb_avg(2,:),'b'); hold on; plot(aDb1(2,:),'r--');
% title("Comparision of g2 averaging for cuff data MPCM004 width=50 cycles")
% legend("g2 Averaged signal","Raw signal")
% xlabel("samples (Time = samples/20)");
% ylabel("aDb")

%% Plotting the frequency spectrum of the averaged and raw signal
Fs = 20;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
signal = adb_avg(2,:);
L = length(signal);     
% Length of signal
t = (0:L-1)*T;  

Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure();
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('FFT of raw DCS 2.5cm TNQT @ 80mmHg signal')
xlabel('f (Hz)')
ylabel('|P1(f)|')


%% Compare the 2.5cm and 1.5cm signals

t = (1:1:length(adb_avg(2,:)))/20;
plot(t,normalize(adb_avg(1,:)),'b'); hold on; plot(t,normalize(adb_avg(2,:)));
title("Comparison Between 1.5cm and 2.5cm DCS signal @TNQT 80mmHg")
ylabel("Normalized Unit");
xlabel("Time(s)");
legend("r_s_d = 1.5cm","r_s_d = 2.5cm");

%% Compare the 2.5cm and 1.5cm signals after filtering

dcs_1lp = lpf(adb_avg(1,:),7,20);
dcs_3lp = lpf(adb_avg(2,:),7,20);

t = (1:1:length(dcs_1lp))/20;
plot(t,normalize(dcs_1lp),'b'); hold on; plot(t,normalize(dcs_3lp));
title("Comparison: Filtered & g2 avg 1.5cm & 2.5cm DCS signal @TNQT 80mmHg")
ylabel("Normalized Unit");
xlabel("Time(s)");
legend("r_s_d = 1.5cm","r_s_d = 2.5cm");

%% Compare the 2.5cm signal with and without g2 averaging

t = (1:1:length(adb_avg(1,:)))/20;
plot(t,normalize(adb_avg(2,:)),'b'); hold on; plot(t,normalize(aDb1(2,:)));
title("Comparing 2.5cm signal with and without g2 averaging @TNQT 80mmHg")
ylabel("Normalized Unit");
xlabel("Time(s)");
legend("g2 averaged","original");

%% making an animation to show the g2 temporal averaging
hfig = figure();

subplot(2,4,[1 2 3 4])
span = 1:100;
t = (1:1:length(aDb1(1,span)))/20;
plot(normalize(aDb1(1,span)))
hold on;
plot((1:1:length(ecg_ad(locs_ECG(1):locs_ECG(6))))/50,normalize(ecg_ad(locs_ECG(1):locs_ECG(6))))

x_pos=[floor(locs_ECG(2)/50),floor((locs_ECG(3))/50)]; %task strat time in minutes
% txt = ["25%", "50%","150%","BSL"];
txt = ["80 mmHg", "160mmHg"];

rectangle('Position',[t(x_pos(2)),0.1*10^-9,length(locs_ECG(2):locs_ECG(3))/50,max(aDb1(1,:))],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none',...
    'LineWidth',3)


% 
% % subplot(2,1,2)
% subplot(2,4,5)
% semilogx(Data_tau,squeeze(Data(8,1,:)))
% subplot(2,4,5)
% semilogx(Data_tau,squeeze(Data(8,1,:)))
% subplot(2,4,5)
% semilogx(Data_tau,squeeze(Data(8,1,:)))
% subplot(2,4,5)
% semilogx(Data_tau,squeeze(Data(8,1,:)))
