% clear all
close all
% 
% load('dcs_1_baseline_forearm.mat')
% load('ECG_DCS_forearm_exp.mat')



dcs =aDb1(2,1:600);
% uf = length(ecg1)/length(dcs_1cm); 
ecg1 = normalize(ecg1);
 % shifting only if it is a DCS signal else not so condition is check the length of the signal
ecg1 = circshift(ecg1,00);
ecg_da = ecg1(1:size(dcs,2)*50);


break_pt = 1:100:length(dcs);
dcs_d = detrend(normalize(dcs),1,break_pt); 
% [dcs_1,filter_obj] = bandpass(dcs,[0.5 6],20,'ImpulseResponse','fir');
dcs_1 = dcs_d;
plot(dcs_1)
time_DCS=0.05*(1:1:size(dcs_1,2));
time_ECG=0.001*(1:1:size(ecg_da,2));

% time_DCS=0.001*(1:1:size(dcs_1,2));
% time_ECG=0.001*(1:1:size(ecg_da,2));

figure()
hold on
plot(time_DCS,dcs_1/max(dcs_1),'r')
plot(time_ECG,ecg_da/max(ecg_da),'b')


% Find peaks

[pks_DCS,locs_DCS]=findpeaks(dcs_1/max(dcs_1),'MinPeakHeight',0.65)
[pks_ECG,locs_ECG]=findpeaks(ecg_da/max(ecg_da),'MinPeakHeight',0.5)

locs_DCS_time=locs_DCS*0.05;
locs_ECG_time=locs_ECG*0.001;

time_shift1=locs_ECG_time(1)-locs_DCS_time(1) %% in s
% Upsampling the DCS singal
x = 1:1:length(dcs_1);
uf = 50;   % Upsampling factor
xq = (1:(1/uf):length(dcs_1)+((uf-1)/uf)); 
dcs_1_interp = interp1(x, dcs_1,xq,'v5cubic');
% dcs_1_interp = dcs_1;

% filtering
% 
windowsize=10; % how many points you want to use (it will depend on your resolution, we were using 10 so it was 3 s window)
wages=ones(1,windowsize)/windowsize;
% wages = window_1;
dcs_1_smooth=dcs_1_interp; % Y is your time course you want to filter, Y_smoth is filtered data
ecg1_smooth=filtfilt(wages,1,ecg_da);
break_pt = 1:2000:size(dcs_1_smooth,2);
[pks_ECG_smooth,locs_ECG_smooth]=findpeaks(ecg1_smooth./max(ecg1_smooth),'MinPeakHeight',0.65,'MinPeakDistance',600);
[pks_DCS_smooth,locs_DCS_smooth]=findpeaks(detrend(dcs_1_smooth,1,break_pt,"omitnan"),'MinPeakHeight',0.25,'MinPeakDistance',600);




%

% [pks_DCS_interp,locs_DCS_interp]=findpeaks(dcs_1_interp/max(dcs_1_interp),'MinPeakHeight',0.8)
% locs_DCS_time_interp=locs_DCS_interp*0.001;
% time_shift2=locs_ECG_time(1)-locs_DCS_time_interp(1) %% in s
% 
% time_shift2_pt=floor(time_shift2/0.001);
% 
% fig1=figure('units','centimeters', 'Position',[2 2 35 13]) %18 width 15 heigh
% subplot(1,2,1)
% hold on
% plot(time_ECG,dcs_1_interp/max(dcs_1_interp),'r');
% plot(time_ECG,ecg_da/max(ecg_da),'b')
% 
% subplot(1,2,2)
% hold on
% temp=dcs_1_interp/max(dcs_1_interp);
% temp2=circshift(temp(1,:),5);
% plot(time_ECG,circshift(temp(1,:),time_shift2_pt),'r');
% plot(time_ECG,ecg_da/max(ecg_da),'b')

%

fig1=figure('units','centimeters', 'Position',[2 2 35 13]) %18 width 15 heigh
hold on
plot(detrend(dcs_1_smooth,1,break_pt,"omitnan"),'r')
plot(locs_DCS_smooth, pks_DCS_smooth,'*k')

plot(ecg1_smooth/max(ecg1_smooth),'b')
plot(locs_ECG_smooth, pks_ECG_smooth,'*k')

% if locs_ECG_smooth(1)>locs_DCS_smooth(1)
%     Difference=locs_ECG_smooth(1)-locs_DCS_smooth(1);
%     shift = floor(1.75*Difference);
% else
%     shift =0;
% end
% Extracting data

Extract=ones(size(pks_ECG_smooth,2)-3,min(diff(locs_ECG_smooth)));
Extract=Extract*NaN;

dcs_1_smooth2=circshift(dcs_1_smooth,-850);
figure()
plot(dcs_1_smooth2);
hold on;
plot(ecg1_smooth);
hold off;
% dcs_1_smooth2=tcd;


for i=2:size(pks_ECG_smooth,2)-2
    locs_ECG_smooth(i)
    locs_ECG_smooth(i)+min(diff(locs_ECG_smooth));
    signal = dcs_1_smooth2(locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1);
%     base_sig = Data(locs_ECG_smooth(1):locs_ECG_smooth(1)+min(diff(locs_ECG_smooth))-1,:,:);
%     sig = base_sig+ sig;
    
  
    
    Extract(i-1,:)=(dcs_1_smooth2(1,locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1));
    

%     Extract(i,:)=signal

end

% % Testing with adding the nan values at the end of tail
% 
% Extract=ones(size(pks_ECG_smooth,2)-3,(max(diff(locs_ECG_smooth))));
% Extract=Extract*NaN;
% 
% dcs_1_smooth2=circshift(dcs_1_smooth,900);
% 
% for i=2:size(pks_ECG_smooth,2)-2
%     locs_ECG_smooth(i)
%     locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))
%     signal = (dcs_1_smooth2(1,locs_ECG_smooth(i):locs_ECG_smooth(i+1)));
% %     signal(size(signal,2):size(Extract))
% 
%     Extract(i-1,1:size(signal,2))=signal;
%     
% 
% %     Extract(i,:)=signal
% 
% end

% Plotting the shifted DCS and ECG signal
% x = (1:1:length(dcs_1_smooth2))/1000;
% plot(x,normalize(dcs_1_smooth2));
% hold on;
% plot(x,normalize(ecg1_smooth));
% xlabel("Time (s)");
% legend("DCS 1cm","ECG")

%    
x = (1:1:length(Extract'))/1000;
figure();
plot(x,(Extract'))
xlabel("Time(s)")
ylabel("aDb value")
title("Marker=ECG R peak, DCS 2.5cm")

% Plot ensemble average

ttle = 'DCS 2.5cm Ensemble Avg';
avg_dcs_1 = ens_avg(Extract,ttle);

%% 
% plot(x,avg_dcs,'b');
% hold on; 
% plot(x,avg_tcd_shift,'r');
% plot(x,avg_bp_shift,'k');
% legend('Aerage DCS 1cm','Average TCD','Average BP')
% xlabel("Time (s)");
% title("Comparision of average plots of TCD, BP, and DCS 1 cm plot")

%% Finding CCP for individual cycles
ccp = zeros(1,length(bp_stack(:,1)));
for i=1:length(bp_stack(:,1))
    p = polyfit(bp_stack(i,:)',dcs_stack(i,:)',1);
    x = -20:0.005:130;
    f = polyval(p,x);
    fprintf("%d",i);
    ccp(i) = mean(x(round(f,2)==0));
end

%% Taking average of the cosecutive 5 cycles to reduce the effect of noisy signal

step_size = 5; 
Extract=zeros(size(pks_ECG_smooth,2)-3,ceil(size(pks_ECG_smooth,2)/step_size),min(diff(locs_ECG_smooth)));
Extract=Extract*NaN;

dcs_1_smooth2=circshift(dcs_1_smooth,775)

for i=1:step_size:size(pks_ECG_smooth,2)-2
    sig = dcs_1_smooth2(1,locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1);
    count = 0;
    for j=1:1:step_size-1
        if i+j>size(pks_ECG_smooth,2)-2
            break;
        end
        locs_ECG_smooth(i+j);
        locs_ECG_smooth(i+j) + min(diff(locs_ECG_smooth));
        signal = dcs_1_smooth2(1,locs_ECG_smooth(i+j):locs_ECG_smooth(i+j)+min(diff(locs_ECG_smooth))-1);
        sig = sig+signal;
        count = count+1;
    end
    sig = sig./count;
      

    Extract(sig,floor(i/step_size)+1,:)=sig;
end

x = (1:1:length(Extract'))/1000;
plot(x,(Extract'))
xlabel("Time(s)")
ylabel("aDb value")
title("Marker=ECG R peak, DCS 2.5cm")

%%
clc;
avg = ccp_measure(ecg_da,dcs_1lp,break_pt)

