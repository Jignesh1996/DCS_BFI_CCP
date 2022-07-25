function [ens_avg_sig,pind,ensemble_curve] = ensemble_avg(ecg, signal, shift,mode)
ecg1 = ecg;
dcs = signal;
shift = shift;
mode = mode;
ecg1 = normalize(ecg1);
ecg1 = circshift(ecg1,0);
ecg_da = ecg1;


break_pt = 1:50:length(dcs);
dcs_d = detrend(normalize(dcs),1,break_pt); 
dcs_1 = dcs;
plot(dcs_1)
time_DCS=0.05*(1:1:size(dcs_1,2));
time_ECG=0.001*(1:1:size(ecg_da,2));

figure()
hold on
plot(time_DCS,dcs_1/max(dcs_1),'r')
plot(time_ECG,ecg_da/max(ecg_da),'b')


% Find peaks

[pks_DCS,locs_DCS]=findpeaks(dcs_1/max(dcs_1),'MinPeakHeight',0.65)
[pks_ECG,locs_ECG]=findpeaks(ecg_da/max(ecg_da),'MinPeakHeight',0.8)

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
% ecg1_smooth=filtfilt(wages,1,ecg_da);
ecg1_smooth=ecg_da;
break_pt = 1:2000:size(dcs_1_smooth,2);
[pks_ECG_smooth,locs_ECG_smooth]=findpeaks(ecg1_smooth./max(ecg1_smooth),'MinPeakHeight',0.65,'MinPeakDistance',600);
[pks_DCS_smooth,locs_DCS_smooth]=findpeaks(detrend(dcs_1_smooth,1,break_pt,"omitnan"),'MinPeakHeight',0.25,'MinPeakDistance',600);



% fig1=figure('units','centimeters', 'Position',[2 2 35 13]) %18 width 15 heigh
% hold on
% plot(detrend(dcs_1_smooth,1,break_pt,"omitnan"),'r')
% plot(locs_DCS_smooth, pks_DCS_smooth,'*k')
% 
% plot(ecg1_smooth/max(ecg1_smooth),'b')
% plot(locs_ECG_smooth, pks_ECG_smooth,'*k')




dcs_1_smooth2=circshift(dcs_1_smooth,shift);
figure()
plot(dcs_1_smooth2);
hold on;
plot(ecg1_smooth);
hold off;

% 
if mode==0
    Extract=ones(size(pks_ECG_smooth,2)-3,min(diff(locs_ECG_smooth)));
    Extract=Extract*NaN;
    for i=2:size(pks_ECG_smooth,2)-2
        locs_ECG_smooth(i)
        locs_ECG_smooth(i)+min(diff(locs_ECG_smooth));
        signal = dcs_1_smooth2(locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1);
            
        Extract(i-1,:)=(dcs_1_smooth2(1,locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1));
        d = Extract(i-1,:);
        PI(i-1) = (max(d,"omitnan")-min(d(1:500),"omitnan"))/mean(d,"omitnan");
    
    %     Extract(i,:)=signal
    
    end
end

% Testing with adding the nan values at the end of tail


% 
% dcs_1_smooth2=circshift(dcs_1_smooth,900);
% 
if mode==1
    Extract=ones(size(pks_ECG_smooth,2)-3,(max(diff(locs_ECG_smooth))));
    Extract=Extract*NaN;
    for i=2:size(pks_ECG_smooth,2)-2
        locs_ECG_smooth(i)
        locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))
        signal = (dcs_1_smooth2(1,locs_ECG_smooth(i):locs_ECG_smooth(i+1)));

        [h_peaks,l_peaks] = findpeaks(signal,'MinPeakHeight',1.5);
        disp(l_peaks)
        disp(length(l_peaks));
    %     signal(size(signal,2):size(Extract))

        if length(l_peaks)<2

            if length(signal)>=1000
                signal = signal(1:1000);
            end
            Extract(i-1,1:size(signal,2))=signal;
            d = Extract(i-1,:);
            PI(i-1) = (max(d(1:700),[],'omitnan')-min(d,[],'omitnan'))./mean(d,'omitnan');

        else
            Extract(i-1,1:size(signal,2))=nan;
        end
    
    %     Extract(i,:)=signal
    
    end
end
pind = mean(PI,"omitnan");
%    
x = (1:1:length(Extract'))/1000;
figure();
plot(x,(Extract'))
xlabel("Time(s)")
ylabel("aDb value")
title("Marker=ECG R peak, DCS 2.5cm")

% The issue with the ensemble averaging with NaNs is that the extra signal
% at the end iduces discontinuities and the signal goes down till the
% maximum time of the cycle. To compensate for that, I'm cutting the signal
% till the minimum point. 

% Plot ensemble average

ttle = 'DCS 2.5cm Ensemble Avg';
[ens_avg_sig,ensemble_curve] = ens_avg(Extract,ttle);

end