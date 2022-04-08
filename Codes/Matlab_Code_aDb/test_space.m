% Checking the upsampling factor for the interp and interp1 functions, if
% it is linear or not
y = interp(dcs_1,50);
x = 1:length(y);
x_d = (1:length(dcs_1));
[pks_u,locs_u] = findpeaks(y, 'MinPeakHeight', 0.5,'MinPeakDist',600,'MinPeakProminence',0.1);  %Determine peaks and Indices
figure()
[pks_d,locs_d] = findpeaks(dcs_1, 'MinPeakHeight', 0.5,'MinPeakProminence',10);
[pks_u1,locs_u1] = findpeaks(dcs_1up, 'MinPeakHeight', 0.5,'MinPeakDist',600,'MinPeakProminence',0.1);
a = locs_u./locs_d;
b = locs_u1./locs_d;
plot(x_d,dcs_1)
hold on
plot(x_d(locs_d),pks_d, '+r')
hold off

%% 
dcs_1res = resample(dcs_1cm,50,1)
plot(dcs_1res)
hold on
% plot(dcs_1up)

%%
xu = (1:1:length(ecg1))/1000;
xd = (1:1:length(dcs_1lp(1:1200)))/20;

%%
ecg_res = resample(ecg1,1,50,0);
% ecg_res = resample(ecg_res,1,5);
te = (1:1:length(ecg_res))/50;
plot(te,ecg_res,'b')
legend("Resampled","Original")
hold on
plot(xu, ecg1,'r')
%%
ini = dcs_1(locs_d(1):locs_d(2));
cyc =zeros(length(pks)-1,length(ini)); 
cyc(1,:)= ini;
x = (1:1:length(ini))/1000;
count = 0;
avg = ini;
for i=1:1:length(locs_d)-1
    count = count+1;
 
    sig = dcs_1(find(xd==round(xu(locs(i)),1)):find(xd==round(xu(locs(i+1)),1)));
    hold on;
    if length(avg) > length(sig)
        sig(length(sig):length(avg)) = 0;
    elseif length(avg)< length(sig)
        sig = sig(1:length(avg));
    end
    avg = avg+sig;
    plot((1:length(sig))/20,sig, 'DisplayName',""+i+"");
%     legend show
    cyc(count,:) = sig;
        
end
xlabel("Time(s)")
ylabel("aDb value")
title("Marker=ECG R peak, DCS 1cm Baseline Forearm")
hold off
avg = avg/count;

figure()
plot(x,avg)

%Plotting the ensemble average
ensavg = mean(cyc,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cyc,[],1)/sqrt(count);                             % Calculate 95% Confidence Intervals         
figure()
plot(x, ensavg, '-r', 'LineWidth',1)
hold on;
plot(x, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
grid
legend('Ensemble Average', '95% Confidence Intervals')

%%resampling

%%

figure()
hold on
plot(xd,dcs_1/max(dcs_1),'r')
plot(xd,ecg_res/max(ecg),'b')

%% Alligning the peaks of the DCS signal
% The idea is to coincide all the cycles by shifting the peaks to coincide
% with eachother to remove the temporal shift.
Extract=ones(size(pks_ECG_smooth,2)-1,min(diff(locs_ECG_smooth)));
Extract=Extract*NaN;

dcs_1_smooth2=circshift(dcs_1_smooth,shift)

for i=1:size(pks_ECG_smooth,2)-1
    locs_ECG_smooth(i)
    locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))
    signal = dcs_1_smooth2(1,locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1);
    base_sig = dcs_1_smooth2(1,locs_ECG_smooth(1):locs_ECG_smooth(1)+min(diff(locs_ECG_smooth))-1);
    sig = signal + (max(base_sig) - max(signal));
    [pks_b,locs_b] = findpeaks(base_sig, 'MinPeakHeight',0.2);
    [pks_s,locs_s] = findpeaks(sig, 'MinPeakHeight',0.2);
%     if locs_b>=locs_s
%         sig = shi%%%%%%%% Start from here. Subtract the difference betweent base signal location from the current signal
   
    if i==5
        break;
    end
    
%     Extract(i,:)=dcs_1_smooth2(1,locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1);

    Extract(i,:)=sig;
end

%%
bp = 1:100:length(dcs_1lp);
dcs_d = detrend(normalize(dcs_1lp),1,bp);
plot(normalize(dcs_1lp));
hold on;
plot(dcs_d);

%% Signal reconstruction using harmonics
step = 50;
final = zeros(1,step);
for i=1:step
    final = final + squeeze(Data(i+800,4,:));
end
final = final/step;
plot(final)

%%
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
signal = ecg1;
L = length(signal);             % Length of signal
t = (0:L-1)*T;  

Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% calculating the inverse fft to get the amplitude of the signal



f = Fs*(0:(L/2))/L;
figure()
plot(f,P1) 
hold on;
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

m = max(P1(100:length(P1)));


[p_hr,l_hr] = findpeaks(P1(100:length(P1)), "MinPeakHeight", m-0.01);
l_hr = l_hr + 98;
f_hr = Fs*(l_hr)/(L);
scatter(f_hr,p_hr);
hold off;

y = Y;
a = abs(ifft(y,40960));
figure();
plot(1:1:length(a),normalize(a))
% figure();
hold on;
plot(1:1:length(a),ecg1(1:length(a)))



%% CrCP calculation using the Frequency method

% Preprocessing the data
ecg_a = data(datastart(1):dataend(1));
bp_a = data(datastart(2):dataend(2));
tcd_a = data(datastart(3):dataend(3));

ecg1 = ecg_a(1:120000);
ecg1 = normalize(ecg1);
ecg1 = lpf(ecg1,5,1000);
    
tcd = tcd_a(1:length(ecg1));
% tcd = normalize(tcd);
tcd = lpf_ffilts(tcd,30,1000);

bp_a = bp_a(1:length(ecg1));
bp_a = lpf(bp_a,3,1000);

% DCS data processing
% aDb1 = hybrid_dcs(Data,Data_tau);
aDb1 = adb_avg;

dcs_1cm = aDb1(1,:).*10^9;
dcs_1lp = lpf_ffilts(dcs_1cm,15,20);
dcs_15 = aDb1(2,:).*10^9;
dcs_15lp = lpf_ffilts(dcs_15,15,20);
dcs_2 = aDb1(3,:).*10^9;
dcs_2lp = lpf_ffilts(dcs_2,15,20);
dcs_25 = aDb1(4,:).*10^9;
dcs_25lp = lpf_ffilts(dcs_25,15,20);

sig = [ecg1];
p_0 = zeros(1,size(sig,1));
p_f = zeros(1,size(sig,1));
for i=1:size(sig,1)
    signal = sig(i,:);
    Fs = 1000; 
    T = 1/Fs; 
    L = length(signal);             % Length of signal
    t = (0:L-1)*T;  
    
    Y = fft(signal);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    P1_ecg = P1;
    figure();
    f = Fs*(0:(L/2))/L;
    plot(f,P1) 
    hold on;
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    m = max(P1(50:length(P1)))
    [p_sig,l_sig] = findpeaks(P1(50:length(P1)), "MinPeakHeight", m-0.03);
    l_sig = l_sig(1);
    p_sig = p_sig(1);
    p_sig = p_sig+P1(1);
    f_sig = Fs*(l_sig+49)/(L);
    scatter(f_sig,p_sig);
    l_hr = l_sig+49;
    
    y = Y(1:l_hr+150);
    % y = Y;
    a = abs(ifft(y(1:200),4096));
    figure();
    plot(a)

end

%Calculating the amplitude at f(hr) for DCS signal
Fs = 20;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
signal = dcs_25lp;
L = length(signal);             % Length of signal
t = (0:L-1)*T;  

Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P1_DCS = P1;

pDCS = P1(1);
pDCS_f = max(P1((f_sig*L/Fs)-20:(f_sig*L/Fs)+20));

%Calculating the amplitude at f(hr) for ABP signal
signal = tcd(1:length(ecg1));
Fs = 1000;   
T = 1/Fs;  
L = length(signal);             % Length of signal
t = (0:L-1)*T;  

Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P1_TCD = P1;
pTCD = P1(1);
% pTCD_f = P1(l_hr);
pTCD_f = max(P1((f_sig*L/Fs)-20:(f_sig*L/Fs)+20));

%Calculating the amplitude at f(hr) for ABP signal
signal = bp_a(1:length(ecg1));
Fs = 1000;   
T = 1/Fs;  
L = length(signal);             % Length of signal
t = (0:L-1)*T;  

Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P1_abp = P1;
pABP = P1(1);
% pABP_f = P1(l_hr);
pABP_f = max(P1((f_sig*L/Fs)-20:(f_sig*L/Fs)+20));

%%
close all;
t = 0:0.01:1/f_sig;

sin_dcs =pDCS+ pDCS_f*sin(2*pi*f_sig*t);
figure()
a = 0:0.05:1;
plot(a,normalize(dcs_15lp(34:54)));
hold on;
plot(t,normalize(sin_dcs))
title("Fundamental DCS")

sin_tcd =pTCD+ pTCD_f*sin(2*pi*f_sig*t);
figure()
x = 0:0.001:1;
plot(x,normalize(tcd(100:1100)));
hold on;
plot(t,normalize(sin_tcd))
title("Fundamental TCD")

pBP =  (max(bp_a)+ 2*min(bp_a))/3;  % MAP pressure 
sin_bp =pBP+ pABP_f*sin(2*pi*f_sig*t);
figure();
plot(x,normalize(bp_a(100:1100)));
hold on;
plot(t,normalize(sin_bp))
title("Fundamental BP")
figure();
scatter(sin_bp,sin_tcd)
title("Scatter plot of TCD vs BP")



p_tcd = robustfit(sin_bp,sin_tcd);
x= -p_tcd(1)/p_tcd(2):1:max(max(sin_bp,[],2));
y = p_tcd(1)+p_tcd(2)*x;
figure();
plot(x,y);
hold on; 
scatter(sin_bp,sin_tcd)
xlabel("ABP (mmHg)")
ylabel("CBFV (cm/s)")
title("Regression TCD vs BP")
hold off;
ccp_tcd = -p_tcd(1)/p_tcd(2);   % simplifying the linear equation y=m*x + c for y=0 will lead to this equation 
fprintf("\nCcCP TCD is = %d\n",ccp_tcd);

p_dcs = robustfit(sin_bp,sin_dcs);
x= -p_dcs(1)/p_dcs(2):1:max(max(sin_bp,[],2));
y = p_dcs(1)+p_dcs(2)*x;
figure();
plot(x,y);
xlabel("ABP (mmHg)")
ylabel("cBFi (aDb)")
hold on; 
scatter(sin_bp,sin_dcs)
title("Regression DCS vs BP")
hold off;
ccp_dcs = -p_dcs(1)/p_dcs(2);   % simplifying the linear equation y=m*x + c for y=0 will lead to this equation 
fprintf("CrCP DCS is = %d\n",ccp_dcs);

%%
pBP =  (max(bp_a)+ 2*min(bp_a))/3;  % MAP pressure
psig = mean(signal);   % mean value of signal

CrCP = pBP - (pABP_f/pTCD_f)*pTCD;


%% Averaging g2 curves to improve the SNR

g2_1_temp=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(Data(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(Data(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(Data(:,4,:)-1); %g2-1 curve generation

% average g2 curve for 4 channels
g2(1,:,:)=(g2_1_temp);
g2(2,:,:)=(g2_2_temp);
g2(3,:,:)=(g2_3_temp);
g2(4,:,:)=(g2_4_temp);

 

dcs = squeeze(g2(1,:,:));
ecg1 = ecg_a;
ecg1 = normalize(ecg1);
if length(dcs)<length(ecg1)  % shifting only if it is a DCS signal else not so condition is check the length of the signal
    ecg1 = circshift(ecg1,-750);
end
ecg1 = ecg1(1:size(dcs,2)*50);


% bp = 1:100:length(dcs);
% dcs_d = detrend(normalize(dcs),1,bp); 
% [dcs_1,filter_obj] = bandpass(dcs,[0.5 6],20,'ImpulseResponse','fir');
dcs_1 = dcs;
plot(dcs_1)
time_DCS=0.05*(1:1:size(dcs_1,2));
time_ECG=0.001*(1:1:size(ecg1,2));

figure()
hold on
plot(time_DCS,dcs_1/max(dcs_1),'r')
plot(time_ECG,ecg1/max(ecg1),'b')


% Find peaks
[pks_ECG,locs_ECG]=findpeaks(ecg1/max(ecg1),'MinPeakHeight',0.5);

locs_DCS_time=locs_DCS*0.05;
locs_ECG_time=locs_ECG*0.001;

time_shift1=locs_ECG_time(1)-locs_DCS_time(1) %% in s
% Upsampling the DCS singal
x = 1:1:length(dcs_1);
uf = 50;   % Upsampling factor
xq = (1:(1/uf):length(dcs_1)+((uf-1)/uf)); 
dcs_1_interp = interp1(x, dcs_1,xq,'makima');
% dcs_1_interp = dcs_1;


%
fig1=figure('units','centimeters', 'Position',[2 2 35 13]); %18 width 15 heigh
hold on
plot(dcs_1_smooth/max(dcs_1_smooth),'r')
plot(locs_DCS_smooth, pks_DCS_smooth,'*k')

plot(ecg1_smooth/max(ecg1_smooth),'b')
plot(locs_ECG_smooth, pks_ECG_smooth,'*k')

if locs_ECG_smooth(1)>locs_DCS_smooth(1)
    Difference=locs_ECG_smooth(1)-locs_DCS_smooth(1);
    shift = floor(1.75*Difference);
else
    shift =0;
end
% Extracting data

Extract=ones(size(pks_ECG_smooth,2)-3,min(diff(locs_ECG_smooth)));
Extract=Extract*NaN;

dcs_1_smooth2=circshift(dcs_1_smooth,00);

% for i=2:size(pks_ECG_smooth,2)-2
%     locs_ECG_smooth(i)
%     locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))
%     signal = dcs_1_smooth2(locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1);
%   
%     for j=1:
%         temp(1,:,:)=mean(g2(1,floor((locs_ECG_smooth(i-2:i+2))/50),:));
%     end
% 
% end

%% Smoothing the g2 curve with the moving average window for individual curve
g2_smooth = zeros(size(g2));
for i=1:size(g2,2)
    g2_smooth(1,i,:) = smooth(g2(1,i,:),5);
    g2_smooth(2,i,:) = smooth(g2(2,i,:),5);
    g2_smooth(3,i,:) = smooth(g2(3,i,:),5);
    g2_smooth(4,i,:) = smooth(g2(4,i,:),5);
    
end

%% Plotting g2 curves for whole signal
for i=1:5
    semilogx(Data_tau,squeeze(g2_smooth(4,i,:)));
    hold on;
end
%%
% for i=1:size(g2,1)
%     Data_smooth(:,i,:) = g2(i,:,:);
% end
%%
mu= [0.1,0.15,0.17]; %cm^-1 baseline absorption coefficient
mus = 10; 
for i=1:length(mu)

    g2(1,:,:)=squeeze(Data(:,1,:)-1); %g2-1 curve generation
    g2(2,:,:)=squeeze(Data(:,2,:)-1); %g2-1 curve generation
    g2(3,:,:)=squeeze(Data(:,3,:)-1); %g2-1 curve generation
    g2(4,:,:)=squeeze(Data(:,4,:)-1); %g2-1 curve generation

    mua = mu(i);

    Channel=4;
    Curve_no=100;
    rho = [1 1.5 2 2.5];

    beta=g2(Channel,Curve_no,1);
    aDb=aDb1(Channel,Curve_no);
    
    g2_fit=gen_DCS_fit(Data_tau,mua,mus,rho(Channel),beta,aDb);

    semilogx(Data_tau,squeeze(Data(Curve_no, Channel,:))-1,'k')
    hold on
    semilogx(Data_tau,g2_fit)
   

end
%% Calculating the residuals for the fit
% Plots the fitted g2 curve
if exist("g2")
    clear g2
end
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10; 

g2(1,:,:)=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2(2,:,:)=squeeze(Data(:,2,:)-1); %g2-1 curve generation
g2(3,:,:)=squeeze(Data(:,3,:)-1); %g2-1 curve generation
g2(4,:,:)=squeeze(Data(:,4,:)-1); %g2-1 curve generation
g2_fit = zeros(size(g2))*NaN;

for i=1:size(g2,1)
    for j=1:size(g2,2)
        Channel=i;
        Curve_no=j;
        rho = [1 1.5 2 2.5]; 
    
        beta=g2(Channel,Curve_no,1);
        aDb=aDb1(Channel,Curve_no);
        
        g2_fit(i,j,:)=squeeze(gen_DCS_fit(Data_tau,mua,mus,rho(Channel),beta,aDb));
    
%         semilogx(Data_tau,squeeze(Data(Curve_no, Channel,:))-1,'k')
%         hold on
%         semilogx(Data_tau,squeeze(g2_fit(i,j,:)),'r')
%         break
    end
end
residual = (g2_fit-g2).^2; % This is a squared difference between the g2 curve and its fitting
res_error = sum(residual,3);

%% Calcluating the moving average of the g2 curve ans smoothing the signal
g2 = Data;
for i=2:size(Data,1)-1
    for j=1:size(Data,2)
        g2(i,j,:) = mean(g2(i-1:i+1,j,:)); %Moving average of 3 points  || This improves the SNR by 1.1 dB for baseline
%         g2(i,j,:) = (g2(i-1,j,:)+g2(i,j,:)*2+g2(i+1,j,:))./4; % Weighted moving average of g2 curve 
    end
end
adb_smooth = hybrid_dcs(g2,Data_tau);

%% Average of the aDb signal

% for i=2:size(aDb1,2)-1
%     adb_avg(1,i-1) = mean(aDb1(1,i-1:i+1));
% end
% plot(adb_avg,'b');
% hold on;
% plot(adb_smooth(1,:),'r');
% plot(aDb1(1,:),'k')

%% g2 average on consecutive cycle
% Advancing the ECG signal to remove the time delay
ecg_ad = circshift(ecg1,-700);
% Finding the R-R peaks of ECG signal

if exist("g2")
    clear g2
end
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10; 

g2 = Data;
g2_avg = zeros(size(g2))*NaN;



[h_pks,l_pks] = findpeaks(normalize(ecg1),"MinPeakHeight",2.5);
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
aDb1 = hybrid_dcs(Data(1:1200,:,:),Data_tau);

%% Analysing different interpolation methods
dcs_1 = dcs_25lp_avg;

x = 1:1:length(s);
uf = 50;   % Upsampling factor
xq = (1:(1/uf):length(dcs_1)+((uf-1)/uf)); 
dcs_25lp_m = interp1(x, dcs_1,xq,'makima');
dcs_25lp_l = interp1(x, dcs_1,xq,'linear');
dcs_25lp_cb = interp1(x, dcs_1,xq,'cubic');

dcs_25lp_pc = interp1(x, dcs_1,xq,'pchip');
dcs_25lp_v5 = interp1(x, dcs_1,xq,'v5cubic');
dcs_25lp_spl = interp1(x, dcs_1,xq,'spline');

t = (0:1:length(s)-1)./20;
tu = (0:1:length(dcs_25lp_l)-1)./1000;

hold on;
plot(t,normalize(s),'b');
plot(tu,normalize(dcs_25lp_m));
plot(tu,normalize(dcs_25lp_spl));
plot(tu,normalize(dcs_25lp_l));
plot(tu,normalize(dcs_25lp_cb));
plot(tu,normalize(dcs_25lp_pc));
plot(tu,normalize(dcs_25lp_v5));
legend("Original",'makima','spline','linear','cubic','pchip','v5')

%% Signal Distortions mesurement metrics
