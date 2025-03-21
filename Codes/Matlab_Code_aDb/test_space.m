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
%% Creating a GIF file from the Data to showcase the correlation between g2 curve, adb signal, 
close all;
% figure('units','normalized','outerposition',[0 0 1 1])
% adb_avg = aDb1;
chan = 1;
span = 1190:1209;

mx = max(max(g2_n(span,chan,:)));
mn = min(min(g2_n(span,chan,:)));

for i=1:1:length(adb_avg(chan,span))-1
    
    subplot(1,2,2);
    x = 0:1:i;
    
    plot(x,adb_avg(chan,span(1):span(1)+i),'Linewidth',1.5)
    axis([0 22 0.9*min(adb_avg(chan,span)) 1.2*max(adb_avg(chan,span))])
    xlabel("Time(s)",'FontSize',11)
    ylabel("CBFi",'FontSize',11);
    title("DCS 1.5cm signal",'FontSize',11)
    hold on;
    plot(x(i)+1,adb_avg(chan,span(1)+i),'o','MarkerFaceColor','r');
    hold off;

    
    Channel=chan;
    Curve_no=span(1)+i;
    rho = [1 2.5];
    
    beta=g2( Curve_no,Channel,1);
    aDb_fit=aDb1(Channel,Curve_no);
    
    g2_fit=gen_DCS_fit(Data_tau,mua,mus,rho(Channel),beta,aDb_fit);
    
%     semilogx(Data_tau,squeeze(g2(1,1,:)),'k')
%     hold on
    
    sub_f = subplot(1,2,1);
%     semilogx(sub_f,Data_tau,squeeze(Data(span(1)+i,chan,:)),'o','MarkerFaceColor',[0 0.447 0.741])
    semilogx(sub_f,Data_tau,squeeze(g2(span(1)+i,chan,:)),'b','Linewidth',1.5)
    hold on;
    semilogx(sub_f,Data_tau,g2_fit,'r')
    hold off;
    title("Autocorrelation curve",'FontSize',11)
    xlabel("Delay time - tau (s)",'FontSize',11)
    legend("g2","g2 fit")
    ylabel("g_2(tau)",'FontSize',11);
    xlim(sub_f,[10^-6 10^-3]);
    ylim(sub_f,[mn mx]);
    if i==1
        gif('g2_vs_CBFi_SP.gif','DelayTime',1/5);
    end
    gif
end

%% Creating GIF for intensity vs g2 curve
close all;
for j=1:length(Data_tau)
    if j==1
        gif('intensity_vs_tau.gif','DelayTime',1/5);
    end
    subplot(1,2,1)
    plot(squeeze(Data(1:200,1,j)))

    subplot(1,2,2)
    semilogx(Data_tau,squeeze(Data(1,1,:)));
    hold on;
    plot(Data_tau(j),squeeze(Data(1,1,j)),'o','MarkerFaceColor','r')
    hold off;
    gif
end

%% Creating the GIF for the presentaiton about the sigal processing
close all;
m = 1;


  hfig = figure();

    picturewidth = 30; % set this parameter and keep it forever
    hw_ratio =0.5; 
    

    j =1;
    x = (1:1:3600)/20;
    plot (x,aDb1(2,1:3600),'b');
    rectangle("position",[ld_pks(j)/20, 2e-9, ld_pks(30+j), 0.6e-8] )
    xlabel("Time(s)")
    % hold on;
    % plot(adb_avg(2,1:3600),'r')
    set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

    j= 0;
     %%
    hfig = figure();
m = 1
    picturewidth = 20; % set this parameter and keep it forever
    hw_ratio =1.1; 
    set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

    for i = 1:1:min_length
        set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        subplot(3,2,[1 2])
        

        d = 4; % number of points to shift back for the

        t = (1:1:length(aDb1(2,[ld_pks(10)-d:ld_pks(40)-d])))/20;
        plot(t,aDb1(2,[ld_pks(10)-d:ld_pks(40)-d]),'k','Linewidth',1.5)
        xlim([0 30])
        hold on;
        disp(ld_pks(10:40)+i)
        scatter(t(ld_pks(10:39)+i+1-ld_pks(10)),aDb1(2,ld_pks(10:39)+i-d),'r','filled')
        ylabel("BFi cm^2/s)")
        xlabel("Time(s)");
        hold off;
        pause(0.1);
        
        subplot(3,2,3)
        ylim([0.98 1.17])
        xlim([1e-6 1e-3])
        semilogx(Data_tau,squeeze(Data(ld_pks(10:39)+i,4,:)));
        ylabel("g_2")
        xlabel("\tau")
        subplot(3,2,4)
        semilogx(Data_tau,squeeze(mean(Data(ld_pks(10:39)+i,4,:))),'k','Linewidth',1.5);
        xlabel("\tau")
        ylim([0.98 1.17])
        xlim([1e-6 1e-3])

    
        subplot(3,2,[5 6])
%         plot(adb_avg(2,[ld_pks(10):ld_pks(39)]))
        plot(t,adb_avg(2,[ld_pks(10)-4:ld_pks(40)-4]),'k','Linewidth',1.5)
        xlim([0 30])
        t = (1:1:length(adb_avg(2,[ld_pks(10):ld_pks(40)])))/20;
    
        hold on;
        scatter(t(ld_pks(10)+i-ld_pks(10)+3),adb_avg(2,ld_pks(11)+i-3),'r','filled')
        ylabel("BFi (cm^2/s)")
        xlabel("Time(s)")
        hold off;
        

        if m==1
            gif('moving average.gif','DelayTime',1/2);
            m = 0;
        end

        gif
    end

    %%
    close all;

   for i = 1:1:min_length
        subplot(3,2,[1 2])

        d = 4; % number of points to shift back for the

        t = 1:1:length(aDb1(2,[ld_pks(10)-d:ld_pks(40)-d]));
        plot(aDb1(2,[ld_pks(10)-d:ld_pks(40)-d]))
        hold on;
        disp(ld_pks(10:40)+i)
        scatter(t(ld_pks(10:39)+i+1-ld_pks(10)),aDb1(2,ld_pks(10:39)+i-d),'r')
%         hold off;
%         pause(0.1);
%         
%         subplot(3,2,3)
%         ylim([0.98 1.17])
%         semilogx(Data_tau,squeeze(Data(ld_pks(10:39)+i,4,:)));
%         subplot(3,2,4)
%         semilogx(Data_tau,squeeze(mean(Data(ld_pks(10:39)+i,4,:))));
%         ylim([0.98 1.17])
%     
%         subplot(3,2,[5 6])
% %         plot(adb_avg(2,[ld_pks(10):ld_pks(39)]))
%         plot(adb_avg(2,[ld_pks(10)-4:ld_pks(39)-4]))
%         t = 1:1:length(adb_avg(2,[ld_pks(10):ld_pks(40)]));
%     
%             hold on;
%             scatter(t(ld_pks(10)+i-ld_pks(10)+3),adb_avg(2,ld_pks(11)+i-3),'r')
%         hold off;
% 
%         if m==1
%             gif('moving average.gif','DelayTime',1/2);
%             m = 0;
%         end
%         gif
    end


        %%
        close all
        plot(adb_avg(2,[ld_pks(10)+i-3:ld_pks(39)+i-3]))
        for i = 1:20
          t = 1:1:length(adb_avg(2,[ld_pks(10):ld_pks(40)]));
            
            hold on;
            scatter(t(ld_pks(10)+i-ld_pks(10)+3),adb_avg(2,ld_pks(11)+i-3),'r')
        end
%%
a = (1:1:length(1:100))/20;
plot(a,squeeze(Data(1:100,1,1)))
xlabel("Time(s)")
ylabel("Intensity")
title("DCS Intensity plot")
%%
close all;
x = (1:1:size(adb_avg,2))/20;
plot(x,normalize(adb_avg(2,:)),'Linewidth',1.5);
xlabel("Time(s)")
ylabel("Blood Pressure and CBFi (Normalized unit)")
hold on;
% plot(x,normalize(adb_avg(4,:)));
x_u = (1:1:size(bp,2))/1000;
plot(x_u,normalize(bp),'Linewidth',1.5);
% plot(x_u,normalize(tcd),'k','Linewidth',1.5)
legend("BFI(r_s_d=1.5cm)", "ABP");

%%
close all;

plot(x_u,bp,x,circshift(adb_avg(1,:).*10^10,10),'Linewidth',1.5)
xlabel("Time(s)")
ylabel("ABP(mmHg) and CBFi*10^1^0(cm^2/s)")
title("ABP & DCS(r_s_d=1.5cm) waveform")
legend("ABP","BFI(r_s_d=1.5cm)");
%% Plotting TCD and BP
close all;
x_u = (1:1:size(bp,2))/1000;
plot(x_u,(bp),x_u,tcd,'Linewidth',1.5);
% hold on;
% plot(x_u,(tcd),'r','Linewidth',1.5)
xlabel("Time(s)");
ylabel("ABP(mmHg) and CBFV(cm/s)")
legend("ABP","TCD");

%% Comparison DCS 1.5cm vs 2.5cm signal
close all;
x = (1:1:size(adb_avg,2))/20;
plot(x,normalize(adb_avg(2,:)));
xlabel("Time(s)")
ylabel("CBFi")
hold on;
% plot(x,normalize(adb_avg(4,:)));
plot(x,normalize(adb_avg(4,:)));
legend("BFI(r_s_d=1.5cm)", "BFI(r_s_d=2.5cm)");

%% Plotting a regression plot
p = robustfit(bp(44938:45915)',dcs_1_smooth2(44938:45915)');
x= -10:1:max(max(bp(1:700),[],2));
y = p(1)+p(2)*x;
figure();
plot(x,y);
hold on; 
scatter(bp(44938:45915)',dcs_1_smooth2(44938:45915)')
hold off;
ccp(i) = -p(1)/p(2);   % simplifying the linear equation y=m*x + c for y=0 will lead to this equation 
fprintf("%d",ccp(i));

%% Data analysis and plotting from CrCP values from MPCM study

CCP = [26.0259	27.3766	60.98740877	14.1122	63.8348	10.3636	22.7256	12.2824	42.2544
20.6735	38.9178	51.72297061	33.6915	61.6646	13.2007	19.6981	12.2742	29.605
27.2814	25.6519	54.67380382	33.1734	51.5185	26.6407	6.2934	11.2485	31.0852
15.2388	24.2512	23.11308583	31.756	37.0627	13.1352	10.2829	-9.1616	31.633
18.7298	18.2728	19.60785249	30.8903	61.8123	11.2472	-2.4666	-13.0602	43.5967
]; %CrCP values table column represents each subjects respectively from 1 to 9 and row represents TCD;DCS:1cm,1.5cm,2cm,2.5cm values

mean_ccp = mean(CCP,2);
std_ccp = std(CCP,1,2);
[rho1,p1]  = corrcoef(CCP(1,:),CCP(2,:));
[rho15,p15]  = corrcoef(CCP(1,:),CCP(3,:));
[rho2,p2]  = corrcoef(CCP(1,:),CCP(4,:));
[rho25,p25]  = corrcoef(CCP(1,:),CCP(5,:));



%% Plotting the signal with the points.

t_res=0.05; % seconds
t_avg=1; % window used for averaging in seconds

t_avg_pt=t_avg/t_res; % window used for averaging in points
time=t_res*(1:1:size(aDb1,2));
Data_time(1,i)=i*t_res;
j=1;
for i=2:t_avg_pt:size(aDb1,2)
    aDb1_avg(1,j)=mean(aDb1(1,i-1:i+t_avg_pt-2));
    aDb1_avg(2,j)=mean(aDb1(2,i-1:i+t_avg_pt-2));
%     time_avg(1,j)=(Data_time(1,i+t_avg_pt-1));
    j=j+1;
end


t = (1:1:length(aDb1(1,:)))/20;

x_pos=[2,3,4,5,6,7,8,9]; %task strat time in minutes
txt = ["10 mmHg","30 mmHg","50 mmHg","70 mmHg","90 mmHg","120 mmHg","150 mmHg","BSL"];
for i=1:size(x_pos,2)
if mod(i,2)==0
    rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,60,max(aDb1(1,:))],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
        'LineWidth',3)
    text(t(x_pos(i)*1200)*1.05,0.9*max(aDb1(1,:)),txt(i));
else
    rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,60,max(aDb1(1,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
        'LineWidth',3)
    text(t(x_pos(i)*1200)*1.05,0.9*max(aDb1(1,:)),txt(i));

end
end

hold on;
plot(t,aDb1','DisplayName','aDb1');
hold on;
plot(t(120*20),aDb1(1,2400),'r*');
title("Tourniquet Pressure Levels vs CBFi");
xlabel("Time(s)");
ylabel("CBFi");
plot(aDb1_avg','DisplayName','aDb1_avg','LineWidth',1.5)

%% Extracting the data from the figures
open('C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\Pulsatility Analysis\DCS_10_S1_BSL.fig');
a = get(gca,'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');


%% Testing the g2 averaging

close all;
if exist("g2")
    clear g2;
end
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10; 

% g2 = Data;

g2_n = Data_all;
strt_time =  [1,190,310,430,540];   %time in seconds
stp_time =  [160,290,410,530,660];    %time in seconds.


% for m=1:length(strt_time)

g2 = Data_all(1:3600,1:4,:);
ecg = ecg_ad(1:180*1000);


[h_pks,l_pks] = findpeaks(normalize(ecg),"MinPeakHeight",2.5,'MinPeakDistance',750);

fig1=figure('units','centimeters', 'Position',[2 2 35 13]); %18 width 15 heigh
hold on
plot(normalize(ecg),'r')
plot(l_pks, h_pks,'*k')
hd_pks = floor(h_pks./50);
ld_pks = floor(l_pks./50);
% g2_avg = zeros(size(g2,1),size(g2,2),size(g2,3));
avg_window_width = 40;
for i=1:size(ld_pks,2)-1
    if size(ld_pks,2)-i < avg_window_width
        avg_window_width = avg_window_width-1;
    end
    min_length = min(diff(ld_pks(i:i+avg_window_width)));
    base_sig = g2(ld_pks(i):ld_pks(i)+min_length,:,:);
   
    for j=1:avg_window_width-1
        base_sig = base_sig + g2(ld_pks(i+j):ld_pks(i+j)+min_length,:,:);
%         adb_1 = hybrid_dcs(base_sig,Data_tau);
       
    end
    g2(ld_pks(i):ld_pks(i)+min_length,:,:) = base_sig./avg_window_width;
end

% end

% for i=1:size(g2,1)
%     g2(i,2,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
% end



adb_avg = standalone_dcs(g2,Data_tau);
% SNR(1,k) = snr(adb_avg(2,:),20);
figure
snr(adb_avg(2,:),20)

% adb_avg = adb_avg; % This is cut out specific portion of the signal to plot
% figure();
adb = standalone_dcs(Data_all(1:3600,1:4,:),Data_tau);
dcs_1lp = lpf(adb_avg(1,:),7,20);
dcs_25lp = lpf(adb_avg(2,:),7,20);
% dcs_25lp_tr = lpf(adb_avg(3,:),7,20);

% aDb1 = aDb1;
% figure();
% snr(adb_avg(2,:),20);
figure();
plot(adb_avg(2,:),'b',"LineWidth",1.5);
hold on; 
plot(adb(2,:),'r');
legend("Ensemle Temporal Averaged","Raw")
% title("Comparision of g2 averaging for cuff data MPCM004 width=50 cycles")
% legend("g2 Averaged signal","Raw signal")
% xlabel("samples (Time = samples/20)");
% ylabel("aDb")

%% Pressure Modulation test
dcs = Data;
g2(1,:,:)=squeeze(dcs(:,1,:)-1); %g2-1 curve generation
g2_1_temp=squeeze(dcs(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(dcs(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(dcs(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(dcs(:,4,:)-1); %g2-1 curve generation

% average g2 curve for large source detector separation
% for i=1:size(g2,2)
%     g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:)+g2_5_temp(i,:)+g2_6_temp(i,:)+g2_7_temp(i,:)+g2_8_temp(i,:))/8;
% end
for i=1:size(g2,2)
    g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
end

%%
aDb = standalone_dcs(Data,Data_tau);
%% Normalize the g2 curve
g2(g2<0) = 0;
%%
g2_0_rl = squeeze(mean(g2(2,1:500,:),2)); 
g2_0_rs = squeeze(mean(g2(1,1:500,:),2)); 
g2_p_rl = squeeze(mean(g2(2,630:1000,:),2)); 
g2_p_rs = squeeze(mean(g2(1,630:1000,:),2)); 

%%
close all;
del_OD_l_p = -log((g2_p_rl')./(g2_0_rl'));
del_OD_s_p = -log((g2_p_rs')./(g2_0_rs'));

sens_fact = abs(del_OD_l_p./del_OD_s_p);
mean(sens_fact)
plot(sens_fact)

%% 
% g2(1,:,:) = Data(:,1,:);
% g2(2,:,:) = Data(:,4,:);
close all;
clear all;
subject = 1;

% filename=strcat(pwd,'\2 layer model\2 layer model\00',num2str(subject),'_CUFF_g2.mat');    
filename = "D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\DCS\20220608 - 1\Data.mat";
load(filename);
g2_t(1,:,:) = Data(:,1,:);
for i=1:size(Data,1)
    g2_t(2,i,:) = (Data(i,2,:)+Data(i,3,:)+Data(i,4,:))/3;
end

% g2_t(1,:,:) = g2(1,:,:)+1;
% g2_t(2,:,:) = g2(4,:,:)+1;

% g2_t(:,1,:) = g2(1,:,:)+1;
% g2_t(:,2,:) = g2(2,:,:)+1;
% g2_t(:,3,:) = g2(3,:,:)+1;
% g2_t(:,4,:) = g2(4,:,:)+1;
data = g2_t;
load('D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Matlab_Code_aDb\2 layer model\data_tau.csv')
tau = data_tau;


time=(1:1:size(g2_t,2))*0.05;

%Since this data loaded is in the g2 format, following code convert it into


windowsize=5;
wagi=ones(1,windowsize)/windowsize;
for r=1:2
    for i=1:size(data,3)
        data_temp=squeeze(data(r,:,:));
        data(r,:,i)=filtfilt(wagi,1,data_temp(:,i));
    end
end
% windowsize=2;
% wagi=ones(1,windowsize)/windowsize;
% for r=1:4
%     for i=1:size(data,3)
%         data_temp=squeeze(data(:,r,:));
%         data(:,r,i)=filtfilt(wagi,1,data_temp(:,i));
%     end
% end
% semilogx(tau,squeeze(g2_t(1,100,:)),'b')
% hold on;
% semilogx(tau,squeeze(data(1,100,:)),'r')

% Following code compares the adb fitting from raw g2 and filtered g2
% adb = hybrid_dcs(g2_t,tau);
% adb_filt = hybrid_dcs(data,tau);
% close all;
% plot(adb(1,:)); hold on; plot(adb_filt(1,:));

% 
bsl = 400:500;
cuff = 1400:1500;
%  bsl = [20 25];
%     bsl = bsl/0.05;
%     cuff = [70 75];
%     cuff = cuff/0.05;

g2_s0=mean(squeeze(data(1,bsl,:)))-1;
g2_sp=mean(squeeze(data(1,cuff,:)))-1;
g2_s0=(g2_s0)/(max(g2_s0));
g2_sp=(g2_sp)/(max(g2_sp));



g2_l0=mean(squeeze(data(2,bsl,:)))-1;
g2_lp=mean(squeeze(data(2,cuff,:)))-1;
g2_l0=(g2_l0)/(max(g2_l0));
g2_lp=(g2_lp)/(max(g2_lp));



DeltaOD_s=real(-log(g2_sp./g2_s0));
DeltaOD_l=real(-log(g2_lp./g2_l0));

Ratio=DeltaOD_l./DeltaOD_s;

nanmean(Ratio(1,20:29))
plot(Ratio)

%%
data = g2_t;

g2_s0=mean(squeeze(data(1,bsl,:)));
g2_sp=mean(squeeze(data(1,cuff,:)));
g2_s0=(g2_s0)/(max(g2_s0));
g2_sp=(g2_sp)/(max(g2_sp));



g2_l0=mean(squeeze(data(2,bsl,:)));
g2_lp=mean(squeeze(data(2,cuff,:)));
g2_l0=(g2_l0)/(max(g2_l0));
g2_lp=(g2_lp)/(max(g2_lp));



DeltaOD_s=real(-log((g2_sp-1)./(g2_s0-1)));
DeltaOD_l=real(-log((g2_lp-1)./(g2_l0-1)));

Ratio_1=DeltaOD_l./DeltaOD_s;

nanmean(Ratio_1(1,20:29))
plot(Ratio_1)


%%

a = [2.7	1.89	0	1.5	0.305
4.02	0.82	0	0.73	0.289
2.9	2.7	0	2.1	0.291
0.9	0.869	0	0.78	0.298
5.39	3.3	0	2.1	0.311
];
a(:,[1,2,4]) = a(:,[1,2,4])/10;

b = [1.03	0.478	0	0.49	0.332
0.81	0.35	0	0.144	0.313
5.5	4.12	0	3.2	0.307
10.7	2.75	0	1.83	0.303
8.1	5.2	0	3.3	0.328
];
xdt = [0.305,0.289,0.291,0.298,0.311];
ydt = [0.332,0.313,0.307,0.303,0.328];
[h,p,ci,stats] = ttest2(xdt,ydt);

xdx = a(:,1);
ydx = b(:,1);

for i = 1:size(a,2)
    x = a(:,i);
    y = b(:,i);
    [h(i),p(i)] = ttest2(x,y,'Alpha',0.05);
end

bfi_pp = mean(a(:,5));
bfi_pp_sd = std(a(:,5));
dfc_pp = mean(b(:,5));
dfc_pp_sd = std(b(:,5));


%%
% a = [0.7
% 0.875621891
% 0.782758621
% 0.8
% 0.612244898
% ];

% a = [2.7
% 4.02
% 2.9
% 0.9
% 5.39
% ]*1e-9;
a = [0.715277778
1.396551724
4.296875
1.221978022
1.087248322
];
% b =[0.512621359
% 0.654320988
% 0.749090909
% 0.665467626
% 0.641975309
% ] ;
% b = [1.03
% 0.81
% 5.5
% 5.56
% 8.1
% ]*1e-8;

b = [0.602678571
0.743068392
0.84057971
0.294117647
0.965085049
];
[h,p] = ttest2(a,b)
v = [a b]
boxplot(v)

%%
a = [63
58
78
64
72
];
b = [40
43
74
56
67
];
c = [4
3
4
7
6
];

[h1,p1] = ttest2(a,c);
[h2,p2] = ttest2(b,c);


%%
load('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\important variables\pulsatility comparison BFi DFC\6.mat')
s1 = sig;
load('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\important variables\pulsatility comparison BFi DFC\8.mat')
s2 = sig;
load('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\important variables\pulsatility comparison BFi DFC\01.mat')
s3 = sig;
load('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\important variables\pulsatility comparison BFi DFC\03.mat')
s4 = sig;
load('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\important variables\pulsatility comparison BFi DFC\4.mat')
s5 = sig;


%%
close all;
plot(s1');
hold on;
plot(s2');
plot(s3');
plot(s4');
plot(s5');

%% 
close all;
bfi = [s1(1,1:1000); s2(1,1:1000); s3(1,1:1000); s4(1,1:1000); s5(1,1:1000)];
bfi_mean = mean(bfi);
bfi_std = std(bfi);
dfc =  [s1(2,1:1000); s2(2,1:1000); s3(2,1:1000); s4(2,1:1000); s5(2,1:1000)];
dfc_mean = mean(dfc);
dfc_std = std(dfc);

plot(bfi_mean,'k')
hold on;
plot(bfi_mean+bfi_std,'k--')
plot(bfi_mean-bfi_std,'k--')
plot(dfc_mean,'r')
plot(dfc_mean+dfc_std,'r--')
plot(dfc_mean-dfc_std,'r--')

%% Plotting a box plot for the pulsatility comparison
close all;

% The a,d are the deltaT for the BFi and CBFi correspondingly. b,e are
% normalized amplitude of S2 and c,f are normalized amplitude of d
% (dicrotic notch) for BFi and CBFi respectively. Yb and Yc are the
% amplitudes of systolic to diastolic peak
yb = [0.602678571
0.743068392
0.84057971
0.294117647
0.965085049
]
yc = [0.715277778
1.396551724
4.296875
1.221978022
1.087248322
]
a = [0.305
0.289
0.295
0.298
0.311
];
b = [0.7
0.875621891
0.782758621
0.8
0.612244898
]
c = [0.555555556
0.629353234
0.724137931
0.711111111
0.38961039
]

d = [0.332
0.313
0.307
0.303
0.328
]
e = [0.512621359
0.654320988
0.749090909
0.665467626
0.641975309
]
f = [0.475728155
0.597530864
0.581818182
0.485611511
0.407407407
]

hfig = figure();
ysd = boxplot([yb,yc],'Labels',{'BFi','CBFi'},'Colors',['b','r'])
set(ysd,'LineWidth', 2);
% Xlabel("\Deltat (s)")
ylabel("\Deltat (s)")
picturewidth = 20; % set this parameter and keep it forever
hw_ratio =0.75; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])



hfig = figure();
dt = boxplot([a,d],'Labels',{'BFi','CBFi'},'Colors',['b','r'])
set(dt,'LineWidth', 2);
% Xlabel("\Deltat (s)")
ylabel("\Deltat (s)")
picturewidth = 20; % set this parameter and keep it forever
hw_ratio =0.75; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])


hfig = figure();
s2 = boxplot([b,e],'Labels',{'BFi','CBFi'},'Colors',['b','r'])
set(s2,'LineWidth', 2);
ylabel("S2 (A.U)")
picturewidth = 20; % set this parameter and keep it forever
hw_ratio =0.75; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

hfig = figure();
d = boxplot([c,f],'Labels',{'BFi','CBFi'},'Colors',['b','r'])
set(d,'LineWidth', 2);
ylabel("d (A.U.)")
picturewidth = 20; % set this parameter and keep it forever
hw_ratio =0.75; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%% aDb calculations
% close all;
% clear all;
% load("D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\Donya\test2\jignesh\thenar\Data.mat")
dcs = Data;
g2(1,:,:)=squeeze(dcs(:,1,:)-1); %g2-1 curve generation
g2_1_temp=squeeze(dcs(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(dcs(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(dcs(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(dcs(:,4,:)-1); %g2-1 curve generation
% g2_5_temp=squeeze(dcs_t(:,1,:)-1);
% g2_6_temp=squeeze(dcs_t(:,2,:)-1);
% g2_7_temp=squeeze(dcs_t(:,3,:)-1);
% g2_8_temp=squeeze(dcs_t(:,4,:)-1);
% average g2 curve for large source detector separation
% for i=1:size(g2,2)
%     g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:)+g2_5_temp(i,:)+g2_6_temp(i,:)+g2_7_temp(i,:)+g2_8_temp(i,:))/8;
% end
for i=1:size(g2,2)
    g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
end

rho = [0.7 2.5]; %source detector separations in cm 
mua = 0.15; %cm^-1 baseline absorption coefficient
mus = 10; %cm^-1 baseline reduced scattering coefficient

tau_values=Data_tau;

for chan=1:size(g2,1)
% for chan=1:1
     for i=1:size(g2,2) 
        rsd=rho(chan);
        g2_temp(i,:)=squeeze(g2(chan,i,:));
        LB = [0];
        UB = [inf];
        Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
        beta= squeeze(g2(chan,i,1)); %0.1568;
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rsd,beta);
        aDb1(chan,i) = FittedParams(1);
    end
end

t = (1:1:length(aDb1))/20;
plot(t,aDb1')
% adb_lpf = (lpf(aDb1,7,50));
hold on;
% plot(t,adb_lpf)
xlabel("TIme (s)");
ylabel("\alphaD_b")
% title("Baselilne 50ms DCS 2.5cm")
title("New Pressure Cuff vs \alphaD_b")
legend("1 cm","2.5 cm")

%% 
close all;
% adb_lpf = adb_lpf(1:6000);
adb_5 = circshift(adb_50,20)
t2 = (1:1:length(adb_lpf))/50;
t5 = (1:1:length(adb_50))/20;

plot(t2, aDb1,'b');
hold on;
plot(t5, adb_50,'r');
xlabel("Time (s)");
ylabel("\alphaD_b")
legend("20ms Res","50ms Res");

%%
Fs = 20;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
signal = aDb1(2,3600:4800);
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