%% Loading the file
% clear all;
% filename=strcat('D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Data\Test Data DCS baseline\brachial baseline ECG BP.mat');
% load(filename)

%% Upscaling the data by 3
ecg_a = data(datastart(1):dataend(1));
bp_a = data(datastart(2):dataend(2));
tcd_a = data(datastart(3):dataend(3));

ecg1 = ecg_a(1:120000);
ecg1 = normalize(ecg1);
ecg1 = lpf(ecg1,5,1000);
    
tcd = tcd_a(1:length(ecg1));
% tcd = normalize(tcd);
tcd = lpf_ffilts(tcd,30,1000);

bp = bp_a(1:length(ecg1));
bp = lpf_ffilts(bp,40,1000);
%% Plotting the frequency spectrum
Fs = 20;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
signal = adb_all(2,:);
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
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% finding the maxima to find the individual signals
y = ecg1;
x = (1:length(ecg1));
[pks,locs] = findpeaks(y, 'MinPeakHeight', 0.5,'MinPeakDist',500,'MinPeakProminence',0.1);  %Determine peaks and Indices
figure()
plot(x,y)
hold on
plot(x(locs),pks, '+r')
hold off
grid

for k1 = 1:numel(locs)-0.1
    yc{k1} = y(locs(k1)-60:locs(k1+1)-60);                            % Define MUAP Frames
    xc{k1} = x(locs(k1)-60:locs(k1+1)-60);
end

figure()
hold all 
for k1 = 1:numel(yc)
   plot(xc{k1}-xc{k1}(1),yc{k1})
end
hold off
grid

% — CALCULATE & PLOT ENSEMBLE AVERAGE —                                                                    
minlen = min(cellfun(@numel, yc));                                     % Minimum Length Of MUAP Records
ens = zeros(minlen, numel(yc));                                        % Preallocate
for k1 = 1:numel(yc)
    ens(:,k1) = yc{k1}(1:minlen);                                      % Trim MUAPs To Shortest Length
end
ensavg = mean(ens,2);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(ens,[],2)/sqrt(numel(yc));                             % Calculate 95% Confidence Intervals
eatv = mean(diff(x))*(0:minlen-1);                                     % Time Vector For Ensemble AVerage
figure()
plot(eatv, ensavg, '-r', 'LineWidth',1)
hold on
plot(eatv, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(eatv, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
grid
xlabel('Time (s)')
ylabel('Amplitude')
title("Ensemble average of ECG signal")
legend('Ensemble Average', '95% Confidence Intervals')

%Repeating the waveform
ensavg = [ensavg; ensavg];
ecg_ens = ensavg*1000;
if length(ecg_ens)<=1600
    ecg_ens(length(ecg_ens):1600) = 1;
else
    ecg_ens = ecg_ens(1:1600);
end
ecg_ens = ecg_ens';
%Saving the variable
% writematrix(ensavg,'D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Results and Plots\output_variables\Marianne\ecg_ens.csv','Delimiter','comma');

%% TCD signal Processing
tcd = tcd_a(1:length(ecg1));
% tcd = normalize(tcd);
tcd = lpf_ffilts(tcd,15,1000);
%% Plotting the TCD signal
minima = islocalmin(tcd,'MinProminence',2,'MinSeparation',950 );
x = 1:length(minima);
plot(x,tcd,x(minima),tcd(minima),'r*');

l = length(tcd);
ini = tcd(200:1149);

cyc =zeros(sum(minima==1)-1,950); 
% cyc(1,:)= ini;
count = 1;
avg = ini;
for i=1:1:length(minima)
    if (minima(i)==1) && (i+949<=length(minima))
        count = count+1;
        plot(tcd(i:i+949))
        hold on
        avg = avg+tcd(i:i+949);
        cyc(count,:) = tcd(i:i+949);
    end
 
end
hold off
avg = avg/count;
x = (1:1:length(avg))/1000;
figure()
plot(avg)

%Plotting the ensemble average
ensavg = mean(cyc,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cyc,[],1)/sqrt(count);                             % Calculate 95% Confidence Intervals         
figure()
plot(x, ensavg, '-r', 'LineWidth',1)
hold on
plot(x, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
grid
legend('Ensemble Average', '95% Confidence Intervals')
xlabel('Time (s)')
ylabel('BFi')
title("Ensemble average of TCD signal")
ensavg = [ensavg ensavg];

tcd_ens = ensavg;
%Saving the variable
% writematrix(ensavg,'D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Results and Plots\output_variables\jig\tcd_ens.csv','Delimiter','comma');

%% Processing the blood pressure data
bp = bp(1:length(ecg1));
bp = lpf(bp,3,1000);
bp = bp(1:length(tcd));
bp = bp;
minima = islocalmin(bp,'MinProminence',20);
x = 1:length(minima);
plot(x,bp,x(minima),bp(minima),'r*');

l = length(bp);
ini = bp(200:999);

cyc =zeros(sum(minima==1)-1,800); 
cyc(1,:)= ini;
count = 1;
avg = ini;
for i=1000:1:length(minima)
    if (minima(i)==1) && (i+799<=length(minima))
        count = count+1;
        plot(bp(i:i+799))
        hold on
        avg = avg+bp(i:i+799);
        cyc(count,:) = bp(i:i+799);
    end
end
hold off
avg = avg/count;
x = (1:1:800)/1000;
figure()
plot(avg)

%Plotting the ensemble average
ensavg = mean(cyc,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cyc,[],1)/sqrt(count);                             % Calculate 95% Confidence Intervals         
figure()
plot(x, ensavg, '-r', 'LineWidth',1)
hold on
plot(x, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
grid
legend('Ensemble Average', '95% Confidence Intervals')
xlabel('Time (s)')
ylabel('BFi')
title("Ensemble average of TCD signal")
ensavg = [ensavg ensavg];

abd_ens = ensavg;

%Saving the variable
% writematrix(ensavg,'D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Results and Plots\output_variables\marianne\bp_ens.csv','Delimiter','comma');
%% Processing the DCS data

%% Loading the data for Standalone DCS system
% filename=strcat('D:\Jignesh\MSc Western Uni\Research MSc\Data\DCS\20211207\Data.mat');
% load(filename)

% g2(1,:,:)=squeeze(Data(:,1,:)-1); %g2-1 curve generation
% g2_2_temp=squeeze(Data(:,2,:)-1); %g2-1 curve generation
% g2_3_temp=squeeze(Data(:,3,:)-1); %g2-1 curve generation
% g2_4_temp=squeeze(Data(:,4,:)-1); %g2-1 curve generation



%This changed for g2 averaging
g2(1,:,:)=squeeze(Data_avg(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(Data_avg(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(Data_avg(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(Data_avg(:,4,:)-1); %g2-1 curve generation

% average g2 curve for large source detector separation
for i=2:size(g2,2)
    g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
end

% aDb calculation
rho = [1 3]; %source detector separations in cm 
mua = 0.1; %cm^-1 baseline absorption coefficient
mus = 10; %cm^-1 baseline reduced scattering coefficient

tau_values=Data_tau;

for chan=1:size(g2,1)
    for i=1:size(g2,2)
        rsd=rho(chan);
        g2_temp(i,:)=squeeze(g2(chan,i,:));
        LB = [0];
        UB = [inf];
        Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
        beta= squeeze(g2(chan,i,1)); %0.1568;
        beta= squeeze(mean(g2(chan,1:500,1))); %0.1568;
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rsd,beta);
        aDb1(chan,i) = FittedParams(1);
    end
end

dcs_1 = aDb1(1,:).*10^9;
dcs_3 = aDb1(2,:).*10^9;


% Filtering the singal
dcs_1lp = lpf(dcs_1,5,20);
dcs_3lp = lpf(dcs_3,5,20);

%% plotting
Chan=1;
aDb=aDb1(chan,1);

rho=3;

ua = mua; %cm^-1 baseline absorption coefficient
us = mus; %cm^-1 baseline reduced scattering coefficient
tau=tau_values;

v=3;
zo = 1/(us); %first point source term in cm
% zo = 1/(ua+us); %first point source term in cm
D = 1/(3*(us)); %diffusion coefficient in cm
% D = 1/(3*(ua+us));%diffusion coefficient in cm
zb = 2*D*(1+0.493)/(1-0.493); %2nd point source term in cm
r1 = sqrt(rho.^2+zo.^2);
r2 = sqrt((zo+2*zb).^2+(rho.^2));
n = 1.4;%refractive index of the medium
c=v/n;
lamda = 780e-7; %786.5e-007;%wavelength of the light in cm
k = (2*pi*n)/lamda;%Wavenumber of light in the medium
k_D = sqrt((3*us*ua)+(6*us^2*k.^2*aDb*tau));
Amp = (3*us)/(4*pi);
% beta = 0.15;
% Field autocorrelation function
G1 = Amp.*((exp(-k_D.*r1)./r1)-(exp(-k_D.*r2)./r2));
g1 = (G1)./max(G1);
g2_1_fit = beta.*(g1.^2);

semilogx(tau_values,squeeze(g2(2,3,:)))
hold on
semilogx(tau_values,g2_1_fit)


%% upsampling the data

% Don't use interp, instead use interp1


%% Processing the hybrid DCS system data

% aDb1 = hybrid_dcs(Data,Data_tau);

% filename=strcat('D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Data\Test Data DCS baseline\20211214-4\','Data.mat');
% load(filename)



%

g2(1,:,:)=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2(2,:,:)=squeeze(Data(:,2,:)-1); %g2-1 curve generation
g2(3,:,:)=squeeze(Data(:,3,:)-1); %g2-1 curve generation
g2(4,:,:)=squeeze(Data(:,4,:)-1); %g2-1 curve generation

% aDb calculations

rho = [1 1.5 2 2.5]; %source detector separations in cm 
mua = 0.1; %cm^-1 baseline absorption coefficient
mus = 10; %cm^-1 baseline reduced scattering coefficient

tau_values=Data_tau;

for chan=1:size(g2,1)
    for i=1:size(g2,2)
        disp(i);
        rsd=rho(chan);
        g2_temp(i,:)=squeeze(g2(chan,i,:));
        LB = [0];
        UB = [inf];
        Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
        beta= squeeze(g2(chan,1:i,1)); %0.1568;
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rsd,beta);
        aDb1(chan,i) = FittedParams(1);
    end
end

Chan=1;
aDb=aDb1(chan,1);

rho=3;

ua = mua; %cm^-1 baseline absorption coefficient
us = mus; %cm^-1 baseline reduced scattering coefficient
tau=tau_values;

v=3;
zo = 1/(us); %first point source term in cm
% zo = 1/(ua+us); %first point source term in cm
D = 1/(3*(us)); %diffusion coefficient in cm
% D = 1/(3*(ua+us));%diffusion coefficient in cm
zb = 2*D*(1+0.493)/(1-0.493); %2nd point source term in cm
r1 = sqrt(rho.^2+zo.^2);
r2 = sqrt((zo+2*zb).^2+(rho.^2));
n = 1.4;%refractive index of the medium
c=v/n;
lamda = 780e-7; %786.5e-007;%wavelength of the light in cm
k = (2*pi*n)/lamda;%Wavenumber of light in the medium
k_D = sqrt((3*us*ua)+(6*us^2*k.^2*aDb*tau));
Amp = (3*us)/(4*pi);
% beta = 0.15;
% Field autocorrelation function
G1 = Amp.*((exp(-k_D.*r1)./r1)-(exp(-k_D.*r2)./r2));
g1 = (G1)./max(G1);
g2_1_fit = beta.*(g1.^2);

semilogx(tau_values,squeeze(g2(2,3,:)))
hold on
semilogx(tau_values,g2_1_fit)

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
dcs_1lp = lpf_ffilts(dcs_1cm,7,20);
dcs_15 = aDb1(2,:).*10^9;
dcs_15lp = lpf_ffilts(dcs_15,7,20);
dcs_2 = aDb1(3,:).*10^9;
dcs_2lp = lpf_ffilts(dcs_2,7,20);
dcs_25 = aDb1(4,:).*10^9;
dcs_25lp = lpf_ffilts(dcs_25,7,20);



%% Shifting the signal to time allign the DCS signal to ECG
break_pt = 1:100:size(tcd,2);
bp_shift = circshift(bp,-1180)
[a,l_bp] = findpeaks(normalize(bp_shift),"MinPeakHeight",1.5,'MinPeakDistance',500) ;
[a,l_tcd] = findpeaks(normalize(detrend(tcd,1,break_pt)),"MinPeakHeight",1);
shift = l_bp(1)-l_tcd(1);
tcd_shift = circshift(tcd,-100)

%% Calculating the Critical closing pressure
close all;
l = 2400;
ccp_dcs  = ccp_measure_tail(ecg1(1:l*50),dcs_2lp(1:l),bp_shift(1:l*50),60);
% ccp_dcs  = ccp_measure(ecg1,dcs_25lp,bp,10);
% close all;
% scatter(1:length(ccp_tcd),ccp_tcd,'red');
% hold on;
% scatter(1:length(ccp_dcs),ccp_dcs,'blue');
plot(ccp_dcs);
ylabel("CrCP (mmHg)");
title("CrCP using Tail, TCD averaged over 25 cycles")
% close all
%%
load ccp_var_stack.mat;
bp_stack = stack(1:floor(length(stack(:,1))/2),:);
% bp_stack = bp_stack.*0.6; % validate it from the paper
sig_stack = stack(floor(length(stack(:,1))/2)+1:length(stack(:,1)),:);
%%
[pks_signal,locs_signal] = max(sig_stack');
[pks_bp,locs_bp] = max(bp_stack');
% bp_stack = bp_stack(:,locs_bp:length(bp_stack(1,:)));
% sig_stack = sig_stack(:,locs_bp:length(sig_stack(1,:)));
ccp = zeros(1,length(bp_stack(:,1)));
for i=1:length(locs_bp)
    bp = bp_stack(i,locs_bp(i)+300:length(bp_stack(1,:)));
    sig = sig_stack(i,locs_bp(i)+300:length(bp_stack(1,:)));
    p = polyfit(bp',sig',1);
    x = -20:0.005:130;
    f = polyval(p,x);
%     plot(1:length(sig),sig,'r');
%     hold on;
%     plot(1:length(sig),bp,'b');
%     hold on;
    figure();
    scatter(bp,sig);
    hold on;
    ccp(i) = mean(x(round(f,2)==0));
    a = 0:1:max(bp)+10;
    f_eval = polyval(p,a);
    scatter(a,f_eval,'r','+');
    title("CCP using DCS from only the tail data");
    xlabel("ABP")
    ylabel("TCD")
end
hold off;
% figure();
% scatter(1:length(ccp),ccp,'blue');
%%
ccp = zeros(1,length(bp_stack(:,1)));
for i=1:length(bp_stack(:,1))
    p = robustfit(bp_stack(i,:)',sig_stack(i,:)');
    x = -20:0.005:130;
    f = p(1)*x+p(2);
    
    ccp(i) = mean(x(round(f,2)==0));
%     fprintf("%d",ccp(i));
end
%% robustfit implementation
p = robustfit(bp_stack(1,:)',sig_stack(1,:)');
x= -10:1:100;
y = p(1)+p(2)*x;
c = -p(1)/p(2);
%% Upsampling the signal using the linear interpolation
x = 1:1:length(dcs_1);
uf = 50;   % Upsampling factor
xq = (1:(1/uf):length(dcs_1)+((uf-1)/uf));
dcs_1up = interp1(x, dcs_1,xq,'linear');

%% Finding the minima to find the starting of the signal
sg_lp_30 = dcs_1up;
sg_lp_30 = normalize(sg_lp_30);
% % 
% minima = islocalmin(sg_lp_30,'MinSeparation',900, 'ProminenceWindow',1,'MinProminence',1 ,'FlatSelection', 'last');
% x = 1:length(minima);
% locs = x(minima)

% y = ecg1;
% x = 1:length(y);
% [pks,locs] = findpeaks(y, 'MinPeakHeight', 0.5,'MinPeakDist',700,'MinPeakProminence',0.1);
% figure();
% plot(x,sg_lp_30,x(locs),pks,'r*');

%% Plotting the data based on the minima
ini = sg_lp_30(locs(1):locs(2));
cyc =zeros(sum(minima==1)-1,length(ini)); 
cyc(1,:)= ini;
count = 2;
avg = ini;
% for i=851:1:length(minima)
%     if (minima(i)==1) && (i+849<=length(minima))
%         count = count+1;
%         plot(sg_lp_30(i:i+849))
%         hold on
%         avg = avg+sg_lp_30(i:i+849);
%          cyc(count,:) = sg_lp_30(i:i+849);
%         
%     end
% end
for i=1200:1:length(minima)
    if (minima(i)==1) && (i+length(ini)<=length(minima))
        count = count+1;
        plot(sg_lp_30(i:i+length(ini)-1))
        hold on
        avg = avg+sg_lp_30(i:i+length(ini)-1);
         cyc(count,:) = sg_lp_30(i:i+length(ini)-1);
        
    end
end
hold off
avg = avg/count;
x = (1:1:length(avg))/1000;
figure()
plot(avg)
cyc = Extract;
x = (1:1:length(cyc))/1000;
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

ensavg = [ensavg ensavg ensavg];
dcs_1a_ens = ensavg;

%Saving the variable
% writematrix(ensavg,'D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Results and Plots\output_variables\marianne\dcs_3cm_ens.csv','Delimiter','comma');





%% Test Space

% y = interp(dcs_1,50);
% x = 1:length(y);
% x_d = (1:length(dcs_1));
% [pks_u,locs_u] = findpeaks(y, 'MinPeakHeight', 0.5,'MinPeakDist',600,'MinPeakProminence',0.1);  %Determine peaks and Indices
% figure()
% [pks_d,locs_d] = findpeaks(dcs_1, 'MinPeakHeight', 0.5,'MinPeakProminence',10);
% a = locs_u./locs_d;
% plot(x_d,dcs_1)
% hold on
% plot(x_d(locs_d),pks_d, '+r')
% hold off

% Plotting the data of different modality.
% % dcs_1a_raw = normalize(interp(dcs_1,50));
% % dcs_1a_raw = dcs_1a_raw';
% % plot(ecg1);
% % hold on
% % % plot(tcd);
% % plot(normalize(dcs_1a));
% % plot(dcs_1a_raw);
% % legend('ECG','DCS 1.5cm LP Filtered','DCS 1.5cm RAW')
% % title("Signal comparision of ECG, and DCS")
% % xlabel('Samples (Time = samples/1000)');
% ini = sg_lp_30(400:1199);
% cyc =zeros(sum(minima==1)-1,800); 
% cyc(1,:)= ini;
% count = 1;
% avg = ini;
% for i=1200:1:length(minima)
%     if (minima(i)==1) && (i+799<=length(minima))
%         count = count+1;
%         plot(sg_lp_30(i:i+799))
%         hold on
%         avg = avg+sg_lp_30(i:i+799);
%          cyc(count,:) = sg_lp_30(i:i+799);
%         
%     end
% end
% hold off
% avg = avg/count;
% x = (1:1:800)/1000;
% figure()
% plot(avg)
% 
% %Plotting the ensemble average

% ensavg = mean(cyc,1);                                                   % Calculate Ensemble Average
% ci95 = 1.96*std(cyc,[],1)/sqrt(count);                             % Calculate 95% Confidence Intervals         
% figure()
% plot(x, ensavg, '-r', 'LineWidth',1)
% hold on;
% plot(x, ensavg+ci95, ':g', 'LineWidth',1.5)
% plot(x, ensavg-ci95, ':g', 'LineWidth',1.5)
% hold off
% grid
% legend('Ensemble Average', '95% Confidence Intervals')
% 
% ensavg = [ensavg ensavg];
% dcs_1a_ens = ensavg;
% 
% %Saving the variable
% writematrix(ensavg,'D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Results and Plots\output_variables\jig\dcs_3cm_ens.csv','Delimiter','comma');
% 
figure();
ini = dcs_1(locs(1):locs(2));
cyc =zeros(length(pks)-1,length(ini)); 
cyc(1,:)= ini;
x = (1:1:length(ini))/1000;
count = 0;
avg = ini;
for i=1:1:length(locs)-1
    count = count+1;
 
    sig = dcs_1(locs(i)-8:locs(i+1)-8);
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
title("Marker=ECG R peak, DCS 1cm Baseline Brachial")
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

% ensavg = [ensavg ensavg];
% dcs_1a_ens = ensavg;
% 
% %Saving the variable
% writematrix(ensavg,'D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Results and Plots\output_variables\jig\dcs_3cm_ens.csv','Delimiter','comma');

%%
% for i=1:20
%     figure()
%     plot(cyc(i,:));
%     hold on
%     plot(cyc_lp(i,:),'--');
%     hold off
% end
ecg_d= zeros(1,length(ecg1)/50)
for i=1:length(ecg_d-51)
    ecg_d(i) = median(ecg1(50i+1:50.*(i+1)));
end
plot(ecg_d)
%% Creating the exportable fil for the windkessel model
final_data = [ecg_ens; abd_ens ;tcd_ens ;dcs_1a_ens ;dcs_3a_ens]
writematrix(final_data,'D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Results and Plots\output_variables\jig\final_file.txt','Delimiter','tab');

%% Plotting the data based on the minima for the upsampled signal
ini = sg_lp_30_u(1:48);
cycu =zeros(sum(minima_u==1)-1,48); 
cycu(1,:)= ini;
count_u = 1;
avg = ini;
for i=49:1:length(minima_u)
    if (minima_u(i)==1) && (i+47<=length(minima_u))
        count_u = count_u+1;
        fprintf("%d\n",i)
        plot(sg_lp_30_u(i:i+47))
        hold on
        avg = avg+sg_lp_30_u(i:i+47);
        cycu(count_u,:) = sg_lp_30_u(i:i+47);
    end
end
hold off
x_u = (1:1:48)/60;
avg = avg/count_u;
figure()
plot(x_u,avg)

%Plotting  the ensemble average
ensavg = mean(cycu,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cycu,[],1)/sqrt(count_u);                             % Calculate 95% Confidence Intervals
figure()
plot(x_u, ensavg, '-r', 'LineWidth',1)
hold on
plot(x_u, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x_u, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
grid
legend('Ensemble Average', '95% Confidence Intervals')
xlabel('Time (s)')
ylabel('aDb *10^9')
title("DCS-1.5cm Upsampled signal")

%% Subplotting the original singal and upscaled signal

subplot(1,2,1)
ensavg = mean(cyc,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cyc,[],1)/sqrt(count);                             % Calculate 95% Confidence Intervals                                    
plot(x, ensavg, '-r', 'LineWidth',1)
hold on
plot(x, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
xlabel('Time (s)')
ylabel('aDb *10^9 ')
title("Ensemble average of original signal")
grid
legend('Ensemble Average', '95% Confidence Intervals')

subplot(1,2,2)
ensavg = mean(cycu,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cycu,[],1)/sqrt(count_u);                             % Calculate 95% Confidence Intervals  
plot(x_u, ensavg, '-r', 'LineWidth',1)
hold on
plot(x_u, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x_u, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
xlabel('Time (s)')
ylabel('aDb *10^9')
title("Ensemble average of Upsampled signal")
grid
legend('Ensemble Average', '95% Confidence Intervals')

%% Any

p = robustfit(CCP(1,:)',CCP(3,:)');
x = -1:0.005:40;
f = p(1)*x+p(2);
plot(x,f);
hold on;
scatter(CCP(1,:),CCP(3,:));

%% plotting g2 curves
for i=1:5:20
    semilogx(Data_tau,squeeze(Data(i,4,:)))
    hold on;
    scatter(Data_tau,squeeze(Data(i,4,:)))
    hold on;
end
