% This code analyzes the pulsatility of the muscles in the arm and the
% thumb muscle (Thinner muscle) in the palm. This uses the ABP signal to
% extract the pulsatility from the DCS signal.

close all;
clear all;
    
load("D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\new pressure cuff\leena test2\Data.mat")
load("D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\new pressure cuff\Leena_tourniquet.mat")
% load("D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\Donya\test2\daniel\thenar\Data.mat")
% load("D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\Donya\test\jignesh_test.mat")


g2(1,:,:)=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2_1_temp=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(Data(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(Data(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(Data(:,4,:)-1); %g2-1 curve generation

for i=1:size(g2,2)
    g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
end

rho = [1 2.5]; %source detector separations in cm 
mua = 0.15; %cm^-1 baseline absorption coefficient
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
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rsd,beta);
        aDb1(chan,i) = FittedParams(1);
    end
end

%%
tsig = (1:1:length(aDb1))/20;
plot(tsig,aDb1','LineWidth',1.5)
xlabel("Time",'FontSize',12,'FontWeight','bold');
ylabel("\alphaD_b", 'FontSize',12,'FontWeight','bold')
legend("1 cm", "2.5 cm")
title("arm", 'FontSize',14,'FontWeight','bold')
%% Assessing the fitting of the data
close all
Channel=2;
for i = 1:100
Curve_no=i;
rho = [1 2.5];

beta=g2(Channel, Curve_no,1);
aDb_fit=aDb1(Channel,Curve_no);

g2_fit=gen_DCS_fit(Data_tau,mua,mus,rho(Channel),beta,aDb_fit);

semilogx(Data_tau,squeeze(g2(Channel,Curve_no,:)),'k')
hold on
semilogx(Data_tau,g2_fit,'r')
end
%% Extracting the ABP signal of interest
strp =60; % Start point of the ABP signal
fs = 200; %sampling frequency of the singal
td = 60; % duration of the singal to be extracted.

abp = data(strp*fs:(strp+td)*fs-1);

figure();
t1 = (1:1:length(abp))/200;
plot(t1,normalize(abp));
hold on;
t2 = (1:1:length(aDb1))/20;
plot(t2,normalize(aDb1'))

%%
close all;
for i=1:1
    figure()
    [ens_avg_sig,pind] = abp_ensemble_avg(abp,aDb1(1,1200:2400),200,-500,0);
end
%% Shifting the signal to make sure the sync with the ABP sigal.
abps = circshift(abp,-470);

figure();
t1 = (1:1:length(abps))/200;
plot(t1,normalize(abps));
hold on;
t2 = (1:1:length(aDb1))/20;
plot(t2,normalize(aDb1'))
%% Extracting the pulsatility using the ABP as a marker

[h_pks,l_pks] = findpeaks(normalize(abps),"MinPeakHeight",1.5,"MinPeakDistance",150);
hd_pks = ceil(h_pks./10);
ld_pks = ceil(l_pks./10);
close all;
% if exist("g2")
%     clear g2;
% end
mua = 0.15; %cm^-1 baseline absorption coefficient
mus = 10; 

    
% 
g2_backup = Data;
ecg = abps;
% g2_backup = g2;
g2_n = g2_backup;

%  strt_time =  [1,32,62];
%  stp_time =  [28,58,90]; 

   
    clear g2;
   
    g2 = g2_n;
    ecg = abps;
    
    
    [h_pks,l_pks] = findpeaks(normalize(ecg),"MinPeakHeight",1.5,"MinPeakDistance",150);
    
    fig1=figure('units','centimeters', 'Position',[2 2 35 13]); %18 width 15 heigh
    hold on
    plot(normalize(ecg),'r')
    plot(l_pks, h_pks,'*k')
    hd_pks = ceil(h_pks./10);
    ld_pks = ceil(l_pks./10);
    % g2_avg = zeros(size(g2,1),size(g2,2),size(g2,3));
    avg_window_width = 5;
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

    g2_n = g2;





adb_avg = standalone_dcs(g2_n,Data_tau);
% SNR(1,k) = snr(adb_avg(2,:),20);


dcs_1lp = (lpf(adb_avg(1,:),7,20));
dcs_25lp =(lpf(adb_avg(2,:),7,20));
% snr(dcs_25lp(8400:10800),20)

figure();
plot(adb_avg(2,:),'b',"LineWidth",1.5);
hold on; 
plot(aDb1(2,:),'r');
legend("Ensemle Temporal Averaged","Raw")



%% Applying the filter in the aDb signal
dcs_1lp = (lpf(aDb1(1,:),4,20));
dcs_25lp =(lpf(aDb1(2,:),4,20));
% dcs_1lp = aDb1(1,:);
% dcs_25lp =aDb1(2,:);

%% Extracting the pulsatility
close all;
[ens_avg_sig,pind] = abp_ensemble_avg(abps,adb_avg(2,:),200,0,0);


