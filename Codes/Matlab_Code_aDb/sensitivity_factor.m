% This code is for calculating the sensitivity factor for the DCS data from
% the hybrid system. make the appropriate changes for the standalone system
% dcs data. 


close all;
clear all;

% load("E:\PARKWOOD_STUDY\DCS\MCI003_PRESSURE\Data.mat");
load("D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Matlab_Code_aDb\2 layer model\2 layer model\001_CUFF_g2.mat");
tau = csvread("D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Matlab_Code_aDb\2 layer model\data_tau.csv");
Data_tau = tau;


%Since this data loaded is in the g2 format, following code convert it into
%the Data format for future use.
Data(:,1,:) = g2(1,:,:)+1;
Data(:,2,:) = g2(2,:,:)+1;
Data(:,3,:) = g2(3,:,:)+1;
Data(:,4,:) = g2(4,:,:)+1;

aDb1 = hybrid_dcs(Data,Data_tau); %hybrid_dcs is a function for calculating the BFi values.
figure();
plot((1:1:length(aDb1))/20,aDb1');
xlabel("Time(s)")
ylabel("aDb")

%% fitting for aDb

rho = [1 1.5 2 2.5]; %source detector separations in cm 
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 8; %cm^-1 baseline reduced scattering coefficient

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
for i=1:10:100
    Channel=2;
    Curve_no=i;
    rho = [1 1.5 2 2.5];
    
    beta=g2(Channel, Curve_no,1);
    aDb_fit=aDb1(Channel,Curve_no);
    
    g2_fit=gen_DCS_fit(Data_tau,mua,mus,rho(Channel),beta,aDb_fit);
    
    semilogx(Data_tau,squeeze(g2(1,1,:)),'k')
    hold on
    semilogx(Data_tau,g2_fit,'r')
    hold on
end
%%
% Data = g2_n;

bsl = 1:550;
pres = 650:1200;


% if exist("g2")
%     clear g2
% end

close all;
g2(1,:,:)=squeeze(Data(:,1,:)); %g2-1 curve generation
g2(2,:,:)=squeeze(Data(:,4,:)); %g2-1 curve generation

g2_p(1,:,:)=squeeze(Data(:,1,:))-1; %g2-1 curve generation
g2_p(2,:,:)=squeeze(Data(:,4,:))-1; %g2-1 curve generation


l1 = size(g2,2);

for i=1:l1
    g2_new = g2;
    g2_new(:,i,:) = g2(:,i,:)-mean(g2(:,i,end-5:end),3);
    g2_new(:,i,:) = (g2_new(:,i,:)./max(g2_new(:,i,:)))/2+1;
end
for i=1:l1
%     g_b=squeeze(g2(1,i,:))';
%     g_b=g_b-mean(g_b(1,end-10:end));
%     g_b=(g_b/max(g_b))/2+1;
%     g_b(1) = round(g_b(1),1);
%     g2(1,i,:) = g_b;

    bsl_art=mean(g2_p(1,bsl,end));

    g_b=squeeze(g2(1,i,:))';
    g_b=g_b-g_b(1,end);
    g_b=(g_b/max(g_b));
%     g_b(1) = round(g_b(1),1);
    g2(1,i,:) = g_b*0.5+1;

 
%     g_p=squeeze(g2(2,i,:))';
%     g_p=g_p-mean(g_p(1,end-10:end));
%     g_p=(g_p/max(g_p))/2+1;
%     g_p(1) = round(g_p(1),1);
%     g2(2,i,:) = g_p;
    bsl_art=mean(g2_p(2,bsl,end));

    g_b=squeeze(g2(2,i,:))';
    g_b=g_b-g_b(1,end);
    g_b=(g_b/g_b(1,1));
%     g_b(1) = round(g_b(1),1);
    g2(2,i,:) = g_b*0.5+1;
end

for i=1:10:l1
    semilogx(Data_tau,squeeze(g2(2,i,:)));
    hold on;
end

g2_temp=round(g2,4);
% g2_temp(g2_temp<=1)=1;

%%
close all;
bsl = 1:550;
pres = 650:1200;
    
g2_sf = zeros(2,2,50);

figure();
semilogx(Data_tau,(mean(squeeze(g2_temp(1,pres,:)))));
hold on;
semilogx(Data_tau,(mean(squeeze(g2_temp(1,bsl,:)))));
semilogx(Data_tau,(mean(squeeze(g2_temp(2,pres,:)))));
semilogx(Data_tau,(mean(squeeze(g2_temp(2,bsl,:)))));
legend("Pressure1","Baseline1","Pressure2","Baseline2")

hold off;

figure();

DeltaOD_sp=real(-log((mean(squeeze(g2_temp(1,pres,:)))-1)./(mean(squeeze(g2_temp(1,bsl,:)))-1)));
DeltaOD_lp=real(-log((mean(squeeze(g2_temp(2,pres,:)))-1)./(mean(squeeze(g2_temp(2,bsl,:)))-1)));

g2_sf(1,1,:) = (mean(squeeze(g2_temp(1,bsl,:))));
g2_sf(1,2,:) = (mean(squeeze(g2_temp(1,pres,:))));
g2_sf(2,1,:) = (mean(squeeze(g2_temp(2,bsl,:))));
g2_sf(2,2,:) = (mean(squeeze(g2_temp(2,pres,:))));

% DeltaFec(1,:)=DeltaOD_sp/(0.3*10^-9);
% DeltaFec(2,:)=DeltaOD_lp/(0.3*10^-9);

DeltaFec(1,:)=(DeltaOD_sp);
DeltaFec(2,:)=(DeltaOD_lp);
% 
% DeltaOD_s = -log((g2_sig(1,:)-1)/(g2(1,:)-1));
% DeltaOD_l = -log((g2_sig(2,:)-1)/(g2(3,:)-1));

ratio = DeltaFec(2,:)./DeltaFec(1,:);
semilogx(tau,DeltaFec(2,:)./DeltaFec(1,:));


ratio_d = ratio(3:end-25);

index=ratio_d==Inf;

sens_fact = nanmean(ratio_d(index==0))

% sens_fact = mean(ratio_d)

figure()

semilogx(tau,DeltaFec(1,:),'r')
hold on
semilogx(tau,DeltaFec(2,:),'b')
legend("Short Rsd","Long Rsd")


%% Calculating Delta Fc
close all

DeltaOD_s = -log(((squeeze(g2_temp(1,:,:)))-1)./(mean(squeeze(g2_temp(1,bsl,:)))-1));
DeltaOD_l =  -log(((squeeze(g2_temp(2,:,:)))-1)./(mean(squeeze(g2_temp(2,bsl,:)))-1));

% DeltaFc = ((DeltaOD_l) - ((DeltaOD_lp./DeltaOD_sp).*(DeltaOD_s)));
% DeltaFc = ((DeltaOD_l) - ((sens_fact).*(DeltaOD_s)));
DeltaFc = ((DeltaOD_l) - ((0.35).*(DeltaOD_s)));
DeltaFc(isinf(DeltaFc)) = nan;
t = (1:1:size(DeltaFc,1))/20;
mean_dfc = nanmean(DeltaFc(:,2:25),2);
mean_dfc_sm = smooth(mean_dfc,0.015,'rlowess');
per_change_dfc = mean_dfc*100;
change_dfc_flip = (mean_dfc_sm-mean(mean_dfc_sm(1:1800))/mean(mean_dfc_sm(1:1800)));
per_change_adb(1,:) = ((aDb1(1,:)-mean(aDb1(1,1:600),2))./mean(aDb1(1,1:600),2))*100;
per_change_adb(2,:) = ((aDb1(2,:)-mean(aDb1(2,1:600),2))./mean(aDb1(2,1:600),2))*100;
per_c_adb_sm(1,:) = smooth(per_change_adb(1,:),0.015);
per_c_adb_sm(2,:) = smooth(per_change_adb(2,:),0.015);
plot(t,mean_dfc); xlabel("Time (s)"); ylabel("\DeltaF_c");
% hold on; 
figure();
plot(t,mean_dfc_sm);
ylim([-1 1])
xlabel("Time (s)")
ylabel("Mean DeltaFc")
title("Mean DeltaFc across tau(2:25)")
figure(); 
plot(t,per_c_adb_sm(1,:));
hold on;
plot(t,per_c_adb_sm(2,:));

xlabel("Time (s)")
hold on;
ylabel("% change")
title("Pecrent Change in aDb wrt Baseline")
legend("aDb 1cm","aDb 2.5cm")
figure();
plot(t, change_dfc_flip);

% csvwrite("2 layer model\009_DFc.csv",mean_dfc_sm)