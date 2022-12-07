close all;
clear all;

filepath = ""  % Path to the file

load(filepath);
%% Loading the file
% load("D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\Pressure_modulation\Aleena\Aleena_Probe_pressure_1\Data.mat");
aDb1 = standalone_dcs(Data,Data_tau);
tau = Data_tau;
%%

g2(1,:,:)=squeeze(Data(:,1,:)); %g2-1 curve generation
g2_1_temp=squeeze(Data(:,1,:)); %g2-1 curve generation
g2_2_temp=squeeze(Data(:,2,:)); %g2-1 curve generation
g2_3_temp=squeeze(Data(:,3,:)); %g2-1 curve generation
g2_4_temp=squeeze(Data(:,4,:)); %g2-1 curve generation

for i=1:size(g2,2)
    g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
end


% Define the span for below conditions.
bsl = 1:600;  % Baseline span
pres = 601:900; %Pressure Modulation span

%% Plotting the mean g2 curves for pressure and baseline for comparision

semilogx(tau,(mean(squeeze(g2(1,pres,:)))-1),tau,(mean(squeeze(g2(1,bsl,:)))-1))
legend("Pressure Modulation","Baseline")
title("g2 comparision");
xlabel("Tau");

%% Making sure that all the values are more than 1 in the g2 curve
g2(g2<=1)=1;
%% Finding the sensitivity factor ratio

g2_temp=round(g2,4);

DeltaOD_sp=-log((mean(squeeze(g2_temp(1,pres,:)))-1)./(mean(squeeze(g2_temp(1,bsl,:)))-1));
DeltaOD_lp=-log((mean(squeeze(g2_temp(2,pres,:)))-1)./(mean(squeeze(g2_temp(2,bsl,:)))-1));
DeltaOD_s = -log(((squeeze(g2_temp(1,:,:)))-1)./(mean(squeeze(g2_temp(1,bsl,:)))-1));
DeltaOD_l =  -log(((squeeze(g2_temp(2,:,:)))-1)./(mean(squeeze(g2_temp(2,bsl,:)))-1));

% DeltaFec(1,:)=DeltaOD_sp/(0.3*10^-9);
% DeltaFec(2,:)=DeltaOD_lp/(0.3*10^-9);

DeltaFec(1,:)=(DeltaOD_sp);
DeltaFec(2,:)=(DeltaOD_lp);
% 
% DeltaOD_s = -log((g2_sig(1,:)-1)/(g2(1,:)-1));
% DeltaOD_l = -log((g2_sig(2,:)-1)/(g2(3,:)-1));

ratio = DeltaFec(2,:)./DeltaFec(1,:);
semilogx(tau,DeltaFec(2,:)./DeltaFec(1,:));
title("Sensitivity Factor Ratio");
xlabel("Tau");
ylabel("S.F. Ratio");

index=ratio==Inf;


sens_fact = nanmean(ratio(index==0))

figure()

semilogx(tau,DeltaFec(1,:),'r')
hold on
semilogx(tau,DeltaFec(2,:),'b')
title("DeltaFec");
xlabel("Tau");

%% Calculating the DeltaFc
DeltaFc = (DeltaOD_l) - (((DeltaOD_lp)./(DeltaOD_sp)).*(DeltaOD_s)); % This  equation also has a multiplicative factor.
