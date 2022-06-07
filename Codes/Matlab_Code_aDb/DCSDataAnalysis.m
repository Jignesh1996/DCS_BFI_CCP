clear all
close all
%--------------------------------------------------------------------------

% Folder name  - provide folder name for which you want to convert the data
% Folder='21.4.27-apnea without cuff';

%--------------------------------------------------------------------------

filename_d=strcat('C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\DCS_TEST\20221205','\Data.mat');
load(filename_d)
dcs = Data;



filename_nd=strcat('C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\DCS_TEST\TRNIRS\20221205','\Data.mat');
load(filename_nd)
dcs_n = Data;

dcs_nt = dcs_n(128:end,:,:);
dcs_t = dcs(121:121+size(dcs_nt,1)-1,:,:);

%%

g2(1,:,:)=squeeze(dcs_t(:,1,:)-1); %g2-1 curve generation
g2_1_temp=squeeze(dcs_t(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(dcs_t(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(dcs_t(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(dcs_t(:,4,:)-1); %g2-1 curve generation
g2_5_temp=squeeze(dcs_nt(:,1,:)-1);
g2_6_temp=squeeze(dcs_nt(:,2,:)-1);
g2_7_temp=squeeze(dcs_nt(:,3,:)-1);
g2_8_temp=squeeze(dcs_nt(:,4,:)-1);
% average g2 curve for large source detector separation
% for i=1:size(g2,2)
%     g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:)+g2_5_temp(i,:)+g2_6_temp(i,:)+g2_7_temp(i,:)+g2_8_temp(i,:))/8;
% end
for i=1:size(g2,2)
    g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
end

for i=1:size(g2,2)
    g2(2,i,:)=( g2_5_temp(i,:)+g2_6_temp(i,:)+g2_7_temp(i,:) + g2_8_temp(i,:))/4;
end
%% aDb calculations

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
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rsd,beta);
        aDb1(chan,i) = FittedParams(1);
    end
end

%% Data plotting

% time resultion - aqusition time used to aquire data

t_res=0.05; % seconds
time=t_res*(1:1:size(aDb1,2));

subplot(2,1,1)

plot(time,aDb1(1,:))
title('{\itr}_{SD}=1 cm')
% set(gca,'xticklabel',{})

subplot(2,1,2)

plot(time,aDb1(2,:))
title('{\itr}_{SD}=3 cm')
xlabel('Time (s)')

%% 
for i=150:1500:size(g2_5_temp,1)
    semilogx(Data_tau,g2_5_temp(i,:))
    title("g2 curves for DCS(hybrid) channel")
    xlabel("Data Tau")
    hold on;
end

%%
semilogx(Data_tau,g2_2_temp(250,:))
hold on;
semilogx(Data_tau,g2_6_temp(250,:))
semilogx(Data_tau,g2_7_temp(250,:))
semilogx(Data_tau,g2_8_temp(250,:))
legend("Chan 5","6","7","8")
title("g2 curves for DCS(hybrid) channel")
xlabel("Data Tau")
