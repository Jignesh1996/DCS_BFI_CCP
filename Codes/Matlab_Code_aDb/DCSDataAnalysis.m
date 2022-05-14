clear all
close all
%--------------------------------------------------------------------------

% Folder name  - provide folder name for which you want to convert the data
Folder='21.4.27-apnea without cuff';

%--------------------------------------------------------------------------

filename=strcat('D:\Jignesh\MSc Western Uni\Research MSc\Data', '\',Folder,'\','Data.mat');
load(filename)

%%

g2(1,:,:)=squeeze(dcs_t(:,1,:)-1); %g2-1 curve generation
g2_1_temp=squeeze(dcs_t(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(dcs_t(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(dcs_t(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(dcs_t(:,4,:)-1); %g2-1 curve generation
g2_5_temp=squeeze(nirs_t(:,1,:)-1); %g2-1 curve generation
g2_6_temp=squeeze(nirs_t(:,2,:)-1); %g2-1 curve generation
g2_7_temp=squeeze(nirs_t(:,3,:)-1); %g2-1 curve generation
g2_8_temp=squeeze(nirs_t(:,4,:)-1); %g2-1 curve generation

% average g2 curve for large source detector separation
for i=1:size(g2,2)
    g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:)+g2_5_temp(i,:)+g2_6_temp(i,:)+g2_7_temp(i,:)+g2_8_temp(i,:))/8;
end
% for i=1:size(g2,2)
%     g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
% end
%% aDb calculations

rho = [1.7 2.7]; %source detector separations in cm 
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
