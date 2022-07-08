clear all
close all
%--------------------------------------------------------------------------

% Folder name  - provide folder name for which you want to convert the data
Folder='21.4.27-apnea without cuff';

%--------------------------------------------------------------------------

filename=strcat('D:\Jignesh\MSc Western Uni\Research MSc\Data\DCS\TR\20211116\50ms\','Data.mat');
load(filename)

%%

g2(1,:,:)=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2(2,:,:)=squeeze(Data(:,2,:)-1); %g2-1 curve generation
g2(3,:,:)=squeeze(Data(:,3,:)-1); %g2-1 curve generation
g2(4,:,:)=squeeze(Data(:,4,:)-1); %g2-1 curve generation

%% aDb calculations

rho = [1 1.5 2 2.5]; %source detector separations in cm 
mua = 0.1; %cm^-1 baseline absorption coefficient
mus = 10; %cm^-1 baseline reduced scattering coefficient
g2 = g2_new;
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

t_res=0.25; % seconds
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

plot(time,aDb1(1,:))
title('{\itr}_{SD}=2 cm')
% set(gca,'xticklabel',{})

subplot(2,2,4)

plot(time,aDb1(2,:))
title('{\itr}_{SD}=2.5 cm')
xlabel('Time (s)')