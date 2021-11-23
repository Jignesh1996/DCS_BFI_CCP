clear all
close all
%--------------------------------------------------------------------------

% Folder name  - provide folder name for which you want to convert the data
Folder='21.4.27-apnea without cuff';

%--------------------------------------------------------------------------

filename=strcat('D:\Jignesh\MSc Western Uni\Research MSc\Data', '\',Folder,'\','Data.mat');
load(filename)

%%

g2(1,:,:)=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(Data(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(Data(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(Data(:,4,:)-1); %g2-1 curve generation

% average g2 curve for large source detector separation
for i=1:size(g2,2)
    g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
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

t_res=0.1; % seconds
time=t_res*(1:1:size(aDb1,2));

figure()

subplot(2,1,1)

plot(time,aDb1(1,:))
title('{\itr}_{SD}=1 cm')
set(gca,'xlim',[0 max(time)])
ylabel('BFi')

subplot(2,1,2)

plot(time,aDb1(2,:))
title('{\itr}_{SD}=3 cm')
xlabel('Time (s)')
set(gca,'xlim',[0 max(time)])
ylabel('BFi')

%% Average Data plotting

% time resultion - aqusition time used to aquire data
t_res=0.1; % seconds
t_avg=5; % window used for averaging in seconds

t_avg_pt=t_avg/t_res; % window used for averaging in points
time=t_res*(1:1:size(aDb1,2));

j=1;
for i=1:t_avg_pt:size(aDb1,2)
    aDb1_avg(1,j)=mean(aDb1(1,i:i+t_avg_pt-1));
    aDb1_avg(2,j)=mean(aDb1(2,i:i+t_avg_pt-1));
    time_avg(1,j)=mean(time(1,i+t_avg_pt/2-1));
    j=j+1;
end

%%
%Ploting avarege time courses 

figure()

subplot(2,1,1)

plot(time_avg,aDb1_avg(1,:))
title('{\itr}_{SD}=1 cm')
set(gca,'xlim',[0 max(time)])
ylabel('BFi')

subplot(2,1,2)

plot(time_avg,aDb1_avg(2,:))
title('{\itr}_{SD}=3 cm')
xlabel('Time (s)')
set(gca,'xlim',[0 max(time)])
ylabel('BFi')

%%
%Ploting avarege time courses as a % changes 

base_period=120; % in seconds
base_period_pt=base_period/(time_avg(1,2)-time_avg(1,1)); %baseline period in points


delta_aDb1_avg(1,:)=100*(aDb1_avg(1,:)/mean(aDb1_avg(1,1:base_period_pt))-1); 
delta_aDb1_avg(2,:)=100*(aDb1_avg(2,:)/mean(aDb1_avg(2,1:base_period_pt))-1); 
%%
figure()

subplot(2,1,1)

plot(time_avg,delta_aDb1_avg(1,:));
title('{\itr}_{SD}=1 cm')
set(gca,'xlim',[0 max(time)])
ylabel('\DeltaBFi (%)')
set(gca,'ylim',[-10 100])
line([0 max(time)], [0 0],'LineStyle',':','Color','k')
line([base_period base_period], [-100 100],'LineStyle','--','Color','r')
line([base_period+300 base_period+300], [-100 100],'LineStyle','--','Color','r')

subplot(2,1,2)

plot(time_avg,delta_aDb1_avg(2,:))
title('{\itr}_{SD}=3 cm')
xlabel('Time (s)')
set(gca,'xlim',[0 max(time)])
ylabel('\DeltaBFi (%)')
set(gca,'ylim',[-10 100])
line([0 max(time)], [0 0],'LineStyle',':','Color','k')
line([base_period base_period], [-100 100],'LineStyle','--','Color','r')
line([base_period+300 base_period+300], [-100 100],'LineStyle','--','Color','r')

%% Find peaks function and ensemble averaging - use this one May 2021
y = aDb1(2,3900:4200);                                                          %set data you want to analyze
x = time(:,3900:4200);

[pks,locs] = findpeaks(y, 'MinPeakHeight', 0.000000003,'MinPeakDist',5);  %Determine peaks and Indices
figure()
plot(x,y)
hold on
plot(x(locs),pks, '+r')
hold off
grid

clear yc2
clear xc2

for k1 = 1:numel(locs)-1
    yc2{k1} = y(locs(k1)-1:locs(k1+1)-1);                            % Define Frames
    xc2{k1} = x(locs(k1)-1:locs(k1+1)-1);
end

figure()
hold all 
for k1 = 1:numel(yc2)
   plot(xc2{k1}-xc2{k1}(1),yc2{k1})
end
hold off
grid

% — CALCULATE & PLOT ENSEMBLE AVERAGE —                                                                    
minlen = min(cellfun(@numel, yc2));                                     % Minimum Length Of Records
ens = zeros(minlen, numel(yc2));                                        % Preallocate
for k1 = 1:numel(yc2)
    ens(:,k1) = yc2{k1}(1:minlen);                                      % Trim  To Shortest Length
end
ensavg = mean(ens,2);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(ens,[],2)/sqrt(numel(yc2));                             % Calculate 95% Confidence Intervals
eatv = mean(diff(x))*(0:minlen-1);                                     % Time Vector For Ensemble AVerage
figure()
plot(eatv, ensavg, '-r', 'LineWidth',1)
hold on
plot(eatv, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(eatv, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
grid
legend('Ensemble Average', '95% Confidence Intervals')
%% grouping all diastolic values together 
[pks_inv,locs_inv] = findpeaks(-y, 'MinPeakDist',5);  %Determine peaks and Indices
pks_inv = -pks_inv
%figure()
%plot(x,-y)
%hold on
%plot(x(locs_inv),pks_inv, '+r')
%hold off
%grid
for k1 = 1:numel(locs_inv)-1
    yc4{k1} = y(locs_inv(k1)-5:locs_inv(k1+1)-5);                            % Define Frames
    xc4{k1} = x(locs_inv(k1)-5:locs_inv(k1+1)-5);
end

figure()
hold all 
for k1 = 1:numel(yc4)
   plot(xc4{k1}-xc4{k1}(1),yc4{k1})
end
hold off
grid

% — CALCULATE & PLOT ENSEMBLE AVERAGE —                                                                    
minlen = min(cellfun(@numel, yc4));                                     % Minimum Length Of Records
ens = zeros(minlen, numel(yc4));                                        % Preallocate
for k1 = 1:numel(yc4)
    ens(:,k1) = yc4{k1}(1:minlen);                                      % Trim  To Shortest Length
end
ensavg = mean(ens,2);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(ens,[],2)/sqrt(numel(yc4));                             % Calculate 95% Confidence Intervals
eatv = mean(diff(x))*(0:minlen-1);                                     % Time Vector For Ensemble AVerage
figure()
plot(eatv, ensavg, '-r', 'LineWidth',1)
hold on
plot(eatv, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(eatv, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
grid
legend('Ensemble Average', '95% Confidence Intervals')

%% Plotting Systolic and Diastolic Values, and combining into one table (1st row is sys, 2nd row is dias)
[minpeak, P] = islocalmin(y, 'MinSeparation',6);
valley = P(minpeak)

figure()
plot(x,y)
hold on
plot(x(locs),pks, '+r')
plot(x,y,x(minpeak),y(minpeak),'b+');
hold off
grid
%CHECK IF SYS is first marker on graph! Prepend the array to make sure
%diastolic value FOLLOWS systolic value 
%pks = [zeros(1,1),pks]; 
pksl = length(pks);
pks_invl = length(pks_inv);
 
pksl
pks_invl
%if pks_inv is greater, padd with the following command: 
pks_padd = [pks,zeros(1,length(pks_inv)-length(pks))];
PI_grid = [pks_padd;pks_inv];

% if pks is greater, padd with the following command: 
%pks_inv_padd=[pks_inv,zeros(1,length(pks)-length(pks_inv))];
%PI_grid = [pks;pks_inv_padd];

PI_grid
PI_grid(PI_grid==0)=NaN
%PI = pks-pks_inv_padd
PI = pks_padd-pks_inv

mean(PI)


