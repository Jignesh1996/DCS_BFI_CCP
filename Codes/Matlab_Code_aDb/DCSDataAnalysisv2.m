% clear all
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

%%Find peaks function
y = aDb1(1,:);
x = time;

[pks,locs] = findpeaks(y, 'MinPeakHeight', 0.000000003,'MinPeakDist',0.5,'MinPeakProminence',0.0000000015);  %Determine peaks and Indices
figure()
plot(x,y)
hold on
plot(x(locs),pks, '+r')
hold off
grid

for k1 = 1:numel(locs)-0.1
    yc{k1} = y(locs(k1)-1:locs(k1+1)-1);                            % Define MUAP Frames
    xc{k1} = x(locs(k1)-1:locs(k1+1)-1);
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
legend('Ensemble Average', '95% Confidence Intervals')

%% Upsampled y and finding its peaks
y = aDb1(1,:);
%Upsampling the y to increase the resolution
y_up = interp(y, 50 );
x_up = interp(x,50);
x = time;


% [pks,locs] = findpeaks(y_up, 'MinPeakHeight', 0.000000003,'MinPeakDist',0.3,'MinPeakProminence',0.0000000008);  %Determine peaks and Indices
[pks,locs] = findpeaks(y_up,'MinPeakDist',0.8,'MinPeakProminence',0.000000001);  %Determine peaks and Indices

figure()

plot(x_up,y_up)
hold on
plot(x_up(locs),pks, '+r')
hold off
grid

for k1 = 1:numel(locs)-0.1
    yc{k1} = y_up(locs(k1)-50:locs(k1+1)-50);                            % Define MUAP Frames
    xc{k1} = x_up(locs(k1)-50:locs(k1+1)-50);
end

figure()
hold all 
for k1 = 1:numel(yc)
   plot(xc{k1}-xc{k1}(1),yc{k1})
end
hold off
grid