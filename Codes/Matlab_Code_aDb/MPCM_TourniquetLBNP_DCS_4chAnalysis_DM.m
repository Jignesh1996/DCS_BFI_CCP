
clear all
close all
%--------------------------------------------------------------------------

% % Folder name  - provide folder name for which you want to convert the data
Folder='20220309 - 5';

%--------------------------------------------------------------------------

% filename=strcat(pwd, '\DCS\',Folder,'\','data.mat');
filename=strcat(pwd, '\DCS\',Folder,'\','data.mat');
load(filename);

%
temp_res=0.05; %in seconds

%% TR was the marker for time-resolved for NEMO system... just 
TR=1;

for i=1:size(Data,1)
    Data_cut(i,:,:)=Data(i,:,:);
    Data_time(1,i)=i*temp_res;
end

%%
g2_1_temp=squeeze(Data_cut(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(Data_cut(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(Data_cut(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(Data_cut(:,4,:)-1); %g2-1 curve generation

% average g2 curve for 4 channels
if TR==0
    g2(1,:,:)=(g2_4_temp);
    g2(2,:,:)=(g2_1_temp+g2_2_temp)/2;
    g2(2,:,:)=(g2_3_temp);
else
    g2(1,:,:)=(g2_1_temp);
    g2(2,:,:)=(g2_2_temp);
    g2(3,:,:)=(g2_3_temp);
    g2(4,:,:)=(g2_4_temp);
end
 
if TR==0 % for stand alone
    semilogx(Data_tau,squeeze((g2(1,1,:))),'r')
    hold on
    semilogx(Data_tau,squeeze((g2(2,1,:))),'b')
    ylabel('g_2')
    legend('1 cm','2.5 cm')
else % for TR
    semilogx(Data_tau,squeeze((g2(1,1,:))),'r')
    hold on
    semilogx(Data_tau,squeeze((g2(2,1,:))),'g')
    semilogx(Data_tau,squeeze((g2(3,1,:))),'b')
    semilogx(Data_tau,squeeze((g2(4,1,:))),'magenta')
    ylabel('g_2')
    legend('1 cm','1.5 cm','2 cm','2.5 cm')
end
%%
%  j=1;
%  for i=1:moving_window:size(g2,1)
%      g2_new(j,:)=mean(g2(i:i+moving_window-1,:));
%      j=j+1;
%  end
%      

%% aDb calculations
%stand alone

if TR==0
    rho = [1 2.5]; %source detector separations in cm 
else 
    rho = [1 1.5 2 2.5]; %source detector separations in cm 
end

mua = 0.1; %cm^-1 baseline absorption coefficient
mus = 10; %cm^-1 baseline reduced scattering coefficient

tau_values=Data_tau;
for r=1:size(rho,2)
    for i=1:size(g2,2)
        g2_temp(i,:)=squeeze(g2(r,i,:))';
        LB = [0.1e-9];
        UB = [10e-8];
        Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
        beta= squeeze((g2(r,i,1:1))); %0.1568;
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rho(1,r),beta);
        aDb1(r,i) = FittedParams(1);
    end
end


%%

figure()

if TR==0
    windowsize=10;
    wagi=ones(1,windowsize)/windowsize;
    aDb1_plot(1,:)=filtfilt(wagi,1,aDb1(1,:));
    aDb1_plot(2,:)=filtfilt(wagi,1,aDb1(2,:));

    subplot(2,1,1)
    hold on

    aDb1_avg=(mean(aDb1(1,:)));
    aDb2_avg=(mean(aDb1(2,:)));
    aDb1_std=(std(aDb1(1,:)));
    aDb2_std=(std(aDb1(2,:)))

    lineProps.col{1} = 'k';
    lineProps.edgestyle = '';
    lineProps.style = '--';
    % rectangle('Position',[30,-79,30,180],'FaceColor',[0.9 .9 .9],'EdgeColor','none',...
    %     'LineWidth',3)
    % plot(Data_time,aDb1(1,:)/mean(aDb1(1,:)),'r','LineWidth',1,'MarkerSize',2)
    plot(Data_time,aDb1_avg,aDb1_std,lineProps,1)
    title('rCBF')
    legend('1 cm')

    subplot(2,1,2)
    hold on
    plot(Data_time,aDb1(2,:)/mean(aDb1(2,:)),'k','LineWidth',1,'MarkerSize',2)
    plot(Data_time,aDb2_avg,aDb2_std,lineProps,1)
    legend('2.5cm')
    set(gca,'ylim',[0.5 1.5])
    xlabel('Time (s)')

else
    
    windowsize=10;
    wagi=ones(1,windowsize)/windowsize;
    aDb1_plot(1,:)=filtfilt(wagi,1,aDb1(1,:));
    aDb1_plot(2,:)=filtfilt(wagi,1,aDb1(2,:));
    aDb1_plot(3,:)=filtfilt(wagi,1,aDb1(3,:));
    aDb1_plot(4,:)=filtfilt(wagi,1,aDb1(4,:));


    subplot(4,1,1)
    hold on

    plot(Data_time/60,aDb1_plot(1,:),'r','LineWidth',1,'MarkerSize',2)
    title('rCBF (aDb)')
    legend('1cm')

    subplot(4,1,2)
    hold on 

    plot(Data_time/60,aDb1_plot(2,:),'m','LineWidth',1,'MarkerSize',2)
    %title('rCBF')
    legend('1.5cm')

    subplot(4,1,3)
    hold on

    plot(Data_time/60,aDb1_plot(3,:),'b','LineWidth',1,'MarkerSize',2)
    %title('rCBF')
    legend('2cm')
    
    subplot(4,1,4)
    hold on
    plot(Data_time/60,aDb1_plot(4,:),'k','LineWidth',1,'MarkerSize',2)
    legend('2.5 cm')
    hold on
    % set(gca,'ylim',[0.5 1.5])
    xlabel('Time (min)')

end
    
%% Data plotting w no filter

% time resultion - aqusition time used to aquire data

t_res=0.05; % seconds
time=t_res*(1:1:size(aDb1,2));

figure()

if TR==0

    subplot(2,1,1)

    plot(time/60,aDb1(1,:))
    title('{\itr}_{SD}=1 cm')
    set(gca,'xlim',[0 max(time/60)])
    ylabel('BFi')

    subplot(4,1,2)

    plot(time/60,aDb1(2,:))
    title('{\itr}_{SD}=2.5 cm')
    set(gca,'xlim',[0 max(time/60)])    
    ylabel('BFi')

else

    subplot(4,1,1)

    plot(time/60,aDb1(1,:))
    title('{\itr}_{SD}=1 cm')
    set(gca,'xlim',[0 max(time/60)])
    ylabel('BFi')

    subplot(4,1,2)

    plot(time/60,aDb1(2,:))
    title('{\itr}_{SD}=1.5 cm')
    set(gca,'xlim',[0 max(time/60)])
    ylabel('BFi')

    subplot(4,1,3)

    plot(time/60,aDb1(3,:))
    title('{\itr}_{SD}=2 cm')
    set(gca,'xlim',[0 max(time/60)])
    ylabel('BFi')

    subplot(4,1,4)

    plot(time/60,aDb1(4,:))
    title('{\itr}_{SD}=2.5 cm')
    set(gca,'xlim',[0 max(time/60)])
    ylabel('BFi')
    xlabel('Time (min)')

end
%% Average Data plotting (no filter)

clear aDb1_avg time_avg

% time resultion - aqusition time used to aquire data
t_res=0.05; % seconds
t_avg=1; % window used for averaging in seconds

t_avg_pt=t_avg/t_res; % window used for averaging in points
time=t_res*(1:1:size(aDb1,2));

j=1;
for i=1:t_avg_pt:size(aDb1,2)
    aDb1_avg(1,j)=mean(aDb1(1,i:i+t_avg_pt-1));
    aDb1_avg(2,j)=mean(aDb1(2,i:i+t_avg_pt-1));
    aDb1_avg(3,j)=mean(aDb1(3,i:i+t_avg_pt-1));
    aDb1_avg(4,j)=mean(aDb1(4,i:i+t_avg_pt-1));
    time_avg(1,j)=(Data_time(1,i+t_avg_pt-1));
    j=j+1;
end

figure()

txt4='LBNP -30 mmHg\rightarrow';
txt5='LBNP -80 mmHg\rightarrow';

subplot(4,1,1)
hold on
x_pos=[2, 2.5, 3, 3.5]; %task strat time in minutes

for i=1:size(x_pos,2)

rectangle('Position',[x_pos(i),0,.25,0.9*max(aDb1_avg(1,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
    'LineWidth',3)
end

plot(time_avg/60,aDb1_avg(1,:),'r','LineWidth',1)
title('{\itr}_{SD}=1 cm')
set(gca,'xlim',[0 max(time/60)])
ylabel('BFi')

text(1.4,0.00000001,txt5,'FontSize',6)

Range(1,1)=0.9*min(aDb1_avg(4,:));
Range(1,2)=1.1*max(aDb1_avg(4,:));

subplot(4,1,2)
hold on

% temp_axaDb=delta_aDb1_avg(1,:)

for i=1:size(x_pos,2)

rectangle('Position',[x_pos(i),0,.25,0.9*max(aDb1_avg(2,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
    'LineWidth',3)
end

plot(time_avg/60,aDb1_avg(2,:),'m','LineWidth',1)
title('{\itr}_{SD}=2.5 cm')
set(gca,'xlim',[0 max(time/60)])
ylabel('BFi')

subplot(4,1,3)
hold on

% temp_axaDb=delta_aDb1_avg(1,:)

for i=1:size(x_pos,2)

rectangle('Position',[x_pos(i),0,.25,0.9*max(aDb1_avg(3,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
    'LineWidth',3)
end

plot(time_avg/60,aDb1_avg(3,:),'b','LineWidth',1)
title('{\itr}_{SD}=2 cm')
set(gca,'xlim',[0 max(time/60)])
ylabel('BFi')

subplot(4,1,4)
hold on

% temp_axaDb=delta_aDb1_avg(1,:)

for i=1:size(x_pos,2)

rectangle('Position',[x_pos(i),0,.25,0.9*max(aDb1_avg(4,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
    'LineWidth',3)
end

plot(time_avg/60,aDb1_avg(4,:),'k','LineWidth',1)
title('{\itr}_{SD}=2.5 cm')
xlabel('Time (s)')
set(gca,'xlim',[0 max(time/60)])
ylabel('BFi')

%% Percent Change (BFi)
% Plotting average time courses as percent change from baseline
% (base_period

base_period=30; % in seconds
base_period_pt=base_period/(time_avg(1,2)-time_avg(1,1));

clear delta_aDb1_avg time_avg_new aDb1_avg_new

% time_avg_new=time_avg(1,Start_pt:End_pt)-time_avg(1,Start_pt-1);
% aDb1_avg_new=aDb1_avg(:,Start_pt:End_pt);

delta_aDb1_avg(1,:)=100*(aDb1_avg(1,:)/mean(aDb1_avg(1,1:base_period_pt))-1); 
delta_aDb1_avg(2,:)=100*(aDb1_avg(2,:)/mean(aDb1_avg(2,1:base_period_pt))-1); 
delta_aDb1_avg(3,:)=100*(aDb1_avg(3,:)/mean(aDb1_avg(3,1:base_period_pt))-1); 
delta_aDb1_avg(4,:)=100*(aDb1_avg(4,:)/mean(aDb1_avg(4,1:base_period_pt))-1); 

figure()

txt4='LBNP -30 mmHg\rightarrow';
txt5='LBNP -80 mmHg\rightarrow';

subplot(4,1,1)
hold on
x_pos=[2, 2.5, 3, 3.5]; %task strat time in minutes

for i=1:size(x_pos,2)

rectangle('Position',[x_pos(i),1.1*min(delta_aDb1_avg(1,:)),.25,1.1*max(delta_aDb1_avg(1,:))-1.1*min(delta_aDb1_avg(1,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
    'LineWidth',3)
end

plot(time_avg/60,delta_aDb1_avg(1,:),'r','LineWidth',1)
title('{\itr}_{SD}=1 cm')
set(gca,'xlim',[0 max(time/60)])
ylabel('\DeltaBFi (%)')

text(1.4,0.00000001,txt5,'FontSize',6)

Range(1,1)=0.9*min(aDb1_avg(4,:));
Range(1,2)=1.1*max(aDb1_avg(4,:));

subplot(4,1,2)
hold on

% temp_axaDb=delta_aDb1_avg(1,:)

for i=1:size(x_pos,2)

rectangle('Position',[x_pos(i),1.1*min(delta_aDb1_avg(2,:)),.25,1.1*max(delta_aDb1_avg(2,:))-1.1*min(delta_aDb1_avg(2,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
    'LineWidth',3)
end

plot(time_avg/60,delta_aDb1_avg(2,:),'m','LineWidth',1)
title('{\itr}_{SD}=2.5 cm')
set(gca,'xlim',[0 max(time/60)])
ylabel('\DeltaBFi (%)')

subplot(4,1,3)
hold on

% temp_axaDb=delta_aDb1_avg(1,:)

for i=1:size(x_pos,2)

rectangle('Position',[x_pos(i),1.1*min(delta_aDb1_avg(3,:)),.25,1.1*max(delta_aDb1_avg(3,:))-1.1*min(delta_aDb1_avg(3,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
    'LineWidth',3)
end

plot(time_avg/60,delta_aDb1_avg(3,:),'b','LineWidth',1)
title('{\itr}_{SD}=2 cm')
set(gca,'xlim',[0 max(time/60)])
ylabel('\DeltaBFi (%)')

subplot(4,1,4)
hold on

% temp_axaDb=delta_aDb1_avg(1,:)

for i=1:size(x_pos,2)

rectangle('Position',[x_pos(i),1.1*min(delta_aDb1_avg(4,:)),.25,1.1*max(delta_aDb1_avg(4,:))-1.1*min(delta_aDb1_avg(4,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
    'LineWidth',3)
end

plot(time_avg/60,delta_aDb1_avg(4,:),'k','LineWidth',1)
title('{\itr}_{SD}=2.5 cm')
xlabel('Time (s)')
set(gca,'xlim',[0 max(time/60)])
ylabel('\DeltaBFi (%)')

%% Percent Change (BFi) FROM TOURNIQUET

clear delta_aDb1_avg time_avg_new aDb1_avg_new T_res

firts_task_st=x_pos(1)*60; % first LBNP in seconds
base_task_st=firts_task_st-base_period;
last_task_end=x_pos(4)*60+30-t_avg; % first LBNP in seconds

Start_pt=find(time_avg==base_task_st);
End_pt=find(time_avg==last_task_end);

time_avg_new=time_avg(1,Start_pt:End_pt)-time_avg(1,Start_pt-1);
T_res=time_avg_new(1,2)-time_avg_new(1,1);
aDb1_avg_new=aDb1_avg(:,Start_pt:End_pt);

delta_aDb1_avg(1,:)=100*(aDb1_avg_new(1,:)/mean(aDb1_avg_new(1,1:base_period_pt))-1); 
delta_aDb1_avg(2,:)=100*(aDb1_avg_new(2,:)/mean(aDb1_avg_new(2,1:base_period_pt))-1); 
delta_aDb1_avg(3,:)=100*(aDb1_avg_new(3,:)/mean(aDb1_avg_new(3,1:base_period_pt))-1); 
delta_aDb1_avg(4,:)=100*(aDb1_avg_new(4,:)/mean(aDb1_avg_new(4,1:base_period_pt))-1); 

figure('units','centimeters', 'Position',[7 1.5 25 18]); %18 height 30 width
txt1='0 mmHg';
txt2='50 mmHg';
txt3='180 mmHg';
txt4='LBNP -30 mmHg\rightarrow';
txt5='LBNP -80 mmHg\rightarrow';

subaxis(4,1,1,'SpacingVert',0.1,'SpacingHoriz',0.1,'MR',0.03, 'ML',0.1,'MT',0.07,'MB',0.11)
hold on

x_pos_new=(x_pos*60-(firts_task_st-base_period))/60;

for i=1:size(x_pos_new,2)

rectangle('Position',[x_pos_new(i),-100,.25,200],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
    'LineWidth',3)
end

plot(time_avg_new/60,delta_aDb1_avg(1,:),'r','LineWidth',1);
title('{\itr}_{SD}=1 cm')
set(gca,'xlim',[T_res/60 max(time_avg_new/60)])
ylabel('\DeltaBFi (%)')
yticks([-100 -75 -50 -25 0 25])
set(gca,'ylim',[-100 35],'Xticklabel',{});
set(gca,'ylim',[-100 35]);

subaxis(4,1,2)
hold on

for i=1:size(x_pos_new,2)

rectangle('Position',[x_pos_new(i),-100,.25,200],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
    'LineWidth',3)
end

plot(time_avg_new/60,delta_aDb1_avg(2,:),'m','LineWidth',1)
title('{\itr}_{SD}=1.5 cm')
set(gca,'xlim',[T_res/60 max(time_avg_new/60)])
ylabel('\DeltaBFi (%)')
yticks([-100 -75 -50 -25 0 25])
set(gca,'ylim',[-100 35],'Xticklabel',{});
% 
% 
subaxis(4,1,3)
hold on

for i=1:size(x_pos_new,2)

rectangle('Position',[x_pos_new(i),-100,.25,200],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
    'LineWidth',3)
end

plot(time_avg_new/60,delta_aDb1_avg(3,:),'b','LineWidth',1);
title('{\itr}_{SD}=2 cm')
set(gca,'xlim',[T_res/60 max(time_avg_new/60)])
ylabel('\DeltaBFi (%)')
yticks([-100 -75 -50 -25 0 25])
set(gca,'ylim',[-100 35],'Xticklabel',{});

subaxis(4,1,4)
hold on

for i=1:size(x_pos_new,2)
    rectangle('Position',[x_pos_new(i),-100,.25,200],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
    'LineWidth',3)
end

plot(time_avg_new/60,delta_aDb1_avg(4,:),'k','LineWidth',1);
title('{\itr}_{SD}=2.5 cm')
xlabel('Time (min)')
xticks([0 1 2 3 4 5 6 7 8])
set(gca,'xlim',[T_res/60 max(time_avg_new/60)])
ylabel('\DeltaBFi (%)')
yticks([-100 -75 -50 -25 0 25])
set(gca,'ylim',[-100 35],'xtick',0.5:0.5:5);

%% Data[ DCS ]
% Ensemble Average for -80 LBNP cycles - DCS

clear time_new
clear delta_aDb1_avg_interp delta_aDb1_temp

Baseline=5; %in seconds
Start_t=x_pos_new(1)*60-Baseline; % in seconds
Period_t=30-t_avg; %in seconds
No_of_rep=4;

% T_res_new=time_new(1,2)-time_new(1,1);
Baseline_pt=Baseline/t_avg;
Period_t_pt=Period_t/t_avg;

Start_t_pt= find(time_avg_new >= Start_t,1,'first');

%DCS analysis
for chan=1:4
    for i=1:No_of_rep
        delta_aDb1_temp(chan,i,:)=delta_aDb1_avg(chan,Start_t_pt+Period_t_pt*(i-1):Start_t_pt+Period_t_pt*(i));
        delta_aDb1_temp(chan,i,:)=delta_aDb1_temp(chan,i,:)-squeeze(mean(delta_aDb1_temp(chan,i,1:Baseline_pt)));
    end
end

figure('units','centimeters', 'Position',[5 2 35 12]); %18 height 30 width


for chan=1:4
    Ch=[1 1.5 2 2.5];

    subaxis(1,4,chan,'SpacingVert',0.01,'SpacingHoriz',0.04,'MR',0.03, 'ML',0.06,'MT',0.09,'MB',0.15)
    hold on
    rectangle('Position',[5,-200,15,350],'FaceColor',[0.9 .9 .9],'EdgeColor','none',...
    'LineWidth',3)

    time_avg_new1=T_res*(1:1:size(delta_aDb1_temp,3));
    DCS_avg=squeeze(mean(delta_aDb1_temp(chan,:,:),2))';
    temp=squeeze(delta_aDb1_temp(chan,:,:));
    DCS_std=std(temp);

    lineProps.col{1} = 'k';
    lineProps.edgestyle = ':';
    lineProps.linewidth = 10;
    mseb(time_avg_new1,DCS_avg,DCS_std,lineProps,1)

    tit=strcat('{\itr}_{SD}=', num2str(Ch(chan)),'cm')
    title(tit)
    
    if chan==1
        ylabel('\DeltaBFi (%)')
    else
        set(gca,'yticklabel',{})
    end
    xlabel('Time (s)')
    
    set(gca,'ylim',[-100 100])
    set(gca,'xlim',[min(time_avg_new1) max(time_avg_new1)])
end
sgtitle('LBNP at -80 mmHg', 'Color', 'Red', 'FontSize', 10, 'EdgeColor', 'k');

%% Saving objects
DCS1_1rep=squeeze(delta_aDb1_temp(1,1,:))
DCS15_1rep=squeeze(delta_aDb1_temp(2,1,:))
DCS2_1rep=squeeze(delta_aDb1_temp(3,1,:))
DCS25_1rep=squeeze(delta_aDb1_temp(4,1,:))

DCS_1rep=horzcat(DCS1_1rep,DCS15_1rep,DCS2_1rep,DCS25_1rep)
%% save object file
save('MPCM009_cuffLBNP_1r.mat','delta_aDb1_temp','DCS_1rep','time_avg_new1')
