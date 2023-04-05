clear all;
close all
%--------------------------------------------------------------------------
% This code is for tourniquet pulsatility study. It has multiple sections
% for preprocessing the data and getting the ensembled averaged pulsatile
% shape of the signal. The preprocessing includes LPF, g2 averaging. 
% Folder name  - provide folder name for which you want to convert the data
% Folder='21.4.27-apnea without cuff';
file_no = 6;
%--------------------------------------------------------------------------
% 
filename_d=strcat('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\DCS\20220608 - '+string(file_no),'\Data.mat');

load(filename_d)
% load("D:\Jignesh\MSc Western Uni\Research\Data\20220706\20220706\Data.mat")
dcs = Data;



filename_nd=strcat('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\DCS_T\20220608 - '+string(2)+'TR','\Data.mat');
load(filename_nd)
% load("D:\Jignesh\MSc Western Uni\Research\Data\20220706\20220706 - TRTP\Data.mat")
dcs_t = Data;

% dcs_nt = dcs_n(128:end,:,:);
% dcs_t = dcs(121:121+size(dcs_nt,1)-1,:,:);
% Data_all = [dcs  dcs_t];

%Loading the ECG file
filename_nd=strcat('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\ECG\20220608_'+string(file_no)+'.mat');
load(filename_nd)

%%
l = 60; % length of signal as ROI in seconds
s = 1; % Signal ROI starting point in seconds;

ecg_a = data(datastart(1):dataend(1));
bp_a = data(datastart(2):dataend(2));
% tcd_a = data(datastart(3):dataend(3));

ecg1 = ecg_a(1:660000);
ecg1 = normalize(ecg1);
% ecg1 = lpf(ecg1,5,1000);
    
% tcd = tcd_a(1:length(ecg1));
% % tcd = normalize(tcd);
% tcd = lpf_ffilts(tcd,30,1000);

bp_a = bp_a(1:length(ecg1));
bp_a = lpf(bp_a,3,1000);

%%
% dcs = Data;
g2(1,:,:)=squeeze(dcs(:,1,:)-1); %g2-1 curve generation
g2_1_temp=squeeze(dcs(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(dcs(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(dcs(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(dcs(:,4,:)-1); %g2-1 curve generation
% g2_5_temp=squeeze(dcs_t(:,1,:)-1);
% g2_6_temp=squeeze(dcs_t(:,2,:)-1);
% g2_7_temp=squeeze(dcs_t(:,3,:)-1);
% g2_8_temp=squeeze(dcs_t(:,4,:)-1);
% average g2 curve for large source detector separation
% for i=1:size(g2,2)
%     g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:)+g2_5_temp(i,:)+g2_6_temp(i,:)+g2_7_temp(i,:)+g2_8_temp(i,:))/8;
% end
for i=1:size(g2,2)
    g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
end

% for i=1:size(g2,2)
%     g2(3,i,:)=( g2_5_temp(i,:)+g2_6_temp(i,:)+g2_7_temp(i,:) + g2_8_temp(i,:))/4;
% end

% for i=1:size(g2,2)
%     g2a(1,i,:)=( normalize(g2_2_temp(i,:))+normalize(g2_3_temp(i,:))+normalize(g2_4_temp(i,:))+ normalize(0.25*g2_5_temp(i,:))+ 0.25*normalize(g2_6_temp(i,:)) ...
%         + 0.25*normalize(g2_7_temp(i,:)) + 0.25*normalize(g2_8_temp(i,:)))/4;
% end
% for i=1:size(g2,2)
%     g2a(1,i,:)=( (g2_2_temp(i,:)./g2_2_temp(1,1))+(g2_3_temp(i,:)./g2_3_temp(1,1))+(g2_4_temp(i,:)./g2_4_temp(1,1))+ (g2_5_temp(i,:)./g2_5_temp(1,1))+ (g2_6_temp(i,:)./g2_6_temp(1,1)) ...
%         + (g2_7_temp(i,:)./g2_7_temp(1,1)) + (g2_8_temp(i,:)./g2_8_temp(1,1)))/7;
% end


%% Processing the g2 curves to test the beta values in fibres of the standalone and hybrid dcs system

g2_temp(:,1,:) = g2(1,:,:);
g2_temp(:,2,:) = g2(2,:,:);
g2_2 = g2(2,:,:);
 for i=1:length(g2_2)
    g_b=squeeze(g2_2(1,i,:))';
    g_b=g_b-mean(g_b(1,end-10:end));
    g_b=(g_b/max(g_b))/2+1;
    g_b(1) = round(g_b(1),1);
    g2_temp(i,2,:) = g_b;
% 
%     bsl_art=mean(g2_p(1,bsl,end));
% 
%     g_b=squeeze(g2(1,i,:))';
%     g_b=g_b-g_b(1,end);
%     g_b=(g_b/max(g_b));
% %     g_b(1) = round(g_b(1),1);
%     g2(1,i,:) = g_b*0.5+1;

 
%     g_p=squeeze(g2(2,i,:))';
%     g_p=g_p-mean(g_p(1,end-10:end));
%     g_p=(g_p/max(g_p))/2+1;
%     g_p(1) = round(g_p(1),1);
%     g2(2,i,:) = g_p;
%     bsl_art=mean(g2_p(2,bsl,end));
% 
%     g_b=squeeze(g2(2,i,:))';
%     g_b=g_b-g_b(1,end);
%     g_b=(g_b/g_b(1,1));
% %     g_b(1) = round(g_b(1),1);
%     g2(2,i,:) = g_b*0.5+1;
 end
 %% Plotting the g2 transform as above
 semilogx(Data_tau, squeeze(g2_temp(10,2,:)));
 hold on;
 semilogx(Data_tau, squeeze(dcs(10,4,:)));
title("g2_transform for low beta values");
xlabel("Tau (s)");


figure();
aDb_r = standalone_dcs(dcs,Data_tau);
aDb_n = standalone_dcs(g2_temp,Data_tau);
figure()
plot(aDb_r(2,:));
hold on;
plot(aDb_n(2,:));
legend("Raw_2.5","Transformed_2.5");
xlabel("time(s)");
ylabel("\alpha_D_b")
%% Crating the average g2 curves for the tourniquet pulsatility experiment to check if both the systems are working correctly
hfig = figure()
avg_g2 = mean(g2(:,1:30,:),2);
fname = "avg_g2_figure"
% a_g2 = mean(Data_all(1:120*20,[1 2 5],:)-1,1);
% avg_g2(:,1,:) = a_g2(1,:,:);
% semilogx(Data_tau,squeeze(avg_g2(1,1,:))/max(squeeze(avg_g2(1,1,:))));
% hold on
% semilogx(Data_tau,squeeze(avg_g2(2,1,:))/max(squeeze(avg_g2(2,1,:))));
% semilogx(Data_tau,squeeze(avg_g2(3,1,:))/max(squeeze(avg_g2(3,1,:))));

% semilogx(Data_tau,squeeze(avg_g2(1,1,:)));
% hold on
semilogx(Data_tau,squeeze(avg_g2(2,1,:)+1));
% semilogx(Data_tau,squeeze(avg_g2(3,1,:)));

% legend("Normalized 1cm","Norm 2.5cm","Norm 2.5cm TR","DCS 1cm","DCS 2.5cm","TR sys DCS 2.5cm ")
title("Avg. Autocorrelatio Curve")
xlabel("Data tau")
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-painters')



%% aDb calculations

% dcs = Data;
g2(1,:,:)=squeeze(dcs(:,1,:)-1); %g2-1 curve generation
g2_1_temp=squeeze(dcs(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(dcs(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(dcs(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(dcs(:,4,:)-1); %g2-1 curve generation
% g2_5_temp=squeeze(dcs_t(:,1,:)-1);
% g2_6_temp=squeeze(dcs_t(:,2,:)-1);
% g2_7_temp=squeeze(dcs_t(:,3,:)-1);
% g2_8_temp=squeeze(dcs_t(:,4,:)-1);
% average g2 curve for large source detector separation
% for i=1:size(g2,2)
%     g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:)+g2_5_temp(i,:)+g2_6_temp(i,:)+g2_7_temp(i,:)+g2_8_temp(i,:))/8;
% end
for i=1:size(g2,2)
    g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
end

rho = [0.7 2.5]; %source detector separations in cm 
mua = 0.05; %cm^-1 baseline absorption coefficient
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
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 500000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rsd,beta);
        aDb1(chan,i) = FittedParams(1);
    end
end

%% Trial code
% rho = [1 1.5 2 2.5]; %source detector separations in cm 
% mua = 0.17; %cm^-1 baseline absorption coefficient
% mus = 10; %cm^-1 baseline reduced scattering coefficient
% 
% tau_values=Data_tau;
% 
% for chan=1:size(g2a,1)
%      for i=1:size(g2a,2) 
%         rsd=rho(chan);
%         g2_temp(i,:)=squeeze(g2a(chan,i,:));
%         LB = [0];
%         UB = [inf];
%         Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
%         beta= squeeze(g2a(chan,i,1)); %0.1568;
%         options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
%         [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rsd,beta);
%         aDb1a(chan,i) = FittedParams(1);
%     end
% end

%%
close all
Channel=2;
Curve_no=10;
rho = [1 2.5];

beta=g2(Channel, Curve_no,1);
aDb_fit=aDb1(Channel,Curve_no);

g2_fit=gen_DCS_fit(Data_tau,mua,mus,rho(Channel),beta,aDb_fit);

semilogx(Data_tau,squeeze(g2(Channel,Curve_no,:)),'k')
hold on
semilogx(Data_tau,g2_fit,'r')

%% plotting the fit for the whole signal
close all
for i=2400:8000
    Channel=2;
Curve_no=i;
rho = [0.7 2.5];

beta=g2(Channel, Curve_no,1);
aDb_fit=aDb1(Channel,Curve_no);

g2_fit=gen_DCS_fit(Data_tau,mua,mus,rho(Channel),beta,aDb_fit);

semilogx(Data_tau,squeeze(g2(Channel,i,:)),'k')
hold on
semilogx(Data_tau,g2_fit,'r')
ylim([0 1])

hold off
pause(0.01)
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
title('{\itr}_{SD}=2.5 cm')
xlabel('Time (s)')

% subplot(3,1,3)
% 
% plot(time,aDb1(2,:))
% title('{\itr}_{SD}=2.5 cm TR system')
% xlabel('Time (s)')

%% 
plot(normalize(aDb1'));
mean_adb = zeros(5,3,660);

%% Fiding the shift by checking for the SNR
corr = zeros([1,2*size(aDb1,2)-1]);
for i=1:100
    adb_25 = aDb1(2,:)+circshift(aDb1(2,:),i-50);
    s(i) = snr(adb_25,20);
    corr(i,:) = xcorr(aDb1(2,:),circshift(aDb1(3,:),i-50));
end
%% Plotting the normalize signals for comparision
close all;
t_res=0.05; % seconds
t_avg=1; % window used for averaging in seconds

t_avg_pt=t_avg/t_res; % window used for averaging in points
time=t_res*(1:1:size(aDb1,2));
% Data_time(1,i)=i*t_res;
%     aDb1 = adb_avg
j=1;
for i=2:t_avg_pt:size(aDb1,2)
    aDb1_avg(1,j)=mean(aDb1(1,i-1:i+t_avg_pt-2));
    aDb1_avg(2,j)=mean(aDb1(2,i-1:i+t_avg_pt-2));
%     aDb1_avg(3,j)=mean(aDb1(3,i-1:i+t_avg_pt-2));
%     time_avg(1,j)=(Data_time(1,i+t_avg_pt-1));
    j=j+1;
end

% mean_adb(file_no/2,:,:) = aDb1_avg;

t = (1:1:length(aDb1(1,:)))/20;

x_pos=[0.1,4,7,11]; %task strat time in minutes
txt = ["BSL","25%", "50%","150%","Rec"];
% txt = ["80 mmHg", "160mmHg","BSL"];
% for i=1:size(x_pos,2)
% if mod(i,2)==0
%     rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,max(aDb1(1,:))],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
%         'LineWidth',3)
%     text(t(x_pos(i)*1200)*1.05,0.9*max(aDb1(1,:)),txt(i));
% else
%     rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,max(aDb1(1,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
%         'LineWidth',3)
%     text(t(x_pos(i)*1200)*1.05,0.9*max(aDb1(1,:)),txt(i));
% 
% end
% end
% x_pos=[5,8]; %task strat time in minutes
% % txt = ["25%", "50%","150%","BSL"];
% txt = ["80 mmHg", "160mmHg"];
% for i=1:size(x_pos,2)
% if mod(i,2)==0
%     rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,240,max(aDb1(1,:))],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
%         'LineWidth',3)
%     text(t(x_pos(i)*1200)*1.05,0.9*max(aDb1(1,:)),txt(i));
% else
%     rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,180,max(aDb1(1,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
%         'LineWidth',3)
%     text(t(x_pos(i)*1200)*1.05,0.9*max(aDb1(1,:)),txt(i));
% 
% end
% end

hold on;
plot(t,aDb1','DisplayName','aDb1');
hold on;
title("CBFi w.r.t Tourniquet Pressure Levels");

xlabel("Time(s)");
ylabel("CBFi");
plot(aDb1_avg','DisplayName','aDb1_avg','LineWidth',1.5)
legend("DCS 1cm","DCS 2.5cm"); 



figure();

for i=1:size(aDb1_avg,1)
    adb_avg_smooth(i,:) = smooth(aDb1_avg(i,:),3,'moving');
end
hold on;

% per_ch_mean(1,:) = 100*(adb_avg_smooth(1,:)./ mean(adb_avg_smooth(1,1:60),2));
% per_ch_mean(2,:) = 100*(adb_avg_smooth(2,:)./ mean(adb_avg_smooth(2,1:60),2));

per_ch_mean(1,:) = 100*((adb_avg_smooth(1,:)-mean(adb_avg_smooth(1,1:60),2))./ mean(adb_avg_smooth(1,1:60),2));
per_ch_mean(2,:) = 100*((adb_avg_smooth(2,:)-mean(adb_avg_smooth(2,1:60),2))./ mean(adb_avg_smooth(2,1:60),2));

% per_ch_mean(3,:) = 100*(adb_avg_smooth(3,:)./ mean(adb_avg_smooth(3,1:100),2));

t = (1:1:length(aDb1(1,:)))/20;

x_pos=[3,5,7,9]; %task strat time in minutes
% txt = ["25%", "50%","150%","BSL"];
% for i=1:size(x_pos,2)
% if mod(i,2)==0
%     rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,max(per_ch_mean(1,:))],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
%         'LineWidth',3)
%     text(t(x_pos(i)*1200)*1.05,0.95*max(per_ch_mean(1,:)),txt(i));
% else
%     rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,max(per_ch_mean(1,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
%         'LineWidth',3)
%     text(t(x_pos(i)*1200)*1.05,0.95*max(per_ch_mean(1,:)),txt(i));
% 
% end
% end
% x_pos=[5,8]; %task strat time in minutes
% txt = ["25%", "50%","150%","BSL"];
% txt = ["80 mmHg", "160mmHg"];
% for i=1:size(x_pos,2)
% if mod(i,2)==0
%     rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,240,max(per_ch_mean(1,:))],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
%         'LineWidth',3)
%     text(t(x_pos(i)*1200)*1.05,0.95*max(per_ch_mean(1,:)),txt(i));
% else
%     rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,180,max(per_ch_mean(1,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
%         'LineWidth',3)
%     text(t(x_pos(i)*1200)*1.05,0.95*max(per_ch_mean(1,:)),txt(i));
% 
% end
% end

plot(per_ch_mean')
title("Percent change in CBFi");

xlabel("Time(s)");
ylabel("CBFi");
legend("DCS 1cm","DCS 2.5cm"); 

%% normalizing the average signal
adb(1,:) = (adb_avg_smooth(1,:) -min(adb_avg_smooth(1,1:50)))./range(adb_avg_smooth(1,1:50))
adb(2,:) = (adb_avg_smooth(2,:) -min(adb_avg_smooth(1,1:50)))./range(adb_avg_smooth(1,1:50))
plot(adb')


%%

hfig = figure();
picturewidth = 15; % set this parameter and keep it forever
hw_ratio = 0.5; % feel free to play with this ratio
t = (1:1:length(adb_avg_smooth))/20;
plot(t,adb_avg_smooth(1,:))
xlabel("Time(s)")
ylabel("BFi")



set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

% set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%% 
% for i=150:1500:size(g2_5_temp,1)
%     semilogx(Data_tau,g2_5_temp(i,:))
%     title("g2 curves for DCS(hybrid) channel")
%     xlabel("Data Tau")
%     hold on;
% end

% semilogx(Data_tau,g2_2_temp(250,:))
% hold on;
% semilogx(Data_tau,g2_6_temp(250,:))
% semilogx(Data_tau,g2_7_temp(250,:))
% semilogx(Data_tau,g2_8_temp(250,:))
% legend("Chan 5","6","7","8")
% title("g2 curves for DCS(hybrid) channel")
% xlabel("Data Tau")

%% % change in the pulsatility of the signal
%This code has been updated for global analysis for all the subjects
%together. Go to next segment run that and then come back.
close all;
per_ch_g_pul = zeros(6,3,13200);
per_change = zeros(3,13200);
for b=1:6
    aDb1  = squeeze(adb(b,:,:));
    for i=1:size(aDb1,1)
        s = aDb1(i,:);
        breakpoints = 1:200:length(s);
        y = detrend(s,1,breakpoints);
    %     plot(y);
        
        % Plotting the mean average smoothing curve with the envelope
        x = (1:1:length(y))*0.05;
        [envHigh, envLow] = envelope(y,100,'rms');
        envMean = (envHigh+envLow)/2;
        env_diff = envHigh - envLow;
        figure();
        plot(x,y,x,envLow,x,envHigh,x,envMean);
        mh = max(envHigh);
        ml = max(envLow);
        mx = max(mh,ml);
    %         
        snr_val(i,1) = db2mag(snr(s(1:3600),20));
        snr_val(i,2) = db2mag(snr(s(3601:6000),20));
        snr_val(i,3) = db2mag(snr(s(6001:8400),20));
        snr_val(i,4) = db2mag(snr(s(8401:10800),20));
        snr_val(i,5) = db2mag(snr(s(10801:end),20));
       
        
        xt = [10,40,135,240,265];
        yt = [mx,mx-1,mx-2,mx-1,mx];
    %     text(xt,yt,str,'Color','red','FontSize',14);
    %     text((max(x)/2)-10, mx,s_str, 'FontSize',16)
        
        
        title("Pulsatility envelope");
        xlabel("Time (s)");
        ylabel("DCS 1cm aDb");
        
    %         figure();
        env_diff_smooth = smooth(env_diff,150);
    %         plot(x,env_diff_smooth);
    %         title("Pulse amplitude normalized "+s_str+ " "+file_list(j));
    %         xlabel("Time (s)");
    %         ylabel("Pulse Amplitude");
        % % Change in the pulsatility
        avg_sig = smooth(s,150);
        per_change(i,:) = (100*(env_diff_smooth')./mean(env_diff_smooth(1:3500))) ;
        per_ch_g_pul(b,i,:) = per_change(i,:);
    %         plot(x,squeeze(per_change(j,i,:)));
    %         title("Percent change in pulse amplitude with pressure cuff " +s_str);
    %         xlabel("Time (s)");
    %         ylabel("% change in Pulse Amplitude");
    end
    % Plotting the changes in each channel in single plot
    figure();
    t = (1:1:length(aDb1(1,:)))/20;
    
    x_pos=[3,5,7,9]; %task strat time in minutes
    txt = ["25%", "50%","150%","BSL"];
    for i=1:size(x_pos,2)
    if mod(i,2)==0
        rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,max(per_change(1,:))],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
            'LineWidth',3)
        text(t(x_pos(i)*1200)*1.05,0.9*max(per_change(1,:)),txt(i));
    else
        rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,max(per_change(1,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
            'LineWidth',3)
        text(t(x_pos(i)*1200)*1.05,0.9*max(per_change(1,:)),txt(i));
    
    end
    end
    
%     x_pos=[5,8]; %task strat time in minutes
%     % txt = ["25%", "50%","150%","BSL"];
%     txt = ["80 mmHg", "160mmHg"];
%     for i=1:size(x_pos,2)
%     if mod(i,2)==0
%         rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,240,max(per_ch_mean(1,:))],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
%             'LineWidth',3)
%         text(t(x_pos(i)*1200)*1.05,0.95*max(per_ch_mean(1,:)),txt(i));
%     else
%         rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,180,max(per_ch_mean(1,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
%             'LineWidth',3)
%         text(t(x_pos(i)*1200)*1.05,0.95*max(per_ch_mean(1,:)),txt(i));
%     
%     end
%     end
    
    hold on;
    l = length(30:length(per_change)-10);
    x_pc = (1:1:l)/20;
    plot(x_pc,squeeze(per_change(1,30:end-10)),'DisplayName','per_change_dcs_1');
    hold on;
    plot(x_pc,squeeze(per_change(2,30:end-10)),'DisplayName','per_change_dcs_25');
    plot(x_pc,squeeze(per_change(3,30:end-10)),'DisplayName','per_change_dcs_25_TR');
    hold off;
    legend("DCS 1cm","DCS 2.5cm","DCS 2.5cm TR"); 
    title("Percent change in pulse amplitude with pressure cuff");
    xlabel("Time (s)");
    ylabel("% change in Pulse Amplitude");
end
%%
% Mean and STD of the percent change in amplitude.
mean_per_change_ampg = squeeze(mean(per_ch_g_pul,1));
std_per_change_ampg = squeeze(std(per_ch_g_pul,0,1));
ind = [1 3];
color = ['b','r'];
hfig = figure();
fname = "Percent_change_pulse_amplitude_tourniquet";
subplot(2,1,1)

t = (1:1:length(mean_per_change_ampg(1,:)))/20;
x_pos=[3,5,7,9]; %task strat time in minutes
txt = ["25%", "50%","150%","BSL"];
curve1 = mean_per_change_ampg(3,:)+std_per_change_ampg(3,:);
for i=1:size(x_pos,2)
if mod(i,2)==0
    rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,1.1*max(curve1)],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
        'LineWidth',3)
    text(t(x_pos(i)*1200)*1.05,0.65*max(curve1),txt(i));
else
    rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,1.1*max(curve1)],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
        'LineWidth',3)
    text(t(x_pos(i)*1200)*1.05,0.65*max(curve1),txt(i));

end
end
hold on;
for i=1:2
%     subplot(2,1,i)
    d = ind(i);
    y = mean_per_change_ampg(d,:);
    x = (1:numel(y))/20;
    std_dev = std_per_change_ampg(d,:);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween,color(i),'FaceAlpha',0.1);
    hold on;
    plot(x, y, 'LineWidth', 2);
    xlabel("Time(s)")
    
end
ylabel("Percent(%)")
legend("1cm-STD","1cm-Mean","2.5cm-STD","2.5cm-Mean");

subplot(2,1,2)
% figure();


x_pos=[3,5,7,9]; %task strat time in minutes
txt = ["25%", "50%","150%","BSL"];
for i=1:size(x_pos,2)
if mod(i,2)==0
    rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,1.1*max(curve1)],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
        'LineWidth',3)
%     text(t(x_pos(i)*1200)*1.05,0.95*max(curve1),txt(i));
else
    rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,1.1*max(curve1)],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
        'LineWidth',3)
%     text(t(x_pos(i)*1200)*1.05,0.95*max(curve1),txt(i));

end
end
hold on;

t = (1:1:length(per_change_g(1,:)));
plot(t, per_change_g(1,:)', 'LineWidth', 1.5);
plot(t, per_change_g(3,:)', 'LineWidth', 1.5);
ylabel("aDb (CBFi)")
xlabel("Time(s)")
title("Percent change in CBFi");
legend("rsd = 1 cm","rsd = 2.5cm")



sgtitle("Percentage Change in Pulse amplitude")
xlabel("Time(s)")

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 1.15; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-painters')
%% Reading all the data for global analysis

% a = [2,4,6,10,12,14,16,18];
a = [2];
dcs_a = zeros(5,size(Data,1),size(Data,2),size(Data,3));
dcs_at = zeros(5,size(Data,1),size(Data,2),size(Data,3));
adb = zeros(5,3,13200);
mean_adb = zeros(5,3,660);
for b = 1:length(a)
    file_no = a(b);
    %--------------------------------------------------------------------------
    % 
    filename_d=strcat('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\DCS\20220608 - '+string(file_no),'\Data.mat');
    
    load(filename_d)
    % load("C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\TP\20220614\Data.mat")
    dcs_a(b,:,:,:) = Data;
    dcs = Data;
    
    
    
%     filename_nd=strcat('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\DCS_T\20220608 - '+string(file_no)+'TR','\Data.mat');
%     load(filename_nd)
%     % load("C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\TP\20220614 -TR\Data.mat")
%     dcs_at(b,:,:,:) = Data;
%     dcs_t = Data;
    
    % dcs_nt = dcs_n(128:end,:,:);
    % dcs_t = dcs(121:121+size(dcs_nt,1)-1,:,:);
%     Data_all = [dcs  dcs_t];

    g2(1,:,:)=squeeze(dcs(:,1,:)-1); %g2-1 curve generation
    g2_1_temp=squeeze(dcs(:,1,:)-1); %g2-1 curve generation
    g2_2_temp=squeeze(dcs(:,2,:)-1); %g2-1 curve generation
    g2_3_temp=squeeze(dcs(:,3,:)-1); %g2-1 curve generation
    g2_4_temp=squeeze(dcs(:,4,:)-1); %g2-1 curve generation
%     g2_5_temp=squeeze(dcs_t(:,1,:)-1);
%     g2_6_temp=squeeze(dcs_t(:,2,:)-1);
%     g2_7_temp=squeeze(dcs_t(:,3,:)-1);
%     g2_8_temp=squeeze(dcs_t(:,4,:)-1);
    % average g2 curve for large source detector separation
    % for i=1:size(g2,2)
    %     g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:)+g2_5_temp(i,:)+g2_6_temp(i,:)+g2_7_temp(i,:)+g2_8_temp(i,:))/8;
    % end
    for i=1:size(g2,2)
        g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
    end
    
%     for i=1:size(g2,2)
%         g2(3,i,:)=( g2_5_temp(i,:)+g2_6_temp(i,:)+g2_7_temp(i,:) + g2_8_temp(i,:))/4;
%     end
    
    % for i=1:size(g2,2)
    %     g2a(1,i,:)=( normalize(g2_2_temp(i,:))+normalize(g2_3_temp(i,:))+normalize(g2_4_temp(i,:))+ normalize(0.25*g2_5_temp(i,:))+ 0.25*normalize(g2_6_temp(i,:)) ...
    %         + 0.25*normalize(g2_7_temp(i,:)) + 0.25*normalize(g2_8_temp(i,:)))/4;
    % end
%     for i=1:size(g2,2)
%         g2a(1,i,:)=( (g2_2_temp(i,:)./g2_2_temp(1,1))+(g2_3_temp(i,:)./g2_3_temp(1,1))+(g2_4_temp(i,:)./g2_4_temp(1,1))+ (g2_5_temp(i,:)./g2_5_temp(1,1))+ (g2_6_temp(i,:)./g2_6_temp(1,1)) ...
%             + (g2_7_temp(i,:)./g2_7_temp(1,1)) + (g2_8_temp(i,:)./g2_8_temp(1,1)))/7;
%     end


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
    adb(b,:,:) = aDb1;

    close all;
    t_res=0.05; % seconds
    t_avg=1; % window used for averaging in seconds
    
    t_avg_pt=t_avg/t_res; % window used for averaging in points
    time=t_res*(1:1:size(aDb1,2));
    Data_time(1,i)=i*t_res;
    j=1;
    for i=2:t_avg_pt:size(aDb1,2)
        aDb1_avg(1,j)=mean(aDb1(1,i-1:i+t_avg_pt-2));
        aDb1_avg(2,j)=mean(aDb1(2,i-1:i+t_avg_pt-2));
%         aDb1_avg(3,j)=mean(aDb1(3,i-1:i+t_avg_pt-2));
    %     time_avg(1,j)=(Data_time(1,i+t_avg_pt-1));
        j=j+1;
    end
    
    mean_adb(b,:,:) = aDb1_avg;


end
%%
close all;
for i=1:6
    figure();
    plot(squeeze(adb(i,1,:)));
    title("subject "+i)
end

%% Global average over all the participants
close all;

mean_adb_cut = mean_adb([1,3,7],:,:);
adb_g_avg = squeeze(mean(mean_adb_cut,1));
adb_g_std = squeeze(std(mean_adb_cut,0,1));
str = ["r_sd = 1cm","r_sd = 2.5cm","r_sd = 2.5cm-TR"];

per_change_g(1,:) = 100*(adb_g_avg(1,:)./mean(adb_g_avg(1,1:170)));
per_change_g(2,:) = 100*(adb_g_avg(2,:)./mean(adb_g_avg(2,1:170)));
per_change_g(3,:) = 100*(adb_g_avg(3,:)./mean(adb_g_avg(3,1:170)));

% per_std_g(1,:) = 100*(adb_g_std(1,:)./mean(adb_g_std(1,1:170)));
% per_std_g(2,:) = 100*(adb_g_std(2,:)./mean(adb_g_std(2,1:170)));
% per_std_g(3,:) = 100*(adb_g_std(3,:)./mean(adb_g_std(3,1:170)));

per_std_g(1,:) = std(mean_adb_cut(1,:));
per_std_g(2,:) = std(per_change_g(2,:));


figure();
t = (1:1:length(aDb1(1,:)))/20;

x_pos=[3,5,7,9]; %task strat time in minutes
txt = ["25%", "50%","150%","BSL"];
for i=1:size(x_pos,2)
if mod(i,2)==0
    rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,max(per_change_g(1,:))],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
        'LineWidth',3)
    text(t(x_pos(i)*1200)*1.05,0.9*max(per_change_g(1,:)),txt(i));
else
    rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,max(per_change_g(1,:))],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
        'LineWidth',3)
    text(t(x_pos(i)*1200)*1.05,0.9*max(per_change_g(1,:)),txt(i));

end
end
hold on;
plot(per_change_g','DisplayName','per_change_g');
title("% Change in CBFi across all participant")
xlabel("Time(s)")
ylabel("% Change")
legend("r_sd = 1cm","r_sd = 2.5cm","r_sd = 2.5cm-TR")

figure();
for j=1:size(adb_g_avg,1)
    subplot(3,1,j);

    y = adb_g_avg(j,:);
    x = 1:numel(y);
    std_dev = adb_g_std(j,:);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];


    t = (1:1:length(aDb1(1,:)))/20;

    x_pos=[3,5,7,9]; %task strat time in minutes
    txt = ["25%", "50%","150%","BSL"];
    for i=1:size(x_pos,2)
    if mod(i,2)==0
        rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,1.1*max(curve1)],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
            'LineWidth',3)
        text(t(x_pos(i)*1200)*1.05,0.95*max(curve1),txt(i));
    else
        rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,1.1*max(curve1)],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
            'LineWidth',3)
        text(t(x_pos(i)*1200)*1.05,0.95*max(curve1),txt(i));
    
    end
    end
    hold on;
    
    fill(x2, inBetween,'k','FaceAlpha',0.2);
    hold on;
    plot(x, y, 'LineWidth', 2);
    ylabel("aDb (CBFi)")
    xlabel("Time(s)")
    title("Mean CBFi across all subjects @"+str(j))
    legend("STD","Mean")
end

%% Finding the mean and std of percent change of the BFi
close all
m_adb = mean_adb([1,2,3,4,6],1:2,:);
per_change = 100*(m_adb-mean(m_adb(:,:,20:170),3))./mean(m_adb(:,:,20:170),3);
plot(squeeze(per_change(:,2,:))')
figure();

% mean_per_change = mean(per_change,1);
mean_per_change = m_adb;

    bsl = 1:181;
    p1 = 181:301;
    p2 = 301:424;
    p3 = 424:540;
    rec = 540:660;

    sm_factor = 15;
    for i=1:size(m_adb,1)
        per_ch_sm(i,1,bsl) = (smooth(mean_per_change(i,1,bsl),sm_factor));
        per_ch_sm(i,1,p1) = (smooth(mean_per_change(i,1,p1),sm_factor));
        per_ch_sm(i,1,p2) = (smooth(mean_per_change(i,1,p2),sm_factor));
        per_ch_sm(i,1,p3) = (smooth(mean_per_change(i,1,p3),sm_factor));
        per_ch_sm(i,1,rec) = (smooth(mean_per_change(i,1,rec),sm_factor));
        per_ch_sm(i,2,bsl) = (smooth(mean_per_change(i,2,bsl),sm_factor));
        per_ch_sm(i,2,p1) = (smooth(mean_per_change(i,2,p1),sm_factor));
        per_ch_sm(i,2,p2) = (smooth(mean_per_change(i,2,p2),sm_factor));
        per_ch_sm(i,2,p3) = (smooth(mean_per_change(i,2,p3),sm_factor));
        per_ch_sm(i,2,rec) = (smooth(mean_per_change(i,2,rec),sm_factor));
    end
plot(squeeze(per_ch_sm(:,1,:))');
% plot(squeeze(mean_per_change)')


% figure();
% plot(normalize(squeeze(mean_per_change))')

figure();
per_change = 100*(per_ch_sm-mean(per_ch_sm(:,:,20:170),3))./mean(per_ch_sm(:,:,20:170),3);
plot(squeeze(per_change(:,1,:))')

figure();
mean_per_change = mean(per_change,1);
std_per_change = std(per_change,1);

    std_p_ch(1,1,:) = smooth(std_per_change(1,1,:),10)
     std_p_ch(1,2,:) = smooth(std_per_change(1,2,:),10)
plot(squeeze(mean_per_change)');
figure()
plot(squeeze(std_p_ch)')
hold on;
plot(squeeze(std_per_change)')

%% Takig average
 close all;
    t_res=0.05; % seconds
    t_avg=1; % window used for averaging in seconds
    
    t_avg_pt=t_avg/t_res; % window used for averaging in points
    time=t_res*(1:1:size(mean_pdfc,1));
%     Data_time(1,i)=i*t_res;

    j=1;
    for i=2:t_avg_pt:size(mean_pdfc,1)
        mpdfc_avg(j) = mean(mean_pdfc(i-1:i+t_avg_pt-2));
        j = j+1
    end
   mpdfc_avg = smooth(mpdfc_avg,6)
    j=1;
    for i=2:t_avg_pt:size(mean_pdfc,1)
        stdpfc_avg(j) = mean(std_pdfc(i-1:i+t_avg_pt-2));
        j = j+1
    end
    stdpfc_avg = smooth(stdpfc_avg,6)
    plot(stdpfc_avg)

%%
mpch =  squeeze(mean_per_change);
hfig = figure();
j = 1;

    y = mpch(j,:);
    x = 1:numel(y);
    std_dev = stdpch(j,:);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];


%     t = (1:1:length(aDb1(1,:)))/20;
    t = (1:1:length(mean_pdfc(:,1)))/20;


    x_pos=[3,5,7,9]; %task strat time in minutes
    txt = ["25%", "50%","150%","Recovery"];
    if j==1
        for i=1:size(x_pos,2)
            if i==1
                text(t(1200)*1.05,0.95*max(curve1),"BSL");
            end
        if mod(i,2)==0
            rectangle('Position',[t(x_pos(i)*1200),-90,120,180],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
                'LineWidth',3)
            text(t(x_pos(i)*1200)*1.05,0.95*max(curve1),txt(i));
        else
            rectangle('Position',[t(x_pos(i)*1200),-90,120,180],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
                'LineWidth',3)
            text(t(x_pos(i)*1200)*1.05,0.95*max(curve1),txt(i));
        
        end
        end
    end
    ylim([-90 60])
    xlim([0 660])
    hold all;
%     plot(x, mpdfc_avg)
    
    fill(x2, inBetween,'k','FaceAlpha',0.2);
    hold all;
    plot(x, y,'k', 'LineWidth', 2,'DisplayName','BFi_1 _c_m');
    ylabel("aDb (CBFi)")
    xlabel("Time(s)")
%     title("Mean CBFi across all subjects @"+str(j))
%     legend("STD","Mean")


j = 2;
    y = mpch(j,:);
    x = 1:numel(y);
    std_dev = stdpch(j,:);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];


%     t = (1:1:length(aDb1(1,:)))/20;

%     x_pos=[3,5,7,9]; %task strat time in minutes
%     txt = ["25%", "50%","150%","Recovery"];
%     if j==1
%         for i=1:size(x_pos,2)
%         if mod(i,2)==0
%             rectangle('Position',[t(x_pos(i)*1200),-80,120,150],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
%                 'LineWidth',3)
%             text(t(x_pos(i)*1200)*1.05,0.95*max(curve1),txt(i));
%         else
%             rectangle('Position',[t(x_pos(i)*1200),-80,120,150],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
%                 'LineWidth',3)
%             text(t(x_pos(i)*1200)*1.05,0.95*max(curve1),txt(i));
%         
%         end
%         end
%     end
    ylim([-90 60])
    xlim([0 660])
    hold all;
%     plot(x, mpdfc_avg)
    
    fill(x2, inBetween,'b','FaceAlpha',0.2);
    hold all;
    plot(x, y,'b', 'LineWidth', 2,'DisplayName','BFi_2_._5 _c_m');
    ylabel("aDb (CBFi)")
    xlabel("Time(s)")
%     title("Mean CBFi across all subjects @"+str(j))
%     legend("STD1","Mean1")



% mpdfc_avg = mean_pdfc;
hold all
y1 = mpdfc_avg';
x1 = 1:numel(y);
std_dev_1 = stdpfc_avg';
curve1 = y1 + std_dev_1;
curve2 = y1 - std_dev_1;
x2 = [x1, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
% plot(x1, mpdfc_avg,'r')
    
fill(x2, inBetween,'r','FaceAlpha',0.2);
hold all;
plot(x1, y1, 'r', 'LineWidth', 2,'DisplayName','CBFi');
ylabel("\Delta (%)")
xlabel("Time(s)")
% title("Mean CBFi across all subjects @"+str(j))
% legend show
% legend("BFi_1 _c_m Mean","BFi_1 _c_m STD","BFi_2.5 _c_m Mean","BFi_2.5 _c_m STD","CBFi Mean","CBFi STD")
% legend("BFi_1 _c_m Mean","BFi_2.5 _c_m Mean","CBFi Mean")

picturewidth = 35; % set this parameter and keep it forever
hw_ratio = 0.5; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%% Plotting the percent change in BFi across subject
for j=[1,2]
    figure();
    y = per_change_g(j,:);
    x = 1:numel(y);
    std_dev = per_std_g(j,:);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    
    
    t = (1:1:length(aDb1(1,:)))/20;
    
    x_pos=[3,5,7,9]; %task strat time in minutes
    txt = ["25%", "50%","150%","BSL"];
    for i=1:size(x_pos,2)
    if mod(i,2)==0
        rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,1.1*max(curve1)],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none',...
            'LineWidth',3)
        text(t(x_pos(i)*1200)*1.05,0.95*max(curve1),txt(i));
    else
        rectangle('Position',[t(x_pos(i)*1200),0.1*10^-9,120,1.1*max(curve1)],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none',...
            'LineWidth',3)
        text(t(x_pos(i)*1200)*1.05,0.95*max(curve1),txt(i));
    
    end
    end
    hold on;
    
    fill(x2, inBetween,'k','FaceAlpha',0.2);
    hold on;
    plot(x, y, 'LineWidth', 2);
    ylabel("aDb (CBFi)")
    xlabel("Time(s)")
    title("Mean CBFi across all subjects @"+str(j))
    legend("STD","Mean")
end

%%
for i=1:7
    figure()
    plot(squeeze(per_change_g(i,[1,2],:))')
end
%% Pulsatility analysis and g2 averaging of the signal

l = 60; % length of signal as ROI in seconds
s = 1; % Signal ROI starting point in seconds;

ecg_a = data(datastart(1):dataend(1));
bp_a = data(datastart(2):dataend(2));
% tcd_a = data(datastart(3):dataend(3));

ecg1 = ecg_a(1:660000);
ecg1 = normalize(ecg1);
ecg1 = lpf(ecg1,5,1000);
    
% tcd = tcd_a(1:length(ecg1));
% % tcd = normalize(tcd);
% tcd = lpf_ffilts(tcd,30,1000);

bp_a = bp_a(1:length(ecg1));
bp_a = lpf(bp_a,3,1000);

%% Finding the shift

close all;
ecg_ad = circshift(ecg1,-700);
x_e = (1:1:length(ecg_ad))/1000;
plot(x_e,normalize(ecg_ad),'b');
hold on;
x_d = (1:1:size(aDb1,2))/20;
plot(x_d,normalize(aDb1(1,:)),'r');
%% Plotting the ensemble average graph
[ens_avg_sig,pind] = ensemble_avg(ecg1(1:90000),aDb1(1,1:1800),600,1);
%% Averaging thr g2 curve
% close all;
% shift = 1:45:800;
% SNR = zeros(1,length(shift));
% bp_shift = circshift(bp_a,0);
% ecg_ad = circshift(ecg1,0); % Advancing the ECG signal to match the DCS signal
%% Finding the R-R peaks of ECG signal

close all;
% if exist("g2")
%     clear g2;
% end
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10; 

% g2 = Data;
% g2_avg = zeros(size(g2))*NaN;
    
% Defining the protocol ltime line


% ecg_1 = ecg_ad(1:300000);   % 5 min baseline data
% Data_1 = Data(1:6000,:,:);
% bp_1 = bp_shift(1:300000);
% ecg_2 = ecg_ad(300001:480000); % 3 min TNQT @80 mmHg
% Data_2 = Data(6001:9600,:,:);
% bp_2 = bp_shift(300001:480000);
% ecg_3 = ecg_ad(480001:720000); % 4 min TNQT @160 mmHg
% Data_3 = Data(9601:14400,:,:);
% bp_3 = bp_shift(480001:720000);
% ecg_4 = ecg_ad(720001:780000); % 1 min baseline
% Data_4 = Data(14401:15600,:,:);
% bp_4 = bp_shift(720001:780000);
% 
g2_backup = Data;
ecg = ecg_ad;
% g2_backup = g2;
g2_n = g2_backup;
strt_time =  [1,180,300,420,540];   %time in seconds
stp_time =  [180,310,420,540,660];    %time in seconds.

%  strt_time =  [1,32,62];
%  stp_time =  [28,58,90]; 

%For subject 5(Farah), as the protocol is a bit different;
% strt_time =  [1,300,480,720];   %time in seconds
% stp_time =  [300,480,720,780];    %time in seconds.
for m=1:length(strt_time)
    disp(m);
    clear g2;
    stop_time = stp_time(m);
    start_time = strt_time(m);
    g2 = g2_n(start_time*20:stop_time*20,:,:);
    ecg = ecg_ad(start_time*1000:stop_time*1000);
    % g2(:,1,:)=squeeze(D(:,1,:)-1); %g2-1 curve generation
    % g2_2_temp=squeeze(D(:,2,:)-1); %g2-1 curve generation
    % g2_3_temp=squeeze(D(:,3,:)-1); %g2-1 curve generation
    % g2_4_temp=squeeze(D(:,4,:)-1); %g2-1 curve generation
    
    % average g2 curve for large source detector separation
    % for i=1:size(g2,1)
    %     g2(i,2,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
    % end
    
    [h_pks,l_pks] = findpeaks(normalize(ecg),"MinPeakHeight",2.5,"MinPeakDistance",750);
    
    fig1=figure('units','centimeters', 'Position',[2 2 35 13]); %18 width 15 heigh
    hold on
    plot(normalize(ecg),'r')
    plot(l_pks, h_pks,'*k')
    hd_pks = ceil(h_pks./50);
    ld_pks = ceil(l_pks./50);
    % g2_avg = zeros(size(g2,1),size(g2,2),size(g2,3));
    avg_window_width = 60;
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

    g2_n(start_time*20:stop_time*20,:,:) = g2;
end

% for i=1:size(g2,1)
%     g2(i,2,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
% end



adb_avg = standalone_dcs(g2_n,Data_tau);
% SNR(1,k) = snr(adb_avg(2,:),20);
figure
snr(adb_avg(2,:),20)

% adb_avg = adb_avg; % This is cut out specific portion of the signal to plot
% figure();
% adb = standalone_tr_dcs(Data_all,Data_tau);
dcs_1lp = (lpf(adb_avg(1,:),7,20));
dcs_25lp =(lpf(adb_avg(2,:),7,20));
% dcs_25lp_tr = lpf(adb_avg(3,:),7,20);

% aDb1 = aDb1;
% figure();
% snr(adb_avg(2,:),20);
figure();
plot(adb_avg(2,:),'b',"LineWidth",1.5);
hold on; 
plot(aDb1(2,:),'r');
legend("Ensemle Temporal Averaged","Raw")
% title("Comparision of g2 averaging for cuff data MPCM004 width=50 cycles")
% legend("g2 Averaged signal","Raw signal")
% xlabel("samples (Time = samples/20)");
% ylabel("aDb")

%% Plotting the ensemble avg of the DCS signal
%send the  mode as wether the averae to be with NaN(1) or without the
%NaN(0)
close all;
%For subjects 1,2,3,6
strt_time =  [10,190,350,440,540];   %time in seconds
stp_time =  [160,290,400,510,660];    %time in seconds.
%For subject 5(Farah), as the protocol is a bit different;
% strt_time =  [10,320,590,720];   %time in seconds
% stp_time =  [280,460,700,760];    %time in seconds.

j =1;
[ens_avg_sig,pind,ensemble_curve] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_25lp(1,strt_time(j)*20:stp_time(j)*20),700,1);
d= ens_avg_sig(1:1000);

PI = (max(d)-min(d))./mean(d);

% ens_avg_sig = ens_avg_sig*10^9;

% filename = "C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\Pulsatility Analysis\variables\sub_6_"+j+"_10.csv";
% csvwrite(filename,ens_avg_sig);

%% Plotting a subplot for the abstract
close all;

sig_shift = 700;

j =1;
[ens_avg_sig,~,ensemble_curve] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_1lp(1,strt_time(j)*20:stp_time(j)*20),sig_shift,1);
a_1 = ensemble_curve;
j =4;
[ens_avg_sig,~,ensemble_curve] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_1lp(1,strt_time(j)*20:stp_time(j)*20),sig_shift,1);
a_2 = ensemble_curve;
j =1;
[ens_avg_sig,~,ensemble_curve] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_25lp(1,strt_time(j)*20:stp_time(j)*20),sig_shift,1);
a_3 = ensemble_curve;
j =4;
[ens_avg_sig,~,ensemble_curve] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_25lp(1,strt_time(j)*20:stp_time(j)*20),sig_shift,1);
a_4 = ensemble_curve;
close all;

hfig = figure();
sgtitle("Comparision of Pulsatility with and without Tourniquet")
x = (1:1:length(d))/1000;
ensavg = a_1(:,1:1000);
subplot(2,2,1);
plot(x, ensavg(2,:), '-r', 'LineWidth',1.5)
hold on;
plot(x, ensavg(1,:), ':g', 'LineWidth',1.5)
plot(x, ensavg(3,:), ':g', 'LineWidth',1.5)
xlabel("Time (s)");
ylabel("BFi_{1 cm}")
title(".   Flow Pulsatility Without Tourniquet (rsd = 1 cm)")
hold off
grid
% legend('Ensemble Average', '95% Confidence Intervals')


ensavg = a_2(:,1:1000);
subplot(2,2,2);
plot(x, ensavg(2,:), '-r', 'LineWidth',1.5)
hold on;
plot(x, ensavg(1,:), ':g', 'LineWidth',1.5)
plot(x, ensavg(3,:), ':g', 'LineWidth',1.5)
title("  Flow Pulsatility With Tourniquet (rsd = 1 cm)")
xlabel("Time (s)");
ylabel("Blood Flow index")
hold off
grid
% legend('Ensemble Average', '95% Confidence Intervals')


ensavg = a_3(:,1:1000);
subplot(2,2,3);
plot(x, ensavg(2,:), '-r', 'LineWidth',1.5)
hold on;
plot(x, ensavg(1,:), ':g', 'LineWidth',1.5)
plot(x, ensavg(3,:), ':g', 'LineWidth',1.5)
title("Flow Pulsatility Without Tourniquet (rsd = 2.5 cm)")
xlabel("Time (s)");
ylabel("Blood Flow index")
hold off
grid



ensavg = a_4(:,1:1000);
subplot(2,2,4);
plot(x, ensavg(2,:), '-r', 'LineWidth',1.5)
hold on;
plot(x, ensavg(1,:), ':g', 'LineWidth',1.5)
plot(x, ensavg(3,:), ':g', 'LineWidth',1.5)
title("Flow Pulsatility With Tourniquet (rsd = 2.5 cm)")
xlabel("Time (s)");
ylabel("Blood Flow index")
hold off
grid
% legend('Ensemble Average', '95% Confidence Intervals')


fname = "Ens_avg_figure_abstract"


% legend("Normalized 1cm","Norm 2.5cm","Norm 2.5cm TR","DCS 1cm","DCS 2.5cm","TR sys DCS 2.5cm ")

xlabel("Time(s)")
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-painters')
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document
legend('Ensemble Average', '95% Confidence Intervals','Interpreter','tex')



%% Plotting a subplot for the abstract
close all;

sig_shift = 300;

set(0, 'DefaultAxesFontname', 'Calisto MT');
font_size = 16;
set(0, 'DefaultAxesFontsize', font_size);
set(0, 'defaulttextfontname', 'Calisto MT');
set(0, 'defaulttextfontsize', font_size);
set(0, 'defaultfigurecolor', 'w')



j =1;
[ens_avg_sig,~,ensemble_curve] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_1lp(1,strt_time(j)*20:stp_time(j)*20),sig_shift ,1);
a_1 = ensemble_curve;
j =4;
[ens_avg_sig,~,ensemble_curve] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_1lp(1,strt_time(j)*20:stp_time(j)*20),sig_shift,1);
a_2 = ensemble_curve;
j =1;
[ens_avg_sig,~,ensemble_curve] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_25lp(1,strt_time(j)*20:stp_time(j)*20),sig_shift,1);
a_3 = ensemble_curve;
j =4;
[ens_avg_sig,~,ensemble_curve] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_25lp(1,strt_time(j)*20:stp_time(j)*20),sig_shift,1);
a_4 = ensemble_curve;
close all;


hfig =figure('units','centimeters', 'Position',[5 5 25 17]) %18 width 15 heigh\
% subaxis(1,1,1,'SpacingVert',0.05,'SpacingHoriz',0.08,'MR',0.05, 'ML',0.12,'MT',0.08,'MB',0.13)


% sgtitle("Comparision of Pulsatility with and without Tourniquet")
x = (1:1:length(d))/1000;
ensavg = a_1(:,1:1000);

subaxis(2,2,1,'SpacingVert',0.12,'SpacingHoriz',0.1,'MR',0.03, 'ML',0.08,'MT',0.07,'MB',0.11)

plot(x, ensavg(2,:), '-r', 'LineWidth',1.5)
hold on;
plot(x, ensavg(1,:), ':b', 'LineWidth',1.5)
plot(x, ensavg(3,:), ':b', 'LineWidth',1.5)
% xlabel("Time (s)");
ylabel("BFi_{1 cm} (cm^2/s)")
mn = min(ensavg(3,:));
mx = max(ensavg(1,:))
ylim([mn*0.93 mx*1.05])
% title(".   Flow Pulsatility Without Tourniquet (rsd = 1 cm)")
title("(a)")
hold off
grid
legend('Ensemble Average', '95% Confidence Intervals','FontSize',11)


ensavg = a_2(:,1:1000);
subaxis(2,2,2)
plot(x, ensavg(2,:), '-r', 'LineWidth',1.5)
hold on;
plot(x, ensavg(1,:), ':b', 'LineWidth',1.5)
plot(x, ensavg(3,:), ':b', 'LineWidth',1.5)
mn = min(ensavg(3,:));
mx = max(ensavg(1,:))
ylim([mn*0.93 mx*1.05])
% title("  Flow Pulsatility With Tourniquet (rsd = 1 cm)")
% xlabel("Time (s)");
% ylabel("Blood Flow index")
title("(b)")
hold off
grid
% legend('Ensemble Average', '95% Confidence Intervals')


ensavg = a_3(:,1:1000);
subaxis(2,2,3);
plot(x, ensavg(2,:), '-r', 'LineWidth',1.5)
hold on;
plot(x, ensavg(1,:), ':b', 'LineWidth',1.5)
plot(x, ensavg(3,:), ':b', 'LineWidth',1.5)
mn = min(ensavg(3,:));
mx = max(ensavg(1,:))
ylim([mn*0.93 mx*1.05])
% title("Flow Pulsatility Without Tourniquet (rsd = 2.5 cm)")
xlabel("Time (s)");
title("(c)")
ylabel("BFi_{2.5 cm} (cm^2/s)")

hold off
grid



ensavg = a_4(:,1:1000);
subaxis(2,2,4);
plot(x, ensavg(2,:), '-r', 'LineWidth',1.5)
hold on;
plot(x, ensavg(1,:), ':b', 'LineWidth',1.5)
plot(x, ensavg(3,:), ':b', 'LineWidth',1.5)
mn = min(ensavg(3,:));
mx = max(ensavg(1,:))
ylim([mn*0.93 mx*1.05])
% title("Flow Pulsatility With Tourniquet (rsd = 2.5 cm)")
xlabel("Time (s)");
title("(d)")
% ylabel("Blood Flow index")
hold off
grid
% legend('Ensemble Average', '95% Confidence Intervals')


fname = "Ens_avg_figure_abstract"


% legend("Normalized 1cm","Norm 2.5cm","Norm 2.5cm TR","DCS 1cm","DCS 2.5cm","TR sys DCS 2.5cm ")

% xlabel("Time(s)")
% picturewidth = 20; % set this parameter and keep it forever
% hw_ratio = 0.65; % feel free to play with this ratio
% set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(hfig,'Position');
% set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% %print(hfig,fname,'-dpdf','-painters','-fillpage')
% print(hfig,fname,'-dpng','-painters')
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document
% legend('Ensemble Average', '95% Confidence Intervals','Interpreter','tex')



%% Finding the s1, s2 and Diastolic features

d = ensemble_curve_dcs25_1(2,:);
hfig = figure();
t = (1:1:length(d))/1000;
[h_s,l_s] = findpeaks(ensemble_curve_dcs25_1(1:800),"MinPeakHeight",0,'NPeaks',5);
[h_m,l_m] = findpeaks(-d,"MinPeakHeight",1);
scatter(l_s/1000,h_s,'LineWidth',2.5,"MarkerFaceColor","r","MarkerEdgeColor","r");
text((l_s/1000)-0.008,h_s-0.1,["S1","S2","d"])
text((l_s/1000)-0.025,h_s-0.2,["(Systolic","(Reflective ","(Diacrotic"])
text((l_s/1000)-0.013,h_s-0.3,[" Peak)"," flow peak)"," notch)"])
hold on;
scatter(l_m/1000,-h_m,'LineWidth',2.5,"MarkerFaceColor","r","MarkerEdgeColor","r")
text((l_m/1000)-0.008,-h_m+0.1,["D", "(Diastolic Peak)"])

plot(t,d,'b','LineWidth',1.5);


fname = "Ens_avg_figure_abstract"


% legend("Normalized 1cm","Norm 2.5cm","Norm 2.5cm TR","DCS 1cm","DCS 2.5cm","TR sys DCS 2.5cm ")
title("Ensemble Average signal (rsd = 2.5 cm)")
xlabel("Time(s)")
ylabel("Normalized Amplitude")
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-painters')
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
text((l_s(2)/1000)-0.008,h_s(2)-0.6,["\Deltat = t_S_1 - t_d"],"Interpreter","tex")
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

%% Plotting the g2 averaged signal to find the deltaFc

%%
figure();
plot((1:1:size(adb_avg,2))/20, adb_avg(2,:),'b',"LineWidth",1.5);
title("Temporal Averaging Comparision");
xlabel("Time(s)")
ylabel("CBFi")
hold on; 
plot((1:1:size(adb_avg,2))/20,adb(2,:),'r');
legend("Ensemle Temporal Averaged","Raw")
%% Plotting the frequency spectrum of the averaged and raw signal
Fs = 20;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
% signal = adb_avg(2,:);
signal = aDb1(2,8000:11000);
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
title('FFT of raw DCS 2.5cm')
xlabel('f (Hz)')
ylabel('|P1(f)|')


%% Compare the 2.5cm and 1.5cm signals

t = (1:1:length(adb_avg(2,:)))/20;
plot(t,normalize(adb_avg(1,:)),'b'); hold on; plot(t,normalize(adb_avg(2,:)));
title("Comparison Between 1.5cm and 2.5cm DCS signal @TNQT 150%")
ylabel("Normalized Unit");
xlabel("Time(s)");
legend("r_s_d = 1.5cm","r_s_d = 2.5cm");

%% Compare the 2.5cm and 1 cm signals after filtering

dcs_1lp = lpf(adb_avg(1,:),7,20);
dcs_3lp = lpf(adb_avg(2,:),7,20);

t = (1:1:length(dcs_1lp))/20;
plot(t,normalize(dcs_1lp),'b'); hold on; plot(t,normalize(dcs_3lp));
title("Comparison: Filtered & g2 avg 1.5cm & 2.5cm DCS signal @TNQT 150%")
ylabel("Normalized Unit");
xlabel("Time(s)");
legend("r_s_d = 1.5cm","r_s_d = 2.5cm");

%% Compare the 2.5cm signal with and without g2 averaging
close all
t = (1:1:length(adb_avg(1,:)))/20;
plot(t,normalize(adb_avg(2,:)),'b'); hold on; plot(t,normalize(adb(2,:)));
title("Comparing 2.5cm signal with and without g2 averaging @TNQT 150%")
ylabel("Normalized Unit");
xlabel("Time(s)");
legend("g2 averaged","original");


%%
figure();
subplot(1,3,1);
semilogx(Data_tau,squeeze(Data(1,2,:)));
xlabel("Data Tau(s)")
title("Autocorrelation Curve");
subplot(1,3,2);
semilogx(Data_tau,squeeze(Data(20,2,:)));
xlabel("Data Tau(s)")
title("Autocorrelation Curve");
subplot(1,3,3);
semilogx(Data_tau,squeeze(Data(40,2,:)));
xlabel("Data Tau(s)")
title("Autocorrelation Curve");


%% Statistical analysis
an_table = zeros(2,5,3);
for i= 1:5
    an_table(:,i,1) = mean(squeeze(per_ch_g_pul(i,[1,3],3700:5500)),2);
    an_table(:,i,2) = mean(squeeze(per_ch_g_pul(i,[1,3],6100:7500)),2);
    an_table(:,i,3) = mean(squeeze(per_ch_g_pul(i,[1,3],8500:10000)),2);
end

csvwrite("anova_matrix_1cm.csv",squeeze(an_table(1,:,:)));
csvwrite("anova_matrix_25mm.csv",squeeze(an_table(2,:,:)));


%%
    g2(1,:,:) = Data(1:floor(length(ecg_ad)/50),1,:);
    g2(2,:,:) = Data(1:floor(length(ecg_ad)/50),2,:);
    g2(3,:,:) = Data(1:floor(length(ecg_ad)/50),3,:);
    g2(4,:,:) = Data(1:floor(length(ecg_ad)/50),4,:);
    
%     g2(:,1,:)=squeeze(D(:,1,:)-1); %g2-1 curve generation
%     g2_2_temp=squeeze(D(:,2,:)-1); %g2-1 curve generation
%     g2_3_temp=squeeze(D(:,3,:)-1); %g2-1 curve generation
%     g2_4_temp=squeeze(D(:,4,:)-1); %g2-1 curve generation
    
%     average g2 curve for large source detector separation
%     for i=1:size(g2,1)
%         g2(i,2,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
%     end
    
    [h_pks,l_pks] = findpeaks(normalize(ecg_ad),"MinPeakHeight",2.5,"MinPeakDistance",750);
    
    fig1=figure('units','centimeters', 'Position',[2 2 35 13]); %18 width 15 heigh
    hold on
    plot(normalize(ecg_ad),'r')
    plot(l_pks, h_pks,'*k')
    hd_pks = floor(h_pks./50);
    ld_pks = floor(l_pks./50);
    % g2_avg = zeros(size(g2,1),size(g2,2),size(g2,3));
    avg_window_width = 10;
    for i=1:size(ld_pks,2)-2
        if size(ld_pks,2)-i < avg_window_width
            avg_window_width = avg_window_width-1;
        end
        min_length = min(diff(ld_pks(i:i+avg_window_width)));
        base_sig = g2(:,ld_pks(i):ld_pks(i)+min_length,:);
       
        for j=1:avg_window_width-1
            base_sig = base_sig + g2(:,ld_pks(i+j):ld_pks(i+j)+min_length,:);
    %         adb_1 = hybrid_dcs(base_sig,Data_tau);
           
        end
        g2(:,ld_pks(i):ld_pks(i)+min_length,:) = base_sig./avg_window_width;
    end

    g2_n(:,1,:) = g2(1,1:floor(length(ecg_ad)/50),:);
    g2_n(:,2,:) = g2(2,1:floor(length(ecg_ad)/50),:);
    g2_n(:,3,:) = g2(3,1:floor(length(ecg_ad)/50),:);
    g2_n(:,4,:) = g2(4,1:floor(length(ecg_ad)/50),:);

    adb = standalone_dcs(g2_n,Data_tau);
%%
    close all;
    dcs_1lp = lpf(adb(1,:),7,20);
    dcs_25lp = lpf(adb(2,:),7,20);

    dcs_up = interp(dcs_1lp,50);
    dcs_up_lp = lpf(dcs_up,2,1500);
    bp_shift = circshift(bp_a,-2950);

    %%
    close all;
    bp_shift = circshift(bp_a,-1200);
    plot(normalize(ecg1)); hold on; 
    plot(normalize(bp_shift));
    %%
    close all;
    l = 186380:192750;
    t = (1:1:length(l))/1000;
    plot(t,dcs_up(l)*10^9);
    hold on;
    plot(t,bp_shift(l));
    xlabel("Time (s)");
    ylabel("aDb")

