close all;
tic;
file_list = ["20220218 - 8\","20220217\","20220217 - 6\","20220218 - 4\","20220223 - 2\","20220223 - 6\","20220304\","20220309 - 5\","20220309 - 7\"];
file_name = ["MPCM002","MPCM003","MPCM004","MPCM005","MPCM006","MPCM007","MPCM008","MPCM009","MPCM010"];

for j=1:length(file_list)

    % Read files;
    dcs_dir = strcat("C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\MPCM Study\DCS\",file_list(j),"Data.mat") ;     % Directory for DCS data
    
    load(dcs_dir);
    data_mat_all(j,:,:,:) = Data;
    aDb1 = hybrid_dcs(Data,Data_tau);
    
    
    % time resultion - aqusition time used to aquire data
    figure();
    sgtitle(file_name(j))
    t_res=1; % seconds
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
    
    plot(time,aDb1(3,:))
    title('{\itr}_{SD}=2 cm')
    % set(gca,'xticklabel',{})
    
    subplot(2,2,4)
    
    plot(time,aDb1(4,:))
    title('{\itr}_{SD}=2.5 cm')
    xlabel('Time (s)')
    
    % Assigning the channels
    dcs_1cm = aDb1(1,:).*10^9;
    dcs_1lp = lpf_ffilts(dcs_1cm,10,20);
    dcs_15 = aDb1(2,:).*10^9;
    dcs_15lp = lpf_ffilts(dcs_15,10,20);
    dcs_2 = aDb1(3,:).*10^9;
    dcs_2lp = lpf_ffilts(dcs_2,10,20);
    dcs_25 = aDb1(4,:).*10^9;
    dcs_25lp = lpf_ffilts(dcs_25,7,20);
    
    % Recombine the DCS signal
    adb_lp = [dcs_1lp;dcs_15lp;dcs_2lp;dcs_25lp];
    
    % Detrending the signal
    str_var = ["DCS 1cm","DCS 1.5cm","DCS 2cm","DSC 2.5cm"];
    if j==1
        per_change = zeros(length(file_list),size(adb_lp,1),size(adb_lp,2));
        snr_val = zeros(length(file_list),4,5);
    end
    for i=1:size(aDb1,1)
        s = adb_lp(i,:);
        s_str = str_var(i);
        breakpoints = 1:50:length(s);
        y = detrend(s,1,breakpoints);
    %     plot(y);
        
        % Plotting the mean average smoothing curve with the envelope
        x = (1:1:length(y))*0.05;
        [envHigh, envLow] = envelope(y,20,'peak');
        envMean = (envHigh+envLow)/2;
        env_diff = envHigh - envLow;
        figure();
        plot(x,y,x,envLow,x,envHigh,x,envMean);
        mh = max(envHigh);
        ml = max(envLow);
        mx = max(mh,ml);
%         
        snr_val(j,i,1) = db2mag(snr(s(1:800),20));
        snr_val(j,i,2) = db2mag(snr(s(800:1200),20));
        snr_val(j,i,3) = db2mag(snr(s(1200:4900),20));
        snr_val(j,i,4) = db2mag(snr(s(4900:5400),20));
        snr_val(j,i,5) = db2mag(snr(s(5400:6000),20));
        str = [strcat("SNR=",string(snr_val(j,i,1))),strcat("SNR=",string(snr_val(j,i,2))),strcat("SNR=",string(snr_val(j,i,3))),strcat("SNR=",string(snr_val(j,i,4))),strcat("SNR=",string(snr_val(j,i,5)))];
        
        xt = [10,40,135,240,265];
        yt = [mx,mx-1,mx-2,mx-1,mx];
        text(xt,yt,str,'Color','red','FontSize',14);
        text((max(x)/2)-10, mx,s_str, 'FontSize',16)
        
        
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
        per_change(j,i,:) = (100*(env_diff_smooth)./mean(env_diff_smooth(1:700))) ;
%         plot(x,squeeze(per_change(j,i,:)));
%         title("Percent change in pulse amplitude with pressure cuff " +s_str);
%         xlabel("Time (s)");
%         ylabel("% change in Pulse Amplitude");
    end
    % Plotting the changes in each channel in single plot
    figure();
    l = length(30:length(per_change)-10);
    x_pc = (1:1:l)/20;
    plot(x_pc,squeeze(per_change(j,1,30:end-10)),'DisplayName','per_change_dcs_1');
    hold on;
    plot(x_pc,squeeze(per_change(j,2,30:end-10)),'DisplayName','per_change_dcs_15');
    plot(x_pc,squeeze(per_change(j,3,30:end-10)),'DisplayName','per_change_dcs_2');
    plot(x_pc,squeeze(per_change(j,4,30:end-10)),'DisplayName','per_change_dcs_25');
    hold off;
    legend("DCS 1cm","DCS 1.5cm","DCS 2cm","DCS 2.5cm"); 
    title("Percent change in pulse amplitude with pressure cuff "+file_name(j));
    xlabel("Time (s)");
    ylabel("% change in Pulse Amplitude");
    break

end
toc;


%% Plotting the average g2 curves for each participants and saving them
close all;
clear all;
tic;
% file_list = ["20220218 - 8\","20220217\","20220217 - 6\","20220218 - 4\",...
%     "20220223 - 2\","20220223 - 6\","20220304\","20220309 - 5\","20220309
%     - 7\"]; % Files for LBNP with tnqt
file_list = ["20220218 - 6\","20220217 - 4\","20220217 - 8\","20220218 - 3\"] % File list for CCA
file_name = ["MPCM002","MPCM003","MPCM004","MPCM005","MPCM006","MPCM007","MPCM008","MPCM009","MPCM010"]';
g2_avg = zeros(length(file_list),3,4,50)*NaN;  % this dimension represents 1: #of files, 2:#of segments averaged,3:#of channels,4:tau points
for j=1:length(file_list)

    % Read files;
    dcs_dir = strcat("C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\MPCM Study\DCS\",file_list(j),"Data.mat") ;     % Directory for DCS data
    
    load(dcs_dir);
    
%     aDb1 = hybrid_dcs(Data,Data_tau);

    g2_1_temp=squeeze(Data(:,1,:)-1); %g2-1 curve generation
    g2_2_temp=squeeze(Data(:,2,:)-1); %g2-1 curve generation
    g2_3_temp=squeeze(Data(:,3,:)-1); %g2-1 curve generation
    g2_4_temp=squeeze(Data(:,4,:)-1); %g2-1 curve generation

    % average g2 curve for 4 channels
    g2(1,:,:)=(g2_1_temp);
    g2(2,:,:)=(g2_2_temp);
    g2(3,:,:)=(g2_3_temp); 
    g2(4,:,:)=(g2_4_temp);

    g2_avg(j,1,:,:) = squeeze(mean(g2(:,1:600,:),2));
    g2_avg(j,2,:,:) = squeeze(mean(g2(:,600:1200,:),2));
    g2_avg(j,3,:,:) = squeeze(mean(g2(:,1800:2400,:),2));
    ttle_str = ["0 mmHg","50 mmHg","180 mmHg"];
    figure();
%     for s=1:size(g2_avg,2)
%         subplot(3,1,s);
%         semilogx(Data_tau,squeeze(g2_avg(j,s,1,:)),'r')
%         hold on
%         semilogx(Data_tau,squeeze(g2_avg(j,s,2,:)),'g')
%         semilogx(Data_tau,squeeze(g2_avg(j,s,3,:)),'b')
%         semilogx(Data_tau,squeeze(g2_avg(j,s,4,:)),'magenta')
%         ylabel('g_2')
%         xlabel("Tau values")
%         legend('1 cm','1.5 cm','2 cm','2.5 cm')
%         title("Average g2 curve for "+ttle_str(s)+" tnqt " +file_name(j))
%     end
    semilogx(Data_tau,squeeze(g2_avg(j,1,1,:)),Color='r',LineWidth=0.9);
    legend("1 cm")
    hold on;
    semilogx(Data_tau,squeeze(g2_avg(j,2,1,:)),'--',Color='r',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,3,1,:)),'-.',Color='r',LineWidth=0.9);
    ylabel('g_2')
    xlabel("tau values")

    semilogx(Data_tau,squeeze(g2_avg(j,1,4,:)),Color='b',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,2,4,:)),'--',Color='b',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,3,4,:)),'-.',Color='b',LineWidth=0.9);
    legend("1cm - 0 mmHg","1 cm - 50 mmHg",'1 cm - 180 mmHg',"2.5cm - 0 mmHg","2.5cm - 50 mmHg",'2.5cm - 180 mmHg');
    title("Average g2 curve for tnqt " +file_name(j))

%     break

end
filename = "avg_g2_MPCM.mat";
save(filename,"g2_avg");

for j = 1:4:size(g2_avg,1)-1
    figure();
    subplot(2,2,1)
    semilogx(Data_tau,squeeze(g2_avg(j,1,1,:)),Color='r',LineWidth=0.9);
    legend("1 cm")
    hold on;
    semilogx(Data_tau,squeeze(g2_avg(j,2,1,:)),'--',Color='r',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,3,1,:)),'-.',Color='r',LineWidth=0.9);
    ylabel('g_2')
    xlabel("tau values")

    semilogx(Data_tau,squeeze(g2_avg(j,1,4,:)),Color='b',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,2,4,:)),'--',Color='b',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,3,4,:)),'-.',Color='b',LineWidth=0.9);
    legend("1cm - 0 mmHg","1 cm - 50 mmHg",'1 cm - 180 mmHg',"2.5cm - 0 mmHg","2.5cm - 50 mmHg",'2.5cm - 180 mmHg');
    title("Average g2 curve for tnqt " +file_name(j))
    
    j=j+1;

    subplot(2,2,2)
    semilogx(Data_tau,squeeze(g2_avg(j,1,1,:)),Color='r',LineWidth=0.9);
    legend("1 cm")
    hold on;
    semilogx(Data_tau,squeeze(g2_avg(j,2,1,:)),'--',Color='r',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,3,1,:)),'-.',Color='r',LineWidth=0.9);
    ylabel('g_2')
    xlabel("tau values")

    semilogx(Data_tau,squeeze(g2_avg(j,1,4,:)),Color='b',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,2,4,:)),'--',Color='b',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,3,4,:)),'-.',Color='b',LineWidth=0.9);
    legend("1cm - 0 mmHg","1 cm - 50 mmHg",'1 cm - 180 mmHg',"2.5cm - 0 mmHg","2.5cm - 50 mmHg",'2.5cm - 180 mmHg');
    title("Average g2 curve for tnqt " +file_name(j))
    
    j=j+1;
    subplot(2,2,3)
    semilogx(Data_tau,squeeze(g2_avg(j,1,1,:)),Color='r',LineWidth=0.9);
    legend("1 cm")
    hold on;
    semilogx(Data_tau,squeeze(g2_avg(j,2,1,:)),'--',Color='r',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,3,1,:)),'-.',Color='r',LineWidth=0.9);
    ylabel('g_2')
    xlabel("tau values")

    semilogx(Data_tau,squeeze(g2_avg(j,1,4,:)),Color='b',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,2,4,:)),'--',Color='b',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,3,4,:)),'-.',Color='b',LineWidth=0.9);
    legend("1cm - 0 mmHg","1 cm - 50 mmHg",'1 cm - 180 mmHg',"2.5cm - 0 mmHg","2.5cm - 50 mmHg",'2.5cm - 180 mmHg');
    title("Average g2 curve for tnqt " +file_name(j))
    
    j=j+1;
    subplot(2,2,4)
    semilogx(Data_tau,squeeze(g2_avg(j,1,1,:)),Color='r',LineWidth=0.9);
    legend("1 cm")
    hold on;
    semilogx(Data_tau,squeeze(g2_avg(j,2,1,:)),'--',Color='r',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,3,1,:)),'-.',Color='r',LineWidth=0.9);
    ylabel('g_2')
    xlabel("tau values")

    semilogx(Data_tau,squeeze(g2_avg(j,1,4,:)),Color='b',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,2,4,:)),'--',Color='b',LineWidth=0.9);
    semilogx(Data_tau,squeeze(g2_avg(j,3,4,:)),'-.',Color='b',LineWidth=0.9);
    legend("1cm - 0 mmHg","1 cm - 50 mmHg",'1 cm - 180 mmHg',"2.5cm - 0 mmHg","2.5cm - 50 mmHg",'2.5cm - 180 mmHg');
    title("Average g2 curve for tnqt " +file_name(j))


end

%% Back up

close all;
tic;
file_list = ["20220218 - 8\","20220217\","20220217 - 6\","20220218 - 4\",...
    "20220223 - 2\","20220223 - 6\","20220304\","20220309 - 5\","20220309 - 7\"]; % Files for LBNP with tnqt
% file_list = ["20220218 - 6\","20220217 - 4\","20220217 - 8\","20220218 - 3\"] % File list for CCA
% file_name = ["MPCM002","MPCM003","MPCM004","MPCM005","MPCM006","MPCM007","MPCM008","MPCM009","MPCM010"];
perchange_val = zeros(length(file_list),4)*NaN;
for j=1:length(file_list)

    % Read files;
    dcs_dir = strcat("C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\MPCM Study\DCS\",file_list(j),"Data.mat") ;     % Directory for DCS data
    
    load(dcs_dir);
    
    aDb1 = hybrid_dcs(Data,Data_tau);
    
    % Plots the fitted g2 curve
    if exist("g2")
        clear g2
    end
    mua = 0.17; %cm^-1 baseline absorption coefficient
    mus = 10; 

    g2(1,:,:)=squeeze(Data(:,1,:)-1); %g2-1 curve generation
    g2(2,:,:)=squeeze(Data(:,2,:)-1); %g2-1 curve generation
    g2(3,:,:)=squeeze(Data(:,3,:)-1); %g2-1 curve generation
    g2(4,:,:)=squeeze(Data(:,4,:)-1); %g2-1 curve generation

    Channel=1;
    Curve_no=1500;
    rho = [1 1.5 2 2.5];

    beta=g2(Channel,Curve_no,1);
    aDb=aDb1(Channel,Curve_no);
    
    g2_fit=gen_DCS_fit(Data_tau,mua,mus,rho(Channel),beta,aDb);

    semilogx(Data_tau,squeeze(Data(Curve_no, Channel,:))-1,'k')
    hold on
    semilogx(Data_tau,g2_fit,'r')
    %
    
    
    % Assigning the channels
    dcs_1cm = aDb1(1,:).*10^9;
    dcs_1lp = lpf_ffilts(dcs_1cm,10,20);
    dcs_15 = aDb1(2,:).*10^9;
    dcs_15lp = lpf_ffilts(dcs_15,10,20);
    dcs_2 = aDb1(3,:).*10^9;
    dcs_2lp = lpf_ffilts(dcs_2,10,20);
    dcs_25 = aDb1(4,:).*10^9;
    dcs_25lp = lpf_ffilts(dcs_25,7,20);
    
    % Recombine the DCS signal
    adb_lp = [dcs_1lp;dcs_15lp;dcs_2lp;dcs_25lp];
    breakpoints = 1:50:length(dcs_1cm);
    for k=1:size(adb_lp,1)
        adb_lp(k,:) = detrend(adb_lp(k,:),1,breakpoints);
    end
    
    % Detrending the signal
    str_var = ["DCS 1cm","DCS 1.5cm","DCS 2cm","DSC 2.5cm"];
    if j==1
        per_change = zeros(length(file_list),size(adb_lp,1),size(adb_lp,2));
        snr_val = zeros(length(file_list),4,5);
    end
    envHigh = zeros(size(aDb1,1),size(aDb1,2));
    envLow = zeros(size(aDb1,1),size(aDb1,2));
    envMean = zeros(size(aDb1,1),size(aDb1,2));
    env_diff = zeros(size(aDb1,1),size(aDb1,2));
    env_diff_smooth = zeros(size(aDb1,1),size(aDb1,2));

    avg_sig = zeros(size(file_name,2),size(aDb1,1),size(aDb1,2));
    for i=1:size(aDb1,1)
        s = adb_lp(i,:);
        s_str = str_var(i);
        
        y = s;
    %     plot(y);
        
        % Plotting the mean average smoothing curve with the envelope
        x = (1:1:length(y))*0.05;
        [envHigh(i,:), envLow(i,:)] = envelope(y,20,'peak');
        envMean(i,:) = (envHigh(i,:)+envLow(i,:))/2;
        env_diff(i,:) = envHigh(i,:) - envLow(i,:);
%         figure();
%         plot(x,y,x,envLow(i),x,envHigh(i),x,envMean(i));
%         mh = max(envHigh(i));
%         ml = max(envLow(i));
%         mx = max(mh,ml);
% %         
        snr_val(j,i,1) = db2mag(snr(s(1:600),20));
        snr_val(j,i,2) = db2mag(snr(s(600:1200),20));
        snr_val(j,i,3) = db2mag(snr(s(1300:4900),20));
        snr_val(j,i,4) = db2mag(snr(s(4900:5400),20));
        snr_val(j,i,5) = db2mag(snr(s(5400:6000),20));

        env_diff_smooth(i,:) = smooth(env_diff(i,:),400,'rlowess');

        avg_sig(j,i,:) = smooth(s,'sgolay');
        per_change(j,i,:) = (100*(env_diff_smooth(i,:))./mean(env_diff_smooth(i,1:600))) ;
        
    end

    % Plotting the signal with SNR
    figure();
    sgtitle(file_name(j))
    t_res=1; % seconds
    time=t_res*(1:1:size(aDb1,2));
    
    subplot(2,2,1)
    i = 1;
    plot(x,adb_lp(i,:),x,envHigh(i,:),x,envMean(i,:),x,envLow(i,:));
    mh = max(envHigh(i,:));
    ml = max(envLow(i,:));
    mx = max(mh,ml);
    str = [strcat("SNR=",string(snr_val(j,i,1))),strcat("SNR=",string(snr_val(j,i,2))),...
        strcat("SNR=",string(snr_val(j,i,3))),strcat("SNR=",string(snr_val(j,i,4))),strcat("SNR=",string(snr_val(j,i,5)))];
    title('{\itr}_{SD}=1 cm')
    xt = [10,40,135,240,265];
    yt = [mx,mx-1,mx-2,mx-1,mx];
    text(xt,yt,str,'Color','red','FontSize',10,"color",'k');
    text((max(x)/2)-10, mx,str_var(i), 'FontSize',12,"color",'b');
    xlabel('Time (s)')
    ylabel("aDb value");
%     set(gca,'xticklabel',{})
    
    subplot(2,2,2)
    i = 2;
    plot(x,adb_lp(i,:),x,envHigh(i,:),x,envMean(i,:),x,envLow(i,:));
    mh = max(envHigh(i,:));
    ml = max(envLow(i,:));
    mx = max(mh,ml);
    str = [strcat("SNR=",string(snr_val(j,i,1))),strcat("SNR=",string(snr_val(j,i,2))),...
        strcat("SNR=",string(snr_val(j,i,3))),strcat("SNR=",string(snr_val(j,i,4))),strcat("SNR=",string(snr_val(j,i,5)))];
    title('{\itr}_{SD}=1.5 cm')
    xt = [10,40,135,240,265];
    yt = [mx,mx-1,mx-2,mx-1,mx];
    text(xt,yt,str,'Color','red','FontSize',10,"color",'k');
    text((max(x)/2)-10, mx,str_var(i), 'FontSize',12,"color",'b');
    xlabel('Time (s)')
    ylabel("aDb value");
    
    subplot(2,2,3)
    
    i = 3;
    plot(x,adb_lp(i,:),x,envHigh(i,:),x,envMean(i,:),x,envLow(i,:));
    mh = max(envHigh(i,:));
    ml = max(envLow(i,:));
    mx = max(mh,ml);
    str = [strcat("SNR=",string(snr_val(j,i,1))),strcat("SNR=",string(snr_val(j,i,2))),...
        strcat("SNR=",string(snr_val(j,i,3))),strcat("SNR=",string(snr_val(j,i,4))),strcat("SNR=",string(snr_val(j,i,5)))];
    title('{\itr}_{SD}=2 cm')
    xt = [10,40,135,240,265];
    yt = [mx,mx-1,mx-2,mx-1,mx];
    text(xt,yt,str,'Color','red','FontSize',10,"color",'k');
    text((max(x)/2)-10, mx,str_var(i), 'FontSize',12,"color",'b');
    xlabel('Time (s)')
    ylabel("aDb value");
    % set(gca,'xticklabel',{})
    
    subplot(2,2,4)
    
    i = 4;
    plot(x,adb_lp(i,:),x,envHigh(i,:),x,envMean(i,:),x,envLow(i,:));
    mh = max(envHigh(i,:));
    ml = max(envLow(i,:));
    mx = max(mh,ml);
    str = [strcat("SNR=",string(snr_val(j,i,1))),strcat("SNR=",string(snr_val(j,i,2))),...
        strcat("SNR=",string(snr_val(j,i,3))),strcat("SNR=",string(snr_val(j,i,4))),strcat("SNR=",string(snr_val(j,i,5)))];
    title('{\itr}_{SD}=2.5 cm')
    xt = [10,40,135,240,265];
    yt = [mx,mx-1,mx-2,mx-1,mx];
    text(xt,yt,str,'Color','red','FontSize',10,"color",'k');
    text((max(x)/2)-10, mx,str_var(i), 'FontSize',12,"color",'b');
    xlabel('Time (s)')
    ylabel("aDb value");

    % Plotting the changes in each channel in single plot
    figure();
    l = length(30:length(per_change)-30);
    x_pc = (1:1:l)/20;
    plot(x_pc,squeeze(per_change(j,1,30:end-30)),'DisplayName','per_change_dcs_1cm ');
    hold on;
    plot(x_pc,squeeze(per_change(j,2,30:end-30)),'DisplayName','per_change_dcs_1.5cm');
    plot(x_pc,squeeze(per_change(j,3,30:end-30)),'DisplayName','per_change_dcs_2cm');
    plot(x_pc,squeeze(per_change(j,4,30:end-30)),'DisplayName','per_change_dcs_2.5cm');
    hold off;
    legend("DCS 1cm","DCS 1.5cm","DCS 2cm","DCS 2.5cm"); 
    title("Percent change in pulse amplitude with pressure cuff "+file_name(j));
    xlabel("Time (s)");
    ylabel("% change in Pulse Amplitude");

    % Calculating the percent change in the signal compared to baseline 30s
    % with 30s data after cuff inflation to 180 mmHg

    perchange_val(j,1) = mean(squeeze(per_change(j,1,1:600)))- mean(squeeze(per_change(j,1,1800:2400)));
    perchange_val(j,2) = mean(squeeze(per_change(j,2,1:600)))- mean(squeeze(per_change(j,2,1800:2400)));
    perchange_val(j,3) = mean(squeeze(per_change(j,3,1:600)))- mean(squeeze(per_change(j,3,1800:2400)));
    perchange_val(j,4) = mean(squeeze(per_change(j,4,1:600)))- mean(squeeze(per_change(j,4,1800:2400)));


%     break

end
toc;
