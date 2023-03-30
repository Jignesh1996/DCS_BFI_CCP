clear all
close all
%%
lambda = 786e-7; % Wavelenght in cm
% tau=load('data_tau.csv');
beta=0.5;

% tau=1e-6:10e-6:5e-3;

%% 2 layered model
close all;
beta=0.5;
LayNo=2;
rho=0.7;
rho2=3;
% l=[0.4:0.1:3.5]; % First layer thicknes in cm
% val = [1:0.1:3.5];
val = [0.6:0.1:2.5];
l = val;
% change = [0.1:0.1:1];
fact= zeros(1,length(val));

% for k = 1:1:length(change)
    for j=1:length(fact)
        
        mua=[0.12 0.16];
        mus=[10 6];
        L = l(j);
       % plot the sensitivity factor wrt the different reduction in the flow
       % values vs perticular L values. so @ fixed L > flow vs sensitivity
       % factor again with different L values.
        F=[1.4e-9 1.4e-8];
        beta1 = 0.45;
        beta2 = 0.1568;
        % F = [mean(aDb1(1,1:600)) mean(aDb1(2,1:600))];
        
        [G1(1,:),g2(1,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta1);
        [G1(2,:),g2(3,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta2);
        
        L = l(j);
        F=[0.7e-9 1.4e-8];
%         F=[(4e-8)*change(k) 2e-9];
        % F = [mean(aDb1(1,600:900)),mean(aDb1(2,600:900))];
        [G1(3,:),g2(2,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta1);
        [G1(4,:),g2(4,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta2);
    
        rsd=[rho rho2];
        mua = 0.16; %cm^-1 baseline absorption coefficient
        mus = 8; %cm^-1 baseline reduced scattering coefficient
        
        for i=1:4
            if i<3
                rsd=rho;
                beta = 0.45;
            else
                rsd=rho2;
                beta = 0.1568;
            end    
                LB = [0.1e-15];
                UB = [10e-6];
            
            Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
%             beta= 0.1568; %0.1568;
            options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
            [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau,g2(i,:)-1,mua,mus,rsd,beta);
            aDb(j,i) = FittedParams(1);
        end
    
        g2_temp=round(g2,4);
    
        DeltaOD_sp=-log((g2_temp(2,:)-1)./(g2_temp(1,:)-1));
        DeltaOD_lp=-log((g2_temp(4,:)-1)./(g2_temp(3,:)-1));
        
        
        DeltaFec(1,:)=DeltaOD_sp/(0.3*10^-9);
        DeltaFec(2,:)=DeltaOD_lp/(0.3*10^-9);
        % 
        % DeltaOD_s = -log((g2_sig(1,:)-1)/(g2(1,:)-1));
        % DeltaOD_l = -log((g2_sig(2,:)-1)/(g2(3,:)-1));
        
        ratio = DeltaFec(2,:)./DeltaFec(1,:);
    %     semilogx(tau,DeltaFec(2,:)./DeltaFec(1,:));
        
        index=ratio==Inf;
        
        
        sens_fact = nanmean(ratio(index==0))
        fact(j) =sens_fact;
        
    
    
    end
    plt(l,fact,"Thickness Vs Sensitivity Factor@ flow reduction","L (Thickness)","Sensitivity Factor","");
%     plot(l,fact,'DisplayName',string(100-(100*change(k)))+'%')
%     legend('-DynamicLegend');
%     xlabel("Length");
%     ylabel("Sensitivity Factor")
%     title("Thickness vs S.F.")
% %     legend(string(100/k))
%     hold all;
% end
%% Calculating the error due to the false thickness assumption
bsl = 1.2;
bsl_fact = fact(3);
for i = 1:1:length(fact)
    error(i) = abs(bsl_fact-fact(i));
end
plt(l,error,"Error","Thickness","Error","");

%% DeltaFc
% deltaFc =  
%% 2 layered model / Change in the SF wrt FLow Reduction and the thickness
close all;
beta=0.5;
LayNo=3;
rho=1;
rho2=2.5;
% l=[0.4:0.1:3.5]; % First layer thicknes in cm
val = [0.8:0.1:3];
l = val;
red = [0.1:0.1:1]; 
fact= zeros(1,length(val));

for k=1:length(red)
    for j=1:length(fact)
        mua=[0.12 0.16];
        mus=[10 6];
        L = l(j);
       % plot the sensitivity factor wrt the different reduction in the flow
       % values vs perticular L values. so @ fixed L > flow vs sensitivity
       % factor again with different L values.
        F=[1.4e-9 1.4e-8];
        beta1 = 0.45;
        beta2 = 0.1568;
        % F = [mean(aDb1(1,1:600)) mean(aDb1(2,1:600))];
        
        [G1(1,:),g2(1,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta1);
        [G1(2,:),g2(3,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta2);
        
        L = l(j);
        F=[(1-red(k))*1.4e-9 1.4e-8];
%         F=[(4e-8)*change(k) 2e-9];
        % F = [mean(aDb1(1,600:900)),mean(aDb1(2,600:900))];
        [G1(3,:),g2(2,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta1);
        [G1(4,:),g2(4,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta2);
    
        rsd=[rho rho2];
        mua = 0.16; %cm^-1 baseline absorption coefficient
        mus = 8; %cm^-1 baseline reduced scattering coefficient
        
        for i=1:4
            if i<3
                rsd=rho;
            else
                rsd=rho2;
            end    
                LB = [0.1e-15];
                UB = [10e-6];
            
            Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
            beta= 0.1568; %0.1568;
            options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
            [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau,g2(i,:)-1,mua,mus,rsd,beta);
            aDb(i) = FittedParams(1);
        end
    
        g2_temp=round(g2,4);
    
        DeltaOD_sp=-log((g2_temp(2,:)-1)./(g2_temp(1,:)-1));
        DeltaOD_lp=-log((g2_temp(4,:)-1)./(g2_temp(3,:)-1));
        
        
        DeltaFec(1,:)=DeltaOD_sp/(0.3*10^-9);
        DeltaFec(2,:)=DeltaOD_lp/(0.3*10^-9);
        % 
        % DeltaOD_s = -log((g2_sig(1,:)-1)/(g2(1,:)-1));
        % DeltaOD_l = -log((g2_sig(2,:)-1)/(g2(3,:)-1));
        
        ratio = DeltaFec(2,:)./DeltaFec(1,:);
    %     semilogx(tau,DeltaFec(2,:)./DeltaFec(1,:));
        
        index=ratio==Inf;
        
        
        sens_fact = nanmean(ratio(index==0))
        fact(k,j) =sens_fact;
        
    
    
    end
    x(k) = (1-red(k))*100;
%     plt(x,fact,"Flow reduction Vs Sensitivity Factor","Flow reduction (percent)","Sensitivity Factor",string(1-red(k)*100));
%     hold all;

end
    plot(val,fact)
    legend(string(100-x))


%%
beta=0.5;
LayNo=2;
rho=0.7;
rho2=3;
mua=[0.12 0.16];
mus=[10 6];
L = 1.2;
beta1 = 0.45;
beta2 = 0.1568;

% plot the sensitivity factor wrt the different reduction in the flow
% values vs perticular L values. so @ fixed L > flow vs sensitivity
% factor again with different L values.
F=[1.4e-9 1.4e-8];

% F = [mean(aDb1(1,1:600)) mean(aDb1(2,1:600))];

[G1(1,:),g2(1,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta1);
[G1(4,:),g2(4,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta2);

%         L = l(j);
F=[0.7e-9 1.4e-8];
% F = [mean(aDb1(1,600:900)),mean(aDb1(2,600:900))];
[G1(2,:),g2(2,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta1);
[G1(5,:),g2(5,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta2);

F=[1.6e-9 1.6e-8];
% F = [mean(aDb1(1,600:900)),mean(aDb1(2,600:900))];
[G1(3,:),g2(3,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta1);
[G1(6,:),g2(6,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta2);
%%
figure()

semilogx(tau,g2(1,:),'r')
hold on
semilogx(tau,g2(2,:),'--r')
semilogx(tau,g2(4,:),'k')
semilogx(tau,g2(5,:),'--k')
semilogx(tau,g2(3,:),'b')
semilogx(tau,g2(6,:),'--b')
ylabel('g_2')
xlabel('Tau (s)')
legend('1cm','1cm_{red}', '2.5cm','2.5cm_{red}','Location','SouthWest')

%%
rsd=[0.7 3];
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 8; %cm^-1 baseline reduced scattering coefficient

for i=1:6
    if i<4 % <4 because ch 1,2,and 3 are for short rsd and 4,5,6 are long rsd
        rsd=rho;
        beta = 0.45;
    else
        rsd=rho2;
        beta = 0.1568;
    end    
        LB = [0.1e-15];
        UB = [10e-6];
        
        Starting = [1e-10]; %[aDb, Beta; cm^2/s, a.u.]
%         beta= 0.1568; %0.1568;
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau,g2(i,:)-1,mua,mus,rsd,beta);
        aDb(i) = FittedParams(1);
end
%%

figure()

bar(1,100*(aDb(2)/aDb(1)-1),'r')
hold on
bar(2,100*(aDb(5)/aDb(4)-1),'k')
ylabel('\DeltaBFi (%)')
xlabel('Tau (s)')
% legend('1cm','1cm_{red}', '2.5cm','2.5cm_{red}','Location','SouthWest')


%%
percent_ch(1) = 100-(mean(aDb1(1,600:900))/mean(aDb1(1,1:600)))*100;
percent_ch(2) = 100-(mean(aDb1(2,600:900))/mean(aDb1(2,1:600)))*100;


%%
g2_temp=round(g2,4);

DeltaOD_sp=-log((g2_temp(2,:)-1)./(g2_temp(1,:)-1));
DeltaOD_lp=-log((g2_temp(5,:)-1)./(g2_temp(4,:)-1));
DeltaOD_s = -log((g2_temp(3,:)-1)./(g2_temp(1,:)-1));
DeltaOD_l = -log((g2_temp(6,:)-1)./(g2_temp(4,:)-1));

DeltaFec(1,:)=DeltaOD_sp/(0.3*10^-9);
DeltaFec(2,:)=DeltaOD_lp/(0.3*10^-9);
% 
% DeltaOD_s = -log((g2_sig(1,:)-1)/(g2(1,:)-1));
% DeltaOD_l = -log((g2_sig(2,:)-1)/(g2(3,:)-1));

ratio = DeltaFec(2,:)./DeltaFec(1,:);
semilogx(tau,DeltaFec(2,:)./DeltaFec(1,:));

index=ratio==Inf;


sens_fact = nanmean(ratio(index==0))

figure()

semilogx(tau,DeltaFec(1,:),'r')
hold on
semilogx(tau,DeltaFec(2,:),'b')

%% Calculating the DeltaFc
DeltaFc = abs(DeltaOD_l) - (abs(DeltaOD_lp./DeltaOD_sp).*abs(DeltaOD_s))

g2_new = g2(6,:);
g2_new(2:35) = g2(6,2:35).*DeltaFc(2:35);
semilogx(tau, g2_new)

%%
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10;
Channel=1;
Curve_no=1;
rho = [1.7 2.7];

beta=g2(Channel, Curve_no,1);
aDb1=aDb1(Channel,Curve_no);

g2_fit=gen_DCS_fit(tau_values,mua,mus,rho(Channel),beta,aDb1);

semilogx(tau_values,squeeze(g2(1,1,:)),'k')
hold on
semilogx(tau_values,g2_fit,'r')

%%
% Checking for the pressure modulation data to make sure that it worked
adb = zeros(22,4,1800);
for i = 1:22
    if i<10
        file_name = strcat("MCI00",string(i),"_PRESSURE_Data.mat");
    elseif i==13
        continue
    else
        file_name = strcat("MCI0",string(i),"_PRESSURE_Data.mat");
    end
    dir = "D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\Parkwood Study\MATfiles\";
    file_path  = strcat(dir,file_name);
    load(file_path);
    aDb = hybrid_dcs(Data(1:1800,:,:),Data_tau);
    adb(i,:,:) = aDb;

end

%% Calculating the SF for real data

clear all;
close all

file_no = 6;
%--------------------------------------------------------------------------
% 
% filename_d=strcat('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\DCS\20220608 - '+string(file_no),'\Data.mat');
% load(filename_d)
% 
% 
% %Loading the ECG file
% filename_nd=strcat('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\ECG\20220608_'+string(file_no)+'.mat');
% load(filename_nd)

filename_d=strcat('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\new pressure cuff\Hypercapnia\Farah\Baseline_HCP_test\Data.mat');
load(filename_d)


%Loading the ECG file
filename_nd=strcat('D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\new pressure cuff\Hypercapnia\Farah\Baseline_HCP_test\Baseline_HCP_Test.mat');
load(filename_nd)


% load("D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Matlab_Code_aDb\2 layer model\2 layer model\008_CUFF_g2.mat");
% tau = csvread("D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Matlab_Code_aDb\2 layer model\data_tau.csv");
% Data_tau = tau;
% load("D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Codes\Matlab_Code_aDb\2 layer model\001_TCD_g2.mat");


if ~exist('Data') && exist("g2")
    Data(:,1,:) = g2(1,:,:)+1;
    Data(:,2,:) = g2(2,:,:)+1;
    Data(:,3,:) = g2(3,:,:)+1;
    Data(:,4,:) = g2(4,:,:)+1;
end
tau = Data_tau;
aDb1 = standalone_dcs(Data,Data_tau);
figure();
plot((1:1:length(aDb1))/20,aDb1');
xlabel("Time(s)")
ylabel("aDb")

%% Pulsatility analysis and g2 averaging of the signal

l = 60; % length of signal as ROI in seconds
s = 1; % Signal ROI starting point in seconds;

ecg_a = data(datastart(1):dataend(1));
% bp_a = data(datastart(2):dataend(2));

% ecg1 = ecg_a(1:660000);
ecg1 = ecg_a(1:6000000);
ecg1 = normalize(ecg1);
ecg1 = lpf(ecg1,5,1000);

% bp_a = bp_a(1:length(ecg1));
% bp_a = lpf(bp_a,3,1000);

%%
% Calculating FT frequency spectrum
close all;
Fs = 20;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
signal = aDb1(2,1:1200);
L = length(signal);             % Length of signal
t = (0:L-1)*T;  

Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot((1:1:length(P1))/Fs,P1)

%%
ecg1 = downsample(ecg1,10);
%% Plotting the ensemble average graph
close all;
shift =450;
[ens_avg_sig,pind] = ensemble_avg(ecg1(1:60000),aDb1(2,1:1200),shift,1);
%% Finding the shift

close all

hfig = figure();
picturewidth = 30; % set this parameter and keep it forever
hw_ratio = 0.5; % feel free to play with this ratio

ecg_ad = circshift(ecg1,-shift);
x_e = (1:1:length(ecg_ad(2900:5750)))/1000;
plot(x_e,normalize(ecg_ad(2900:5750)),'b');
hold on;
x_d = (1:1:size(aDb1(2,58:115),2))/20;
plot(x_d,normalize(aDb1(1,58:115)),'r');
xlabel("Time(s)")
ylabel("A.U.")
legend("ECG","BFi")




set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%% Finding the R-R peaks of ECG signal

close all;
% if exist("g2")
%     clear g2;
% end
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10; 

    
% 
g2_backup = Data;
% ecg = ecg_ad;
% g2_backup = g2;
g2_n = g2_backup;
% strt_time =  [1,180,300,420,540];   %time in seconds
% stp_time =  [180,300,420,540,660];    %time in seconds.

 strt_time =  [1,60,120,240,540];
 stp_time =  [60,120,240,540,600]; 

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
    
    
    [h_pks,l_pks] = findpeaks(normalize(ecg),"MinPeakHeight",2.5,"MinPeakDistance",750);
    
    fig1=figure('units','centimeters', 'Position',[2 2 35 13]); %18 width 15 heigh
    hold on
    plot(normalize(ecg),'r')
    plot(l_pks, h_pks,'*k')
    hd_pks = ceil(h_pks./50);
    ld_pks = ceil(l_pks./50);
    % g2_avg = zeros(size(g2,1),size(g2,2),size(g2,3));
    avg_window_width = 20;
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
figure()


% adb_avg = adb_avg; % This is cut out specific portion of the signal to plot
% figure();
% adb = standalone_tr_dcs(Data_all,Data_tau);
dcs_1lp = (lpf(adb_avg(1,:),7,20));
dcs_25lp =(lpf(adb_avg(2,:),7,20));
% snr(dcs_25lp(8400:10800),20)
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

% %% Finding the outliers from the singal and replacing them withe mean of
%  % adjacent values..
%  data = g2_n;
% close all;
% otlr_1 = isoutlier(adb_avg(1,:),"movmedian",[3 3]);
% t = 1:1:length(otlr_1);
% plot(t,adb_avg(1,:));
% hold on;
% ind = find(otlr_1);
% scatter(ind,adb_avg(1,otlr_1));
% hold off;
% % Removing the outliers from the g2 curves and refitting them
% for i=2:1:length(ind)
%     data(ind(i),:,:) = mean(data(ind(i)-3:ind(i)+3,:,:),1);
% end
% adb_new = standalone_dcs(data,Data_tau);
% figure()
% plot(adb_avg(1,:),'b');
% hold on;
% plot(adb_new(1,:),'r')
% 
% dcs_1lp = (lpf(adb_new(1,:),7,20));
% dcs_25lp =(lpf(adb_new(2,:),7,20));


%% Making the figure for the g2 average illustration
% close all;
% hfig = figure();
% picturewidth = 30; % set this parameter and keep it forever
% hw_ratio = 0.5; % feel free to play with this ratio
% % sgtitle();
% for i=3:6
%     ylim([0.98 1.165])
%     subplot(1,4,i-2)
%     semilogx(Data_tau,squeeze(Data(ld_pks(i)+3,4,:)))
%     xlabel("\tau (s)");
%     ylabel("g_2(\tau)")
% end
% set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document
% 
% % set(findall(hfig,'-property','Box'),'Box','off') % optional
% % set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(hfig,'Position');
% set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% 
% hfig = figure();
% picturewidth = 6; % set this parameter and keep it forever
% hw_ratio = 2; 
% semilogx(Data_tau,squeeze(mean(Data([ld_pks(3:33)+3],4,:))));
% title("Averaged g_2(\tau)");
% xlabel("\tau(s)");
% ylabel("g_2(\tau)")
% set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document
% 
% % set(findall(hfig,'-property','Box'),'Box','off') % optional
% % set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(hfig,'Position');
% set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%% plotting the signal for the presentation
close all;

strt_p = 1560;
% strt_p = 1*15;
% stop_p = 2*60*20;
stop_p = 1700;

hfig = figure()
fname = "BFi TP"

sig = aDb1;

t = (strt_p:stop_p)/20;
plot(t,(sig(2,strt_p:stop_p)'))
ylim([min(sig(2,strt_p:stop_p))*0.8, max(sig(2,strt_p:stop_p))*1.2 ])
% legend("DCS 1cm","DCS 2.5 cm")
legend("DCS [r_s_d = 2.5 cm]")
xlabel("Time (s)")
ylabel("\alphaDb (cm^2/s)")
title("BFi")

% legend("Normalized 1cm","Norm 2.5cm","Norm 2.5cm TR","DCS 1cm","DCS 2.5cm","TR sys DCS 2.5cm ")
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-painters')


% %%
% close all;
% shift = 450;
% [test_1,pind] = ensemble_avg(ecg1(460000:520000),dcs_1lp(460*20:520*20),shift,1);
% [test_2,pind] = ensemble_avg(ecg1(460000:520000),dcs_25lp(460*20:520*20),shift,1);
% plot(normalize(test_1),'b'); hold on; plot(normalize(test_2),'r');
% 
% % %% Finding the single average pulse to calculate the aDb abd DFC using SF
% % %  % method rather than for the whole signal.
% % % 
% % close all;
% % % if exist("g2")
% % %     clear g2;
% % % end
% % mua = 0.15; %cm^-1 baseline absorption coefficient
% % mus = 8; 
% % 
% %     
% % % 
% % g2_backup = Data;
% % g2_n = g2_backup;
% % strt_time =  [1];   %time in seconds
% % stp_time =  [180];    %time in seconds.
% % 
% % % strt_time =  [1,180,300,420,540];   %time in seconds
% % % stp_time =  [180,300,420,540,660];    %time in seconds.
% % 
% % %  strt_time =  [1,32,62];
% % %  stp_time =  [28,58,90]; 
% % 
% % %For subject 5(Farah), as the protocol is a bit different;
% % % strt_time =  [1,300,480,720];   %time in seconds
% % % stp_time =  [300,480,720,780];    %time in seconds.
% % for m=1:length(strt_time)
% %     disp(m);
% %     clear g2;
% %     stop_time = stp_time(m);
% %     start_time = strt_time(m);
% %     g2 = g2_n(start_time*20:stop_time*20,:,:);
% %     ecg = ecg_ad(start_time*1000:stop_time*1000);
% %     
% %     
% %     [h_pks,l_pks] = findpeaks(normalize(ecg),"MinPeakHeight",2.5,"MinPeakDistance",750);
% %     
% %     fig1=figure('units','centimeters', 'Position',[2 2 35 13]); %18 width 15 heigh
% %     hold on
% %     plot(normalize(ecg),'r')
% %     plot(l_pks, h_pks,'*k')
% %     hd_pks = ceil(h_pks./50);
% %     ld_pks = ceil(l_pks./50);
% %     % g2_avg = zeros(size(g2,1),size(g2,2),size(g2,3));
% %     avg_window_width = 50;
% %     min_length = min(diff(ld_pks(1:avg_window_width)));
% %     for i=1:avg_window_width
% %         
% %         base_sig = g2(ld_pks(i):ld_pks(i)+min_length,:,:);
% %        
% %         for j=1:avg_window_width-1
% %             base_sig = base_sig + g2(ld_pks(i+j):ld_pks(i+j)+min_length,:,:);
% %     %         adb_1 = hybrid_dcs(base_sig,Data_tau);
% %            
% %         end
% %         g2_avg(:,:,:) = base_sig./avg_window_width;
% %     end
% % end
% % 
% % 
% % %% Finding the g2 averaged signal to find DeltaFc
% % % if exist("g2")
% % %     clear 'g2'
% % % end
% % close all;
% % g2_s(:,1,:) = g2_avg(:,1,:)-1;
% % for i=1:size(g2_avg,1)
% %     g2_s(i,2,:)=( g2_avg(i,2,:)+g2_avg(i,3,:)+g2_avg(i,4,:))/3;
% % end
% % g2_a(:,2,:) = g2(:,2,:) -1;
% % 
% % % adb = standalone_dcs(g2_avg,Data_tau);
% % % figure()
% % % plot(adb');
% % 
% % rho = [1 2.5]; %source detector separations in cm 
% % mua = 0.17; %cm^-1 baseline absorption coefficient
% % mus = 8; %cm^-1 baseline reduced scattering coefficient
% % 
% % tau_values=Data_tau;
% % 
% % for chan=1:size(g2_s,2)
% %      for i=1:size(g2_s,1) 
% %         rsd=rho(chan);
% %         g2_temp(i,:)=squeeze(g2_s(i,chan,:));
% %         LB = [0];
% %         UB = [inf];
% %         Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
% %         beta= squeeze(g2_avg(i,chan,1)); %0.1568;
% %         options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
% %         [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rsd,beta);
% %         adb(chan,i) = FittedParams(1);
% %     end
% % end
%%
% Data = g2_n;
if exist("g2")
    clear g2
end
close all;
% g2_n = data;
bsl = 1:550;
pres = 650:1200;

% span = 170:197;
span = 1:size(g2_n,1);
g2(1,:,:)=squeeze(g2_n(span,1,:));
g2_1_temp=squeeze(g2_n(span,1,:)); %g2-1 curve generation
g2_2_temp=squeeze(g2_n(span,2,:)); %g2-1 curve generation
g2_3_temp=squeeze(g2_n(span,3,:)); %g2-1 curve generation
g2_4_temp=squeeze(g2_n(span,4,:)); %g2-1 curve generation

for i=1:size(g2,2)
    g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
end
data = g2;

    % 
% windowsize=2;
% wagi=ones(1,windowsize)/windowsize;
% for r=1:2
%     for i=1:size(data,3)
%         data_temp=squeeze(data(r,:,:));
%         data(r,:,i)=filtfilt(wagi,1,data_temp(:,i));
%     end
% end

g2_temp =round(data,4);
    
% g_b=squeeze(g2(1,100,:))';
% g_b=g_b-mean(g_b(1,end-5:end));
% g_b=(g_b/max(g_b))/2+1;
% 
% g_p=squeeze(g2(1,800,:))';
% g_p=g_p-mean(g_p(1,end-5:end));
% g_p=(g_p/max(g_p))/2+1;

% bsl = 1:550;
% 
% l1 = size(g2,2);
% 
% for i=1:l1
%     g2_new = g2;
%     g2_new(:,i,:) = g2(:,i,:)-mean(g2(:,i,end-5:end),3);
%     g2_new(:,i,:) = (g2_new(:,i,:)./max(g2_new(:,i,:)))/2+1;
% end
% for i=1:l1
%     g_b=squeeze(g2(1,i,:))';
%     g_b=g_b-mean(g_b(1,end-10:end));
%     g_b=(g_b/max(g_b))/2+1;
%     g_b(1) = round(g_b(1),1);
%     g2(1,i,:) = g_b;
% % 
% %     bsl_art=mean(g2_p(1,bsl,end));
% % 
% %     g_b=squeeze(g2(1,i,:))';
% %     g_b=g_b-g_b(1,end);
% %     g_b=(g_b/max(g_b));
% % %     g_b(1) = round(g_b(1),1);
% %     g2(1,i,:) = g_b*0.5+1;
% 
%  
%     g_p=squeeze(g2(2,i,:))';
%     g_p=g_p-mean(g_p(1,end-10:end));
%     g_p=(g_p/max(g_p))/2+1;
%     g_p(1) = round(g_p(1),1);
%     g2(2,i,:) = g_p;
% %     bsl_art=mean(g2_p(2,bsl,end));
% % 
% %     g_b=squeeze(g2(2,i,:))';
% %     g_b=g_b-g_b(1,end);
% %     g_b=(g_b/g_b(1,1));
% % %     g_b(1) = round(g_b(1),1);
% %     g2(2,i,:) = g_b*0.5+1;
% end
% 
% for i=1:10:l1
%     semilogx(Data_tau,squeeze(g2(2,i,:)));
%     hold on;
% end

% g2_temp=round(g2,4);
% g2_temp(g2_temp<=1)=1;

%%
close all;
bsl = 1:1200;
pres = 2450:3500;
    
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

DeltaOD_sp=abs(-log((mean(squeeze(g2_temp(1,pres,:)))-1)./(mean(squeeze(g2_temp(1,bsl,:)))-1)));
DeltaOD_lp=abs(-log((mean(squeeze(g2_temp(2,pres,:)))-1)./(mean(squeeze(g2_temp(2,bsl,:)))-1)));

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


ratio_d = ratio(10:end-20);

index=ratio_d==Inf;

% sens_fact = nanmean(ratio_d(index==0))

sens_fact = mean(ratio_d)

figure()

semilogx(tau,DeltaFec(1,:),'r')
hold on
semilogx(tau,DeltaFec(2,:),'b')
legend("Short Rsd","Long Rsd")

%% Plotting the mean g2 curves for pressure and baseline for comparision
close all;
% semilogx(tau,(mean(squeeze(g2_temp(1,pres,:)))-1),tau,(mean(squeeze(g2_temp(1,bsl,:)))-1))
% semilogx(tau,((squeeze(g2_temp(1,100,:)))),'r');
semilogx(tau,((squeeze(g2_temp(1,100,:)))),tau,((squeeze(g2_temp(1,3900,:)))))
legend("Baseline", "Pressure Modulation")
title("Autocorrelation Functions");
xlabel("\tau (s)");
ylabel("g_2(\tau)")

%%
close all;
hfig = figure();
plot((1:1:200)/20,Data(1:200,1,1),'linewidth',1.5)
xlabel("Time (s)");
ylabel("Intensity")

picturewidth = 15; % set this parameter and keep it forever
hw_ratio = 0.5; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
%% Calculating the weighting factor


%% calculating the flow values

g2_1_new = mean(g2(:,1,:),2);

    rho1 = 10; %mm;
    rho2 = 25;
                    
                  
    corr1 = squeeze(g2_1_new(1,:,:));
    corr2 = squeeze(g2_1_new(2,:,:));
    
    tau1 = Data_tau;
    tau2 = Data_tau;
         
    [beta1,Bspt1,Bept1,spt1,ept1] = PointsSC2019(tau1,corr1,rho1);
    [beta2,Bspt2,Bept2,spt2,ept2] = PointsSC2019(tau2,corr2,rho2);
                    
                  
    beta1=corr1(1,1);
    beta2=corr2(1,1);


        LB.three = [0 0 0];
        UB.three = [inf inf inf];
        Starting = [5e-9 5e-8 0.5]; %[Fec, Fc, top layer thickness; cm^2/s, cm^2/s, cm]
       
        options = optimset('Display','final','TolX',1e-11,'MaxIter', 1000, 'MaxFunEvals', 1000);
        
        FittedParams = fminsearchbnd(@MD_DCS_N2H_BETA,Starting,LB.three,UB.three,options,...
            tau1,tau2,corr1,corr2,ua,us,rho1,rho2,spt1,spt2,ept1,ept2,l1,beta1-1,beta2-1);
        Fscalp = FittedParams(1);
%         Fskull = FittedParams(2);
        Fbrain = FittedParams(2);

        f0c = Fbrain;
        f0ec = Fscalp;
        
F = [f0ec f0c];
%%
tau = Data_tau;

l1 = 0.8; %[cm]
l2 = 0.6;

%Optical Properties
% amp = 1000;
ua1 = 0.15; %[cm^-1]
us1 = 10; %[cm^-1]
ua2 = 0.2; %[cm^-1]
us2 = 8; %[cm^-1]

ua = [ua1 ua2];
us = [us1 us2];

ua1 = 0.17; %[cm^-1]
us1 = 10; %[cm^-1]
ua2 = 0.2; %[cm^-1]
us2 = 8; %[cm^-1]
ua3 = 0.2; %[cm^-1]
us3 = 8; %[cm^-1]
close all;

    f0c = 1.44*10^-8;  % Assumed from the literature, baker pressure paper
    f0ec = 1.4*10^-9; % Assumed from the literature, baker pressure paper
    delFc = f0c*10^-5; 
    fc_denom = f0c + delFc/2;
    fc_num = f0c - delFc/2;
    rho_l = 2.5;
    tau1 = Data_tau;
    wf = (2/delFc)*(log(G1fun(rho_l,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,f0ec,0,fc_num) ...
        ./G1fun(rho_l,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,f0ec,0,fc_denom)));
    
    
    semilogx(Data_tau, (1./wf))
    hold on


%% Multiplying the DeltaFc with 1/wf
% DeltaFc_f =(1./wf).*((DeltaOD_l) - ((0.55).*(DeltaOD_s)));
% 
% DeltaFc_f(isinf(DeltaFc_f)) = nan;
% t = (1:1:size(DeltaFc_f,1))/20;
% % mean_dfc_f = nanmean(DeltaFc_f(:,2:25),2)./f0c;
% mean_dfc_f = nanmean(DeltaFc_f(:,2:25),2)
% mean_dfc_sm_f = smooth(mean_dfc_f,0.015);
% 
% plot(mean_dfc_sm_f,'b'); hold on; plot(mean_dfc_sm,'r')

%% Calculating Delta Fc
close all
if exist("mean_dfc")
    clear mean_dfc
    clear mean_dfc_sm;
    clear  per_c_adb_sm;
    clear per_change_adb;
    clear per_change_dfc;

end

bsl = 100:2500;


% sf = [0.05:0.05:1.5];
% for i=1:1:length(sf)   
    DeltaOD_s = real(-log(((squeeze(g2_temp(1,:,:)))-1)./(mean(squeeze(g2_temp(1,bsl,:)))-1)));
    DeltaOD_l =  real(-log(((squeeze(g2_temp(2,:,:)))-1)./(mean(squeeze(g2_temp(2,bsl,:)))-1)));
    
    DeltaFc = (1./wf).*((DeltaOD_l) - ((0.8).*(DeltaOD_s)));

%     DeltaFc = DeltaFc/f0c;
%     DeltaFc = ((DeltaOD_l) - ((0.78).*(DeltaOD_s)));
%     DeltaFc = (1./wf).*((DeltaOD_l) - ((sf(i)).*(DeltaOD_s)));
    DeltaFc(isinf(DeltaFc)) = nan;
    t = (1:1:size(DeltaFc,1))/20;
    mean_dfc = nanmean((DeltaFc(:,2:25)),2);
    mean_dfc_f = lpf(mean_dfc,7,20);
    mean_dfc_sm = smooth(mean_dfc,0.0085,'rloess');
    minimum_val(i) = min(mean_dfc_sm);
%     per_change_dfc = mean_dfc_sm*100;
    per_change_dfc = (mean_dfc_sm-mean(mean_dfc_sm(300:1800)))/mean(mean_dfc_sm(300:1800));
%     per_change_dfc = (mean_dfc_sm)/mean(mean_dfc_sm(300:3200));
%     per_change_dfc = (mean_dfc_sm/f0c);
%     per_change_adb(1,:) = ((aDb1(1,:)-mean(aDb1(1,1:600),2))./mean(aDb1(1,1:600),2))*100;
%     per_change_adb(2,:) = ((aDb1(2,:)-mean(aDb1(2,1:600),2))./mean(aDb1(2,1:600),2))*100;
%     per_c_adb_sm(1,:) = smooth(per_change_adb(1,:),0.05);
%     per_c_adb_sm(2,:) = smooth(per_change_adb(2,:),0.05);
% end
plot(t,real(mean_dfc)); xlabel("Time (s)"); ylabel("\DeltaF_c"); title("\DeltaF_c pulsatility")
% hold on; 
figure(); 
plot(mean_dfc_sm)
figure();
plot(t,per_change_dfc);
% ylim([-1 1])
xlabel("Time (s)")
ylabel("Mean DeltaFc")
title("Mean DeltaFc across tau(2:25)")
%%
hfig = figure(); 
t_adb = (1:1:length(per_ch_mean))
plot(t_adb,per_ch_mean(1,:)-100);
hold on;
plot(t_adb,per_ch_mean(2,:)-100);
plot(t,per_change_dfc)

xlabel("Time (s)")
hold on;
ylabel("% change")
title("Pecrent Change in aDb wrt Baseline")
legend("aDb 1cm","aDb 2.5cm","\DeltaF_c")

picturewidth = 25; % set this parameter and keep it forever
hw_ratio = 0.5; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
% print(hfig,fname,'-dpng','-painters')

% figure();
% plot(t, change_dfc_flip);

% csvwrite("2 layer model\009_DFc.csv",mean_dfc_sm)
%% Plotting the ensemble average graph
[ens_avg_sig,pind] = ensemble_avg(ecg1(1:90000),mean_dfc(1:1800),500,1);

%%
% Calculating FT frequency spectrum
close all;
Fs = 20;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
signal = dcs_1lp(4801:6000);
L = length(signal);             % Length of signal
t = (0:L-1)*T;  

Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot((1:1:length(P1))/Fs,P1)

%% Pulsatility in DeltaFc
close all;
%For subjects 1,2,3,6
strt_time =  [20,190,350,480,540];   %time in seconds
stp_time =  [60,280,400,510,660];    %time in seconds.
%For subject 5(Farah), as the protocol is a bit different;
% strt_time =  [10,320,590,720];   %time in seconds
% stp_time =  [280,460,700,760];    %time in seconds.
shift = 1200;


j =1;
[ens_avg_dcs1_1,pind,ensemble_curve_dcs1_1] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_1lp(strt_time(j)*20:stp_time(j)*20),shift,1);
[ens_avg_dcs25_1,pind,ensemble_curve_dcs25_1] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_25lp(strt_time(j)*20:stp_time(j)*20),shift,1);
[ens_avg_dfc_1,pind,ensemble_curve_dfc_1] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),mean_dfc(strt_time(j)*20:stp_time(j)*20),shift,1);

% j=4;
% [ens_avg_dcs1_2,pind,ensemble_curve_dcs1_2] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_1lp(strt_time(j)*20:stp_time(j)*20),shift,1);
% [ens_avg_dcs25_2,pind,ensemble_curve_dcs25_2] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),dcs_25lp(strt_time(j)*20:stp_time(j)*20),shift,1);
% [ens_avg_dfc_2,pind,ensemble_curve_dfc_2] = ensemble_avg(ecg1(strt_time(j)*1000:stp_time(j)*1000),mean_dfc(strt_time(j)*20:stp_time(j)*20),shift,1);

%% Normalizing the waveforms
close all;
t = (1:1:length(ens_avg_dcs25_1))/1000;
bfi = ens_avg_dcs25_1;
dfc = ens_avg_dfc_1;
bfi_n = (bfi-min(bfi(1:300)))./(max(bfi(1:500))-min(bfi(1:300)));
dfc_n = (dfc-min(dfc(1:300)))./(max(dfc(1:500))-min(dfc(1:300)));
sig = [bfi_n;dfc_n];
% plot(bfi_n,'k');
% hold on; 
% plot(dfc_n,'r');
% 
hfig = figure();
fname = "Pulsatility Comparison"

% plot(t,circshift(bfi_n(1:1000),-30),'k');
t = (1:1:length(bfi_n))/1000;
plot(t,(bfi_n),'k');
hold on;
% plot(t,circshift(dfc_n(1:1000),-30),'r');
plot(t,dfc_n,'r');
xlabel("Time(s)")
ylabel("A.U.")

legend("BFi","\DeltaF_c")

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.75; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-painters')
%%
close all;


hfig = figure()
fname = "BFi TP"
sgtitle("Comparison of Pulsatility")
t = (1:1:length(ens_avg_dcs1_1))/1000;

% subplot(1,3,1)
% plot(t,ens_avg_dcs1_1)
% % ylim([-2 2])
% title("DCS 1cm Baseline")

t = (1:1:length(ens_avg_dcs25_1))/1000;
subplot(1,2,1)
plot(t,ens_avg_dcs25_1 -min(ens_avg_dcs25_1))
% plot(t,bfi_n)
xlabel("Time (s)")
ylabel("\alphaDb (Normalized)")
% ylim([-2 2])
title("BFi_2_._5_c_m")

t = (1:1:length(ens_avg_dfc_1))/1000;
subplot(1,2,2)
plot(t,ens_avg_dfc_1 - min(ens_avg_dfc_1 ))
% plot(t,dfc_n)
xlabel("Time (s)")
ylabel("\DeltaF_c (Normalized)")
% ylim([-2 2])
title("\DeltaF_c Baseline")

figure()
t = (1:1:length(ens_avg_dcs25_1))/1000;
subplot(1,2,1)
% plot(t,ens_avg_dcs25_1 -nanmean(ens_avg_dcs25_1))
plot(t,normalize(ens_avg_dcs25_1)-min(normalize(ens_avg_dcs25_1)))
% plot(t,ens_avg_dcs25_1)
xlabel("Time (s)")
ylabel("\alphaDb (cm^2/s)")
% ylim([-2 2])
title("DCS 2.5cm Baseline")

t = (1:1:length(ens_avg_dfc_1))/1000;
subplot(1,2,2)
% plot(t,ens_avg_dfc_1 - nanmean(ens_avg_dfc_1 ))
plot(t,normalize(ens_avg_dfc_1)- min(normalize(ens_avg_dfc_1 )))
% plot(t,ens_avg_dfc_1)
xlabel("Time (s)")
ylabel("\DeltaF_c")
% ylim([-2 2])
title("\DeltaF_c Baseline")

% t = (1:1:length(ens_avg_dcs1_2))/1000;
% subplot(2,3,4)
% plot(t,ens_avg_dcs1_2)
% % ylim([-2 2])
% title("DCS 1cm 150% Tourniquet inflation")
% 
% subplot(2,3,5)
% t = (1:1:length(ens_avg_dcs25_2))/1000;
% plot(t,ens_avg_dcs25_2)
% % ylim([-2 2])
% title("DCS 2.5cm 150% Tourniquet inflation")
% 
% t = (1:1:length(ens_avg_dfc_2))/1000;
% subplot(2,3,6)
% plot(t,ens_avg_dfc_2)
% % ylim([-2 2])
% title("\DeltaF_c 150% Tourniquet inflation")




% legend("DCS 1cm","DCS 2.5 cm")
% xlabel("Time (s)")
% ylabel("Normalized Unit")
% title("BFi")

% legend("Normalized 1cm","Norm 2.5cm","Norm 2.5cm TR","DCS 1cm","DCS 2.5cm","TR sys DCS 2.5cm ")
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.75; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-painters')

%%
a = (nanmean(ens_avg_dfc_1)-ens_avg_dfc_1)./nanmean(ens_avg_dfc_1);
plot(a)
%%
close all;


hfig = figure()
fname = "BFi Comparision"
sgtitle("Comparison of Pulastility")

t = (1:1:length(ens_avg_dcs25_1))/1000;
subplot(1,2,1)
plot(t,ens_avg_dcs25_1)
% ylim([-2 2])
title("DCS 2.5cm Baseline")

t = (1:1:length(ens_avg_dfc_1))/1000;
subplot(1,2,2)
plot(t,ens_avg_dfc_1)
% ylim([-2 2])
title("\DeltaF_c Baseline")

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',14) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-painters')


%% Subtract the 1cm signal from 2.5cm signal
close all;
s_sub = -1*ens_avg_dcs25_1+ens_avg_dcs1_1(1:length(ens_avg_dcs25_1));

figure();
subplot(1,4,1);
t = (1:1:length(ens_avg_dcs1_1))/1000;
plot(t,ens_avg_dcs1_1)
title("DCS 1cm bsl");

subplot(1,4,2);
t = (1:1:length(ens_avg_dcs25_1))/1000;
plot(t,ens_avg_dcs25_1)
title("DCS 2.5cm bsl");

subplot(1,4,3)
t = (1:1:length(s_sub))/1000;
plot(t,s_sub)
title("DCS(2.5)-DCS(1)");

subplot(1,4,4)
t = (1:1:length(ens_avg_dfc_1))/1000;
plot(t,ens_avg_dfc_1)
title("Mean DFc");

%% 
close all;
for i=1:5:100
    semilogx(Data_tau,squeeze(g2_sub(i,1,:)))
    hold on;
end




%% Subtracting the g2 curves of DCS 1cm from 2.5cm and fitting for adb

g2_t1(:,1,:) = g2_n(:,1,:);
for i=1:size(g2_n,1)
    g2_t1(i,2,:)=( g2_n(i,2,:)+g2_n(i,3,:)+g2_n(i,4,:))/3;
end

g2_sub = g2_t1(:,1,:)-g2_t1(:,2,:);
semilogx(Data_tau,squeeze(g2_t1(5,1,:)),'b')
hold on;
semilogx(Data_tau,squeeze(g2_t1(5,2,:)),'r')
semilogx(Data_tau,squeeze(g2_sub(5,1,:)),'k')
%% fit for adb values
rho = [2.5]; %source detector separations in cm 
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 8; %cm^-1 baseline reduced scattering coefficient

tau_values=Data_tau;

for chan=1:size(g2_sub,2)
     for i=1:size(g2_sub,1) 
        rsd=rho(chan);
        g2_sub1(i,:)=squeeze(g2_sub(i,chan,:));
        LB = [0];
        UB = [inf];
        Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
        beta= squeeze(g2_sub(i,chan,1)); %0.1568;
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_sub1(i,:),mua,mus,rsd,beta);
        adb_sub(chan,i) = FittedParams(1);
    end
end

%%
plot(adb_sub);
%% Percent change in aDb
close all;
adb_sm  = normalize(aDb1(1,:));
adb_lg =  normalize(aDb1(2,:));

per_adb_sm = ((adb_sm- mean(adb_sm(1:1800)))/mean(adb_sm(1:1800)))*100;
per_adb_lg = ((adb_lg- mean(adb_lg(1:1800)))/mean(adb_lg(1:1800)))*100;

plot(per_adb_sm,'b'); hold on; plot(per_adb_lg,'r');

%% Assessing Optical Densities for long and short rsd
close all;
mean_od_l = nanmean(DeltaOD_l(:,2:15),2);
mean_od_s = nanmean(DeltaOD_s(:,2:15),2);
mean_od_l_sm = smooth(mean_od_l,0.05);
mean_od_s_sm = smooth(mean_od_s,0.05);
plot(t,mean_od_s_sm,'r'); hold on; plot(t,mean_od_l_sm,'b');plot(t,(mean_od_l_sm-0.45*mean_od_s_sm),'k')
xlabel("Time(s)")
ylabel("OD")
title("OD across time")
legend("r_s_d = 1 cm","r_s_d = 2.5 cm","OD_{long}- (ratio)*OD_{short}")
%%
sf = [0.1:0.4:5];
for i = 1:length(sf)
    DeltaFc = ((DeltaOD_l) - ((sf(i)).*(DeltaOD_s)));
    DeltaFc(isinf(DeltaFc)) = nan;
    t = (1:1:size(DeltaFc,1))/20;
    mean_dfc = nanmean(DeltaFc(:,2:25),2);
    mean_dfc_sm(i,:) = squeeze(smooth(mean_dfc,0.05));
    per_change_dfc(i,:) = mean_dfc_sm(i,:)*100;

    % plot(mean_dfc); hold on; 
%     plot(t,mean_dfc_sm);
%     xlabel("Time (s)")
%         ylabel("Mean DeltaFc")
%     title("Mean DeltaFc across tau(2:25)")
%     hold on;
figure()
    plot(t, per_change_dfc(i,:));
    xlabel("Time (s)")
    ylabel("% change")
    title("Pecrent Change in aDb wrt Baseline")
    hold on;
    legend(string(sf(i)))
end

%%



%%
tic
for i = 1:1
    standalone_dcs(Data,Data_tau);
end
toc
%%
g_b=squeeze(g2(1,100,:))';
g_b=g_b-mean(g_b(1,end-5:end));
g_b=(g_b/max(g_b))/2+1;

g_p=squeeze(g2(1,800,:))';
g_p=g_p-mean(g_p(1,end-5:end));
g_p=(g_p/max(g_p))/2+1;

% hold on
semilogx(Data_tau, g_b)
hold on
semilogx(Data_tau, g_p)

%% Reading all the files at once to calculate the mean sensitivity
f_ind = [1,3,5,9,11,15,17];
ratio = zeros(length(f_ind),50);
DeltaFc = zeros(length(f_ind),1800,50);
for j = 1:length(f_ind)

    file_name = strcat("D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\DCS\20220608 - ",string(f_ind(j)),"\Data.mat");
%     load("D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\DCS\20220608 - 9\Data.mat");
    load(file_name);
%     aDb1 = standalone_dcs(Data,Data_tau);
    tau = Data_tau;
    
    
    close all;
    g2(1,:,:)=squeeze(Data(:,1,:)); %g2-1 curve generation
    g2_1_temp=squeeze(Data(:,1,:)); %g2-1 curve generation
    g2_2_temp=squeeze(Data(:,2,:)); %g2-1 curve generation
    g2_3_temp=squeeze(Data(:,3,:)); %g2-1 curve generation
    g2_4_temp=squeeze(Data(:,4,:)); %g2-1 curve generation
    
    for i=1:size(g2,2)
        g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
    end
    
    % g_b=squeeze(g2(1,100,:))';
    % g_b=g_b-mean(g_b(1,end-5:end));
    % g_b=(g_b/max(g_b))/2+1;
    % 
    % g_p=squeeze(g2(1,800,:))';
    % g_p=g_p-mean(g_p(1,end-5:end));
    % g_p=(g_p/max(g_p))/2+1;
    
    l1 = size(g2,2);
    
    % for i=1:l1
    %     g2_new = g2;
    %     g2_new(:,i,:) = g2(:,i,:)-mean(g2(:,i,end-5:end),2);
    %     g2_new(:,i,:) = (g2_new(:,i,:)./max(g2_new(:,i,:)))/2+1;
    % end
    for i=1:l1
        g_b=squeeze(g2(1,i,:))';
        g_b=g_b-mean(g_b(1,end-10:end));
        g_b=(g_b/max(g_b))/2+1;
        g2(1,i,:) = g_b;
    
    
        g_p=squeeze(g2(2,i,:))';
        g_p=g_p-mean(g_p(1,end-10:end));
        g_p=(g_p/max(g_p))/2+1;
        g2(2,i,:) = g_p;
    end
    
    for i=1:100:l1
        semilogx(Data_tau,squeeze(g2(1,i,:)));
        hold on;
    end
    
    close all;
    bsl = 400:600;
    pres = 800:1000;
    
    g2_temp=round(g2,4);
    g2_temp(g2_temp<=1)=1;
    
    figure();
    semilogx(Data_tau,(mean(squeeze(g2_temp(1,pres,:)))));
    hold on;
    semilogx(Data_tau,(mean(squeeze(g2_temp(1,bsl,:)))));
    semilogx(Data_tau,(mean(squeeze(g2_temp(2,pres,:)))));
    semilogx(Data_tau,(mean(squeeze(g2_temp(2,bsl,:)))));
    legend("Pressure1","Baseline1","Pressure2","Baseline2")
    
    hold off;
    
    figure();
    
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
    
    ratio(j,:) = DeltaFec(2,:)./DeltaFec(1,:);
    rat_smooth(j,:) = smooth(ratio(j,:),0.4,'rloess');
    semilogx(tau,DeltaFec(2,:)./DeltaFec(1,:));

    DeltaFc(j,:,:) = ((DeltaOD_l) - (rat_smooth(j,:).*(DeltaOD_s)));
    
    
    ratio_d = ratio(3:end-5);
    
    index=ratio_d==Inf;
    
    % sens_fact = nanmean(ratio_d(index==0))
    
    sens_fact(j) = mean(ratio_d);
    
    figure()
    
    semilogx(tau,DeltaFec(1,:),'r')
    hold on
    semilogx(tau,DeltaFec(2,:),'b')
    legend("Short Rsd","Long Rsd")



end

%%
close all
for i=1:size(ratio,1)
%     figure();
    semilogx(Data_tau,ratio(i,:)','DisplayName','rat_smooth'); hold on;  
    legend("Smooth","Raw")
end
figure();
semilogx(Data_tau,mean(rat_smooth));
avg_ratio = mean(rat_smooth);
avg_sf = mean(avg_ratio(2:35));

%%
avg_sf_all = mean(ratio(:,2:35),2)
sd = std(avg_sf_all)
%%
close all
for i=1:100:size(DeltaFc,2)
    semilogx(Data_tau,squeeze(DeltaFc(1,i,:)));
    hold on;
end

%%
close all
G1 =  G1_fun_2_layer(2.5,Data_tau,ua,us,0.5,F);
g1 = G1./max(G1);
g2 = beta1*(g1.^2);
semilogx(Data_tau,g2);
hold on;
G1 =  G1_fun_2_layer(2.5,Data_tau,ua,us,1.5,F);
g1 = G1./max(G1);
g2 = beta1*(g1.^2);
semilogx(Data_tau,g2);

%%  Parkwood study analysis
f_ind = [1:3,5,7,8,10:12,14:18,21,22];
   
ratio = zeros(length(f_ind),50);
for j = 1:length(f_ind)

     if f_ind(j)<10
        file_name = strcat("MCI00",string(f_ind(j)),"_PRESSURE_Data.mat");
    else
        file_name = strcat("MCI0",string(f_ind(j)),"_PRESSURE_Data.mat");
    end
    dir = "D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\Parkwood Study\MATfiles\";
    file_path  = strcat(dir,file_name);
    load(file_path);

%     file_name = strcat("D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\DCS\20220608 - ",string(f_ind(j)),"\Data.mat");
%     load(file_name);

    tau = Data_tau;
    
    
    close all;
    g2(1,:,:)=squeeze(Data(:,1,:)); %g2-1 curve generation
    g2_1_temp=squeeze(Data(:,1,:)); %g2-1 curve generation
    g2_2_temp=squeeze(Data(:,2,:)); %g2-1 curve generation
    g2_3_temp=squeeze(Data(:,3,:)); %g2-1 curve generation
    g2_4_temp=squeeze(Data(:,4,:)); %g2-1 curve generation
    
%     for i=1:size(g2,2)
%         g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
%     end
    g2(2,:,:) = g2_4_temp;
    
    l1 = size(g2,2);
    
    for i=1:l1
        g_b=squeeze(g2(1,i,:))';
        g_b=g_b-mean(g_b(1,end-10:end));
        g_b=(g_b/max(g_b))/2+1;
        g2(1,i,:) = g_b;
    
    
        g_p=squeeze(g2(2,i,:))';
        g_p=g_p-mean(g_p(1,end-10:end));
        g_p=(g_p/max(g_p))/2+1;
        g2(2,i,:) = g_p;
    end
    
    for i=1:100:l1
        semilogx(Data_tau,squeeze(g2(1,i,:)));
        hold on;
    end
    
    close all;
    bsl = 100:550;
    pres = 700:1100;
    
    g2_temp=round(g2,4);
    g2_temp(g2_temp<=1)=1;
    
    figure();
    semilogx(Data_tau,(mean(squeeze(g2_temp(1,pres,:)))));
    hold on;
    semilogx(Data_tau,(mean(squeeze(g2_temp(1,bsl,:)))));
    semilogx(Data_tau,(mean(squeeze(g2_temp(2,pres,:)))));
    semilogx(Data_tau,(mean(squeeze(g2_temp(2,bsl,:)))));
    legend("Pressure1","Baseline1","Pressure2","Baseline2")
    
    hold off;
    
    figure();
    
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
    
    ratio(j,:) = DeltaFec(2,:)./DeltaFec(1,:);
    rat_smooth(j,:) = smooth(ratio(j,:),0.4,'rloess');
    semilogx(tau,DeltaFec(2,:)./DeltaFec(1,:));

%     DeltaFc(j,:,:) = ((DeltaOD_l) - (rat_smooth(j,:).*(DeltaOD_s)));
    
    
    ratio_d = ratio(j,3:end-15);
    
    index=ratio_d==Inf;
    
    % sens_fact = nanmean(ratio_d(index==0))
    
    sens_fact(j) = mean(ratio_d);
    
    figure()
    
    semilogx(tau,DeltaFec(1,:),'r')
    hold on
    semilogx(tau,DeltaFec(2,:),'b')
    legend("Short Rsd","Long Rsd")



end

%% Assessing the changes in the g2 for pressure modulation and SF study
g2(1,:,:) = Data(:,1,:);
for i=1:size(Data,1)
    g2(2,i,:) = (Data(i,2,:)+Data(i,3,:)+Data(i,4,:))/3;
end
 
%% Plotting the adb signal for comparision

plot(adb_avg(:,[1:3631,6014:8576,10897:end])')