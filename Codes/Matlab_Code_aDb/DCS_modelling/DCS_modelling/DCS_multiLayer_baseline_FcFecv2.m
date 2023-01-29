fclose('all'); % Close all open files
clear all; 
close all;
clc;

%%

% tau=load('data_tau.csv')
% load('BFi.mat')
% 
% % data_cut= [510 570; 30 90; 30 90; 30 90; 270 330; 270 330; 30 90; 510 570; 30 90];
% % data_cut2=[270 330; 90 150; 90 150; 90 150; 90 150; 90 150; 90 150; 90 150; 90 150];
% 
% for subject=1:9
%     filename=strcat(pwd,'\2 layer model\00',num2str(subject),'_TCD_g2.mat')    
%     load(filename);
% 
%     baseline=60; % in seconds
%     baseline_nop=baseline/0.05; % in points
% 
%     g2_baseline=squeeze(mean(g2(:,1:baseline_nop,:),2));
% 
%     filename=strcat(pwd,'\2 layer model\00',num2str(subject),'_CUFF_g2.mat')    
%     load(filename);
% 
%     g2_cuff=squeeze(mean(g2(:,baseline_nop+1:2*baseline_nop,:),2));
% 
%     mua = OpticalProp(subject,1); %cm^-1 baseline absorption coefficient
%     mus = OpticalProp(subject,2); %cm^-1 baseline reduced scattering coefficient
% 
%     LB = [0.1e-9];
%     UB = [10e-8];
%     Starting = 1e-9; %[aDb, Beta; cm^2/s, a.u.]
%     beta= g2_baseline(1,1); %0.1568;
%     options = optimset('Display','final','TolX',1e-10,'MaxIter',2000, 'MaxFunEvals', 2000);
%     [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau,g2_baseline(1,:),mua,mus,1,beta);
%     aDb1 = FittedParams(1);
%     
%     beta= g2_cuff(1,1); %0.1568;
%     options = optimset('Display','final','TolX',1e-10,'MaxIter',2000, 'MaxFunEvals', 2000);
%     [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau,g2_cuff(1,:),mua,mus,1,beta);
%     aDb2 = FittedParams(1);
%     
%     DeltaFec(subject)=(aDb2/aDb1);
%     dFec=DeltaFec;
%     
%     if subject==1
%         l1 = 0.6; %[cm]
%         l2 = 0.4;
%     elseif subject==2
%         l1 = 0.6; %[cm]
%         l2 = 0.8;
%     elseif subject==3
%         l1 = 0.6; %[cm]
%         l2 = 0.8;    
%     elseif subject==4
%         l1 = 0.6; %[cm]
%         l2 = 0.8;    
%     elseif subject==5
%         l1 = 0.6; %[cm]
%         l2 = 0.3;   
%     elseif subject==6
%         l1 = 0.8; %[cm]
%         l2 = 0.4;     
%     elseif subject==7
%         l1 = 0.6; %[cm]
%         l2 = 0.3;     
%     elseif subject==8
%         l1 = 0.4; %[cm]
%         l2 = 0.5;
%     elseif subject==9
%         l1 = 0.7; %[cm]
%         l2 = 0.5;    
%         
%     end
%     % Optical Properties
% 
%     ua1 = OpticalProp(subject,1); %[cm^-1]
%     us1 = OpticalProp(subject,2); %[cm^-1]
%     ua2 = 0.6*OpticalProp(subject,1); %[cm^-1]
%     us2 = OpticalProp(subject,2); %[cm^-1]
%     ua3 = OpticalProp(subject,1); %[cm^-1]
%     us3 = OpticalProp(subject,2); %[cm^-1]
%     % ua1 = 0.17; %[cm^-1]
%     % us1 = 8; %[cm^-1]
%     % ua2 = 0.2; %[cm^-1]
%     % us2 = 0.8; %[cm^-1]
%     % % ua3 = 0.2; %[cm^-1]
%     % us3 = 8; %[cm^-1]
% 
% 
%     rho1 = 25; %mm
%     rho2 = 25; %mm;
%                        
%     corr1 = 1+squeeze(g2_baseline(4,:));
%     corr2 = 1+squeeze(g2_cuff(4,:));
%     
%     tau1 = tau;
%     tau2 = tau;
% 
%     beta1=corr1(1,1)-1;
%     beta2=corr2(1,1)-1;
% 
% 
%     LB = [1e-10 0 1e-10 0.5 0.5];
%     UB = [ 1e-7 0 1e-7  0.8 0.8];
%     Starting = [1e-9 0 1e-9 0.6 0.6]; %[Fec, Fc, top layer thickness; cm^2/s, cm^2/s, cm]
% 
%     options = optimset('Display','final','TolX',1e-5,'MaxIter', 2000, 'MaxFunEvals', 2000);
%         
%     FittedParams = fminsearchbnd(@MD_DCS_N2H_BETA_Fc_v2,Starting,LB,UB,options,...
%             tau1,tau2,corr1,corr2,ua1,us1,ua2,us2,ua3,us3,rho1,rho2,beta1,beta2,dFec);
%           
%     Fscalp(subject) = FittedParams(1);
%     Fskull(subject) = FittedParams(2);
%     Fbrain(subject) = FittedParams(3);
% %     DeltaFec(subject) = FittedParams(3);
% close all
% 
% Fbrain./Fscalp
% end

%%

adb = standalone_dcs(Data,Data_tau);
figure()
plot(adb');
%%
close all;


g2(1,:,:)=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2_1_temp=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(Data(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(Data(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(Data(:,4,:)-1); %g2-1 curve generation


% 
g2_baseline=squeeze(mean(g2(:,1:2800,:),2));

g2_cuff=squeeze(mean(g2(:,3800:4800,:),2));

% mua = OpticalProp(subject,1); %cm^-1 baseline absorption coefficient
% mus = OpticalProp(subject,2); %cm^-1 baseline reduced scattering coefficient
mua = 0.17;
mus = 8;
tau = Data_tau;
LB = [0.1e-9];
UB = [10e-8];
Starting = 1e-9; %[aDb, Beta; cm^2/s, a.u.]
beta= g2_baseline(1,1); %0.1568;
options = optimset('Display','final','TolX',1e-10,'MaxIter',2000, 'MaxFunEvals', 2000);
[FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau,g2_baseline(1,:),mua,mus,1,beta);
aDb1 = FittedParams(1);

beta= g2_cuff(1,1); %0.1568;
options = optimset('Display','final','TolX',1e-10,'MaxIter',2000, 'MaxFunEvals', 2000);
[FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau,g2_cuff(1,:),mua,mus,1,beta);
aDb2 = FittedParams(1);

DeltaFec=(aDb2/aDb1);
dFec=DeltaFec;




roi = 3700:4800;

g2_short = g2(1,roi,:);
for i=1:size(g2,2)
    g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
end
g2_long = g2(2,roi,:);
l1 = 0.7; %[cm]
l2 = 0.4;

% Optical Properties

%     ua1 = OpticalProp(subject,1); %[cm^-1]
%     us1 = OpticalProp(subject,2); %[cm^-1]
%     ua2 = 0.6*OpticalProp(subject,1); %[cm^-1]
%     us2 = OpticalProp(subject,2); %[cm^-1]
%     ua3 = OpticalProp(subject,1); %[cm^-1]
%     us3 = OpticalProp(subject,2); %[cm^-1]
ua1 = 0.17; %[cm^-1]
us1 = 8; %[cm^-1]
ua2 = 0.2; %[cm^-1]
us2 = 0.8; %[cm^-1]
ua3 = 0.2; %[cm^-1]
us3 = 8; %[cm^-1]


rho1 = 10; %mm distance for short channel               
rho2 = 25; %mm distance for long channle
                   
corr1 = 1+squeeze(g2_short);
corr2 = 1+squeeze(g2_long);


tau1 = Data_tau;
tau2 = Data_tau;


beta1=corr1(1,1)-1;
beta2=corr2(1,1)-1;


LB = [1e-10 0 1e-15 0.2 0.2];
UB = [ 1e-7 0 1e-5  0.8 0.8];
Starting = [1e-9 0 1e-9 0.6 0.6]; %[Fec, Fc, top layer thickness; cm^2/s, cm^2/s, cm]

options = optimset('Display','final','TolX',1e-5,'MaxIter', 2000, 'MaxFunEvals', 2000);
    
FittedParams = fminsearchbnd(@MD_DCS_N2H_BETA_Fc_v2,Starting,LB,UB,options,...
        tau1,tau2,corr1,corr2,ua1,us1,ua2,us2,ua3,us3,rho1,rho2,beta1,beta2,dFec);
      
Fscalp = FittedParams(1);
Fskull = FittedParams(2);
Fbrain = FittedParams(3);
%     DeltaFec(subject) = FittedParams(3);
close all

Fbrain./Fscalp

%%
% load('Delta_Fc.mat')
% 
% DeltaFc_temp=mean(DeltaFc(:,2:30),2);

for subject=1:9
    [G1_first] = G1fun(rho1/10,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,Fscalp(subject),.01*Fscalp(subject),(Fbrain(subject)-(0.1*Fbrain(subject)/2)));
    [G1_second] = G1fun(rho1/10,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,Fscalp(subject),.01*Fscalp(subject),(Fbrain(subject)+(0.1*Fbrain(subject)/2)));
    
    G_ratio=(G1_first(1,1:30)./G1_second(1,1:30));
    dFc(subject)=sum((2./(0.1*Fbrain(subject)))*(log(G_ratio)));
end

l1=0.8;
l2=0.6;
for subject=1:9
    [G1_first] = G1fun(rho1/10,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,0.16*Fbrain(subject),.01*Fscalp(subject),(Fbrain(subject)-(0.1*Fbrain(subject)/2)));
    [G1_second] = G1fun(rho1/10,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,0.16*Fbrain(subject),.01*Fscalp(subject),(Fbrain(subject)+(0.1*Fbrain(subject)/2)));
    
    G_ratio=(G1_first(1,1:30)./G1_second(1,1:30));
    dFc2(subject)=sum((2./(0.1*Fbrain(subject)))*(log(G_ratio)));
end

%%

l1=0.8;
l2=0.2;

for subject=1
    
    Fbrain_mod=1e-8;
    Fscalp_mod=2e-9;
    
    delta_value=0.00001;
    
    [G1_first]  = G1fun(rho1/10,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,Fscalp_mod,0.01*Fscalp_mod,(Fbrain_mod-(delta_value*Fbrain_mod/2)));
    [G1_second] = G1fun(rho1/10,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,Fscalp_mod,0.01*Fscalp_mod,(Fbrain_mod+(delta_value*Fbrain_mod/2)));
   
    semilogx(tau,G1_first)
    hold on
    semilogx(tau,G1_second)
    
    
    G_ratio=(G1_first(1,1:30)./G1_second(1,1:30));
    dFc3(subject)=sum((2./(delta_value*Fbrain_mod))*(log(G_ratio)));
end

% % Field autocorrelation function
% g1_first = (G1_first)./max(G1_first);
% g2_1_first = (beta1)*(g1_first.^2);
% g2_first = g2_1_first + 1;
% 
% semilogx(tau, g2_baseline(4,:)+1,'k')
% hold on
% semilogx(tau, g2_first,'r')


% n = 60; % average every n values
% for chan=1:4
%     for tau=1:50
%         a=squeeze(g2_1_temp(chan,:,tau));
%         g2_1(chan,tau,:) = arrayfun(@(i) mean(a(i:i+n-1)),1:n:length(a)-n+1)'; % the averaged vecto
%     end
% end
% clear g2
% 
% %
% 
% rsd=[1 1.5 2 2.5];
% 
% % mua = 0.1; %cm^-1 baseline absorption coefficient
% % mus = 10; %cm^-1 baseline reduced scattering coefficient
% 
% mua = OpticalProp(subject,1); %cm^-1 baseline absorption coefficient
% mus = OpticalProp(subject,2); %cm^-1 baseline reduced scattering coefficient
% 
% temp_res=10;
% 
% baseline_t=10; % time in seconds
% baseline=baseline_t/temp_res;
% tau_values=data_tau;
% 
% for r=1:4
%         for i=1:size(g2_1,3)
%             g2_temp(i,:)=squeeze(g2_1(r,:,i));
%             LB = [0.1e-9];
%             UB = [10e-8];
%             rho=rsd(r);
%             Starting = 1e-9; %[aDb, Beta; cm^2/s, a.u.]
%             beta= g2_1(r,1,i); %0.1568;
%             % beta= squeeze(g2_new(1,1)); %0.1568;
%             options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
%             [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rho,beta);
%             aDb1(r,i) = FittedParams(1);
%         end
% end
%             
% 
% filename=strcat(pwd,'\2 layer model\00',num2str(subject),'_CUFF_g2.mat')    
% load(filename);
% 
% g2_2_temp=g2(:,data_cut2(subject,1)/0.05+1:data_cut2(subject,2)/0.05,:);
% 
% n = 60; % average every n values
% for chan=1:4
%     for tau=1:50
%         a=squeeze(g2_2_temp(chan,:,tau));
%         g2_1(chan,tau,:) = arrayfun(@(i) mean(a(i:i+n-1)),1:n:length(a)-n+1)'; % the averaged vecto
%     end
% end
% clear g2
% 
% rsd=[1 1.5 2 2.5];
% 
% % mua = 0.1; %cm^-1 baseline absorption coefficient
% % mus = 10; %cm^-1 baseline reduced scattering coefficient
% 
% mua = OpticalProp(subject,1); %cm^-1 baseline absorption coefficient
% mus = OpticalProp(subject,2); %cm^-1 baseline reduced scattering coefficient
% 
% temp_res=10;
% 
% baseline_t=10; % time in seconds
% baseline=baseline_t/temp_res;
% tau_values=data_tau;
% 
% for r=1:4
%         for i=1:size(g2_1,3)
%             g2_temp(i,:)=squeeze(g2_1(r,:,i));
%             LB = [0.1e-9];
%             UB = [10e-8];
%             rho=rsd(r);
%             Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
%             beta= g2_1(r,1,i); %0.1568;
%             % beta= squeeze(g2_new(1,1)); %0.1568;
%             options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
%             [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rho,beta);
%             aDb2(r,i) = FittedParams(1);
%         end
% end
%     
% BFi_all(subject,:,:)=aDb1;
% BFi_all2(subject,:,:)=aDb2;

% end

%%

% Y1=squeeze(mean(g2_1_temp(1,1:100,:)))
% Y1=Y1/max(Y1)*0.5
% semilogx(tau_values,Y1,'k')
% hold on
% Y2=squeeze(mean(g2_1_temp(2,1:100,:)))
% Y2=Y2/max(Y2)*0.5
% semilogx(tau_values,Y2,'g')
% 
% Y3=squeeze(mean(g2_1_temp(3,1:100,:)))
% Y3=Y3/max(Y3)*0.5
% semilogx(tau_values,Y3,'b')
% 
% Y4=squeeze(mean(g2_1_temp(4,1:100,:)))
% Y4=Y4/max(Y4)*0.5
% semilogx(tau_values,Y4,'magenta')
%% Analyze INDIVIDUAL
% Determine the Blood Flow Index for each Channel
% tau=data_tau;
% time=(1:1:20)*3;
% 
% 
% fig1=figure('units','centimeters', 'Position',[2 2 35 15]) %18 width 15 heigh
% for i=1:4
%     subaxis(2,4,i,'SpacingVert',0.07,'SpacingHoriz',0.05,'MR',0.05, 'ML',0.05,'MT',0.08,'MB',0.14)
%     plot(time, squeeze(BFi_all(:,i,:))');
%     set(gca,'xticklabel',{})
%     if i==1
%         ylabel('aDb_{1cm}')
%     end
%     
%     subaxis(2,4,i+4)
%     plot(time, squeeze(BFi_all2(:,i,:))');
%     xlabel('Time (s)')
%     
%     if i==1
%         ylabel('aDb_{2cm}')
%     end
%     
% end
% 
% 
% %% Analyze INDIVIDUAL
% % Determine the Blood Flow Index for each Channel
% tau=data_tau;
% 
% for subject=1:9
%     for chan=1:4
%         X=squeeze(BFi_all(subject,chan,:))';
%         BFi_all_scaled(subject,chan,:)=X/mean(X(1,6:10));
%         X=squeeze(BFi_all2(subject,chan,:))';
%         BFi_all_scaled2(subject,chan,:)=X/mean(X(1,6:10));
%     end
% end
% 
% 
% 
% fig1=figure('units','centimeters', 'Position',[2 2 35 15]) %18 width 15 heigh
% for i=1:4
%     subaxis(2,4,i,'SpacingVert',0.07,'SpacingHoriz',0.05,'MR',0.05, 'ML',0.06,'MT',0.08,'MB',0.14)
%     plot(time, squeeze(BFi_all_scaled(:,i,:))');
%     set(gca,'xticklabel',{})
%     if i==1
%         ylabel('aDb_{no cuff}')
%     end
%     rsd=[1,1.5,2, 2.5]
%     tit=strcat('r_{SD}=',num2str(rsd(i)),'cm');
%     title(tit)
%     
%     subaxis(2,4,i+4)
%     plot(time, squeeze(BFi_all_scaled2(:,i,:))');
%     xlabel('Time (s)')
%     
%     if i==1
%         ylabel('aDb_{cuff}')
%     end
%     
% end
% 
% %%
% %% Analyze INDIVIDUAL
% % Determine the Blood Flow Index for each Channel
% tau=data_tau;
% 
% for subject=1:9
%     for chan=1:4
%         X=squeeze(BFi_all(subject,chan,:))';
%         BFi_all_scaled(subject,chan,:)=100*(X/mean(X(1,6:10))-1);
%         X=squeeze(BFi_all2(subject,chan,:))';
%         BFi_all_scaled2(subject,chan,:)=100*(X/mean(X(1,6:10))-1);
%     end
% end
% 
% 
% 
% fig1=figure('units','centimeters', 'Position',[2 2 35 15]) %18 width 15 heigh
% for i=1:4
%     subaxis(2,4,i,'SpacingVert',0.07,'SpacingHoriz',0.05,'MR',0.05, 'ML',0.08,'MT',0.08,'MB',0.14)
%     plot(time, squeeze(BFi_all_scaled(:,i,:))');
%     set(gca,'xticklabel',{},'ylim',[-100 100])
%     if i==1
%         ylabel('aDb_{no cuff}')
%     end
%         rsd=[1,1.5,2, 2.5]
%     tit=strcat('r_{SD}=',num2str(rsd(i)),'cm');
%     title(tit)
%     
%     subaxis(2,4,i+4)
%     plot(time, squeeze(BFi_all_scaled2(:,i,:))');
%     xlabel('Time (s)')
%     set(gca,'ylim',[-100 100])
%     if i==1
%         ylabel('aDb_{cuff}')
%     end
%     
% end
% 
% %% Analyze INDIVIDUAL
% % Determine the Blood Flow Index for each Channel
% tau=data_tau;
% 
% for subject=1:9
%     for chan=1:4
%         X=squeeze(BFi_all(subject,chan,:))';
%         BFi_all_scaled(subject,chan,:)=100*(X/mean(X(1,6:10))-1);
%         X=squeeze(BFi_all2(subject,chan,:))';
%         BFi_all_scaled2(subject,chan,:)=100*(X/mean(X(1,6:10))-1);
%     end
% end
% 
% 
% 
% 
% fig1=figure('units','centimeters', 'Position',[2 2 40 12]) %18 width 15 heigh
% for i=1:4
%     
%     mean_BFi=squeeze(mean(BFi_all_scaled(:,i,:),1));
%     std_BFi=std(squeeze(BFi_all_scaled(:,i,:)));
% 
%     mean_BFi2=squeeze(mean(BFi_all_scaled2(:,i,:),1));
%     std_BFi2=std(squeeze(BFi_all_scaled2(:,i,:)));
%     
%     
%     subaxis(1,4,i,'SpacingVert',0.07,'SpacingHoriz',0.02,'MR',0.02, 'ML',0.06,'MT',0.085,'MB',0.15)
%     rectangle('Position', [30 -61 15  122], 'Facecolor', [0.9 0.9 0.9]);
%     lineProps.col{1} = 'k';
%     lineProps.edgestyle = ':';
%     mseb(time,mean_BFi,std_BFi,lineProps,1)
%     lineProps.col{1} = 'r';
%     lineProps.edgestyle = ':';
%     mseb(time,mean_BFi2,std_BFi2,lineProps,1) 
%     
%     rsd=[1,1.5,2, 2.5]
%     tit=strcat('r_{SD}=',num2str(rsd(i)),'cm');
%     title(tit)
%     
%     set(gca,'ylim',[-60 60],'xlim',[10 60])
%     if i==1
%         ylabel('\DeltaBFi (%)')
%         
%     else
%         set(gca,'yticklabel',{})
%     end
%     xlabel('Time (s)')   
%    
%     if i==4
%         legend('Cuff OFF','Cuff ON')
%     end
%     
% end
% %%
% 
% BFI=squeeze(BFi_all_scaled(:,1,:));
% BFI(11:19,:)=squeeze(BFi_all_scaled(:,4,:));
% 
% BFI2=squeeze(BFi_all_scaled2(:,1,:));
% BFI2(11:19,:)=squeeze(BFi_all_scaled2(:,4,:));
% 
% %%
% tau=data_tau;
% 
% for subject=1:9
%     for chan=1:4
%         X=squeeze(BFi_all(subject,chan,:))';
%         BFi_all_scaled(subject,chan,:)=100*(X/mean(X(1,6:10))-1);
%         X=squeeze(BFi_all2(subject,chan,:))';
%         BFi_all_scaled2(subject,chan,:)=100*(X/mean(X(1,6:10))-1);
%     end
% end
% 
% fig1=figure('units','centimeters', 'Position',[2 2 45 18]) %18 width 15 heigh
% 
% for i=1:4
%     for subject=1:9
%     subaxis(4,9,subject+9*(i-1),'SpacingVert',0.01,'SpacingHoriz',0.02,'MR',0.02, 'ML',0.06,'MT',0.04,'MB',0.15)
%     rectangle('Position', [30 -71 15  142], 'Facecolor', [0.9 0.9 0.9]);
%     hold on
%     plot(time,squeeze(BFi_all_scaled(subject,i,:)),'k','LineWidth',2)
%     plot(time,squeeze(BFi_all_scaled2(subject,i,:)),'r','LineWidth',2)
%     rsd=[1,1.5,2, 2.5]
% 
%     set(gca,'ylim',[-70 70],'xlim',[20 60])
%     if subject==1&&i==1
%         ylabel('\DeltaBFi_{1cm} (%)')
%     elseif subject==1&&i==2
%         ylabel('\DeltaBFi_{1.5cm} (%)')
%     elseif subject==1&&i==3
%         ylabel('\DeltaBFi_{2cm} (%)')
%     elseif subject==1&&i==4
%         ylabel('\DeltaBFi_{2.5cm} (%)')
% %         set(gca,'yticklabel',{})
%     end
%     if subject>1
%         set(gca,'yticklabel',{})
%     end
% %         tit=strcat('r_{SD}=',num2str(rsd(i)),'cm');
% %     title(tit)
%     if i<4
%        set(gca,'xticklabel',{})
%     else
%         xlabel('Time (s)')
%     end
%     
% 
%        legend('OFF','ON')
% 
%     
%     end
%     
% end
% 
% %%
% 
