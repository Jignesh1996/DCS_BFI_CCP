fclose('all'); % Close all open files
clear all; 
close all;
clc;

% load('data_tau.csv')

load("D:\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\TP_Study\DCS\20220608 - 5\Data.mat");
tau = Data_tau;

%% Obtaining the avg g2 curves
g2(1,:,:)=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2_1_temp=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(Data(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(Data(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(Data(:,4,:)-1); %g2-1 curve generation
% average g2 curve for large source detector separation
% for i=1:size(g2,2)
%     g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:)+g2_5_temp(i,:)+g2_6_temp(i,:)+g2_7_temp(i,:)+g2_8_temp(i,:))/8;
% end
for i=1:size(g2,2)
    g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
end

g2_1_new = mean(g2(:,1,:),2);
%% Analyze INDIVIDUAL
% Determine the Blood Flow Index for each Channel
% tau=1e-6*(1:1:40);
tau = Data_tau;

l1 = 0.65; %[cm]
l2 = 0.45; %[cm]

%Optical Properties
% amp = 1000;
ua1 = 0.15; %[cm^-1]
us1 = 8; %[cm^-1]
ua2 = 0.08; %[cm^-1]
us2 = 8; %[cm^-1]
ua3 = 0.2; %[cm^-1]
us3 = 8; %[cm^-1]

%% 

    rho1 = 10; %mm;
    rho2 = 25;
                    
                  
    corr1 = 1+squeeze(g2_1_new(1,:,:));
    corr2 = 1+squeeze(g2_1_new(2,:,:));
    
    tau1 = Data_tau;
    tau2 = Data_tau;
         
    [beta1,Bspt1,Bept1,spt1,ept1] = PointsSC2019(tau1,corr1,rho1);
    [beta2,Bspt2,Bept2,spt2,ept2] = PointsSC2019(tau2,corr2,rho2);
                    
                  
    beta1=corr1(1,1);
    beta2=corr2(1,1);


    %% Multi-Layer MD Model Fitting
    % LB.three = [0 0 0 0];
    % UB.three = [inf 0 inf 2];
    % Starting = [1e-8 0 1e-8 0.1]; %[Fec, Fc, top layer thickness; cm^2/s, cm^2/s, cm]
    
%     if window==1
        LB.three = [0 0 0];
        UB.three = [inf 0 inf ];
        Starting = [5e-9 0 5e-9]; %[Fec, Fc, top layer thickness; cm^2/s, cm^2/s, cm]
%     else
%         LB.three = [0.8*FittedParams(1) 0 0.8*FittedParams(3)];
%         UB.three = [1.2*FittedParams(1) 0 1.2*FittedParams(3)];
%         Starting = [FittedParams(1) 0 FittedParams(3)]; %[Fec, Fc, top layer thickness; cm^2/s, cm^2/s, cm]
%     end
       
        options = optimset('Display','final','TolX',1e-11,'MaxIter', 1000, 'MaxFunEvals', 1000);
        
        FittedParams = fminsearchbnd(@MD_DCS_N2H_BETA,Starting,LB.three,UB.three,options,...
            tau1,tau2,corr1,corr2,ua1,us1,ua2,us2,ua3,us3,rho1,rho2,spt1,spt2,ept1,ept2,l1,l2,beta1-1,beta2-1);
        Fscalp = FittedParams(1);
        Fskull = FittedParams(2);
        Fbrain = FittedParams(3);
        
%         Results_FSCBR(window,1) = Fscalp;
%         Results_FSCBR(window,2) = Fbrain;
        
        
        


%     sub=sub+1;

%%
figure(1)

    plot(Results_FSCBR(:,2)./Results_FSCBR(1,2),'r')
    hold on
    plot(Results_FSCBR(:,1)./Results_FSCBR(1,1)','g')
