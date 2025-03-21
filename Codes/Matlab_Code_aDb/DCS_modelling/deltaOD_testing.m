fclose('all'); % Close all open files
clear all; 
close all;
clc;

%%
close all;

% tau=load('data_tau.csv');
% load('BFi.mat');
for i = 1:50
l1=0.55;
l2=0.45;

ua1 = 0.13; %[cm^-1]
us1 = 9; %[cm^-1]
ua2 = 0.6*0.13; %[cm^-1]
us2 = 9; %[cm^-1]
ua3 = 0.05*i; %[cm^-1]
us3 = 1; %[cm^-1]

rho1 = 10; %mm
rho2 = 25; %mm;

beta=0.45;

% for i=1
Fbrain=2e-8;
Fscalp=0.5e-8;
% Fscalp=Fbrain/6;

% noiseSignal = rand(1,50)/200;
noiseSignal = 0;

[G1_first] = G1fun(rho1/10,tau,ua1,us1,ua2,us2,ua3,us3,l1,l2,Fscalp,0.01*Fscalp,Fbrain);
G1_first(1,:)=G1_first(1,:)/G1_first(1,1);
G1(1,:)=beta*(G1_first(1,:).^2)+noiseSignal;

[G1_first] = G1fun(rho1/10,tau,ua1,us1,ua2,us2,ua3,us3,l1,l2,0.2*Fscalp,0.01*Fscalp,Fbrain);
G1_first(1,:)=G1_first(1,:)/G1_first(1,1);
G1a(1,:)=beta*(G1_first(1,:).^2)+noiseSignal;


beta=0.45;

% noiseSignal = rand(1,50)/20;

noiseSignal =0;

[G1_first] = G1fun(rho2/10,tau,ua1,us1,ua2,us2,ua3,us3,l1,l2,Fscalp,0.01*Fscalp,Fbrain);
G1_first(1,:)=G1_first(1,:)/G1_first(1,1);
G2(1,:)=beta*(G1_first(1,:).^2)+noiseSignal;

[G1_first] = G1fun(rho2/10,tau,ua1,us1,ua2,us2,ua3,us3,l1,l2,0.2*Fscalp,0.01*Fscalp,1*Fbrain);
G1_first(1,:)=G1_first(1,:)/G1_first(1,1);
G2a(1,:)=beta*(G1_first(1,:).^2)+noiseSignal;

DeltaOD_1=-log(G1a(1,:)./G1(1,:));
DeltaOD_3=-log(G2a(1,:)./G2(1,:));
DeltaFec=DeltaOD_3./DeltaOD_1;

ratio(i) = nanmean(DeltaFec(1,1:30))

semilogx(tau,G1)
hold on
semilogx(tau,G1a)
semilogx(tau,G2)
semilogx(tau,G2a)
legend("1cm Baseline","1cm Pressure","2.5cm Baseline","2.5cm Pressure")
hold off;

figure();
G1_diff_sim =  G1a-G1;
G2_diff_sim = G2a - G2;
semilogx(Data_tau, G1_diff_sim);
hold on;
semilogx(Data_tau, G2_diff_sim);
legend("Change at short distance", "Change at long distance")
hold off;


end
%%
r=[1 2.5];

for i=1:2
    aDb = [1E-8 1e-8];
    us=10;
    ua=0.1;
    rho=r(i);

    v=3;
    zo = 1/(us+ua); %first point source term in cm
    D = 1/(3*(us+us)); %diffusion coefficient in cm
    zb = 2*D*(1+0.493)/(1-0.493); %2nd point source term in cm
    r1 = sqrt(rho^2+zo^2);
    r2 = sqrt((zo+2*zb)^2+(rho^2));
    n = 1.4;%refractive index of the medium
    c=v/n;
    lamda = 852e-7; %786.5e-007;%wavelength of the light in cm
    k = (2*pi*n)/lamda;%Wavenumber of light in the medium
    k_D = sqrt((3*us*ua)+(6*us^2*k^2*aDb(i)*tau));
    Amp = (3*us)/(4*pi);

    G1 = Amp*((exp(-k_D*r1)./r1)-(exp(-k_D*r2)./r2));
    g1 = (G1)./max(G1);
    g2_1_fit(i,:) = beta*(g1.^2);

end

semilogx(tau,g2_1_fit','DisplayName','g2_1_fit')
%%
% tau1 = tau;
% tau2 = tau;
% 
% beta1=corr1(1,1)-1;
% beta2=corr2(1,1)-1;
% 
% 
%     LB =       [1e-10 0 1e-10 ];
%     UB =       [1e-7  0 1e-7  ];
%     Starting = [1e-9  0 1e-9  ]; %[Fec, Fc, top layer thickness; cm^2/s, cm^2/s, cm]
% 
%     options = optimset('Display','final','TolX',1e-10,'MaxIter', 2000, 'MaxFunEvals', 2000);
%         
%     FittedParams = fminsearchbnd(@MD_DCS_N2H_BETA_FcFec_2023,Starting,LB,UB,options,...
%             tau1,tau2,corr1,corr2,ua1,us1,ua2,us2,ua3,us3,rho1,rho2,l1,l2,beta1,beta2,DeltaFec(subject));
% 
% %     LB =       [1e-10 0 1e-10 0.6 0.6];
% %     UB =       [1e-7  0 1e-7  0.8 0.6];
% %     Starting = [1e-9  0 1e-9  0.6 0.6]; %[Fec, Fc, top layer thickness; cm^2/s, cm^2/s, cm]
% %     
% %         FittedParams = fminsearchbnd(@MD_DCS_N2H_BETA_FcFec_2023_fitd,Starting,LB,UB,options,...
% %             tau1,tau2,corr1,corr2,ua1,us1,ua2,us2,ua3,us3,rho1,rho2,beta1,beta2,DeltaFec(subject));    
% % %           
%     Fscalp(subject) = FittedParams(1);
%     Fskull(subject) = FittedParams(2);
%     Fbrain(subject) = FittedParams(3);
% %     d_fit(1,subject) = FittedParams(4);
% %     d_fit(2,subject) = FittedParams(5);
% % %     DeltaFec(subject) = FittedParams(3);
% close all
% mean(Fbrain)
% mean(Fscalp)
% mean(Fbrain)/mean(Fscalp)
% end
% 
% 
% plot(Fbrain,'r')
% hold on
% plot(Fscalp,'k')
% 
% %%
% round(Fscalp*1e8,2)
% BF(:,1)=ans
% round(Fbrain*1e8,2)
% BF(:,2)=ans
% %%
% Fscalp2=Fscalp'*1e8
% Fbrain2=Fbrain'*1e8
% 
% % %%
% % plot(corr1)
% % hold on
% % plot(corr2)
% 
% %%
% for subject=1:9
%     
% % for subject=1:9
%     [G1_first] = G1fun(rho1/10,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,Fscalp(subject),.01*Fscalp(subject),(Fbrain(subject)-(0.01*Fbrain(subject)/2)));
%     [G1_second] = G1fun(rho1/10,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,Fscalp(subject),.01*Fscalp(subject),(Fbrain(subject)+(0.01*Fbrain(subject)/2)));
%     
%     G_ratio=(G1_first(1,1:30)./G1_second(1,1:30));
%     dFc(subject)=sum((2./(0.01*Fbrain(subject)))*(log(G_ratio)));
% % end
% 
% end
% % 
% % dFc2=dFc'
% % 
% % %%
% % Fbrain(1:9)=1e-8;
% % Fscalp(1:9)=1e-9;
% % l1=mean([4 5 5 6 5 5 5 6 6])/10;
% % l2=mean([8 6 7 7 7 7 8 5 5])/10;
% % 
% % for subject =1:9
% % 
% % for subject=1:9
% %     [G1_first] = G1fun(rho1/10,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,Fscalp(subject),.01*Fscalp(subject),(Fbrain(subject)-(0.01*Fbrain(subject)/2)));
% %     [G1_second] = G1fun(rho1/10,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,Fscalp(subject),.01*Fscalp(subject),(Fbrain(subject)+(0.01*Fbrain(subject)/2)));
% %     
% %     G_ratio=(G1_first(1,1:30)./G1_second(1,1:30));
% %     dFc(subject)=sum((2./(0.01*Fbrain(subject)))*(log(G_ratio)));
% % end
% % 
% % end

%% converting the Data to g1 curves

semilogx(Data_tau, mean(squeeze(g1(1,1:30,:))));
hold on;
semilogx(Data_tau, mean(squeeze(g1(2,1:30,:))));
semilogx(Data_tau, mean(squeeze(g1(1,30:45,:))));
semilogx(Data_tau, mean(squeeze(g1(2,30:45,:))));
legend(string(rho1)+" mm baseline",string(rho2)+" mm baseline",string(rho1)+" mm pressure",string(rho2)+" mm Pressure")

% calculating the difference in the g1 curves at baseline and the pressure

g1_1_diff = mean(squeeze(g1(1,30:45,:))) - mean(squeeze(g1(1,1:30,:)));
g1_25_diff = mean(squeeze(g1(2,30:45,:))) - mean(squeeze(g1(2,1:30,:)));
figure()
semilogx(Data_tau, g1_1_diff)
hold on;
semilogx(Data_tau, g1_25_diff)
