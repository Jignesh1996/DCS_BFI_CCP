clear all
close all
%%
lambda = 786e-7; % Wavelenght in cm
% tau=load('data_tau.csv');
beta=0.5;

% tau=1e-6:10e-6:5e-3;
%%
load("C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\Pressure_modulation\Aleena\Aleena_Probe_pressure_1\Data.mat");
aDb1 = standalone_dcs(Data,Data_tau);
tau = Data_tau;
%% 2 layered model
close all;
beta=0.5;
LayNo=2;
rho=0.7;
rho2=3;
% l=[0.4:0.1:3.5]; % First layer thicknes in cm
% val = [1:0.1:3.5];
val = [1.2];
l = val;
% change = [0.1:0.1:1];
fact= zeros(1,length(val));

% for k = 1:1:length(change)
    for j=1:length(fact)
        mua=[0.17 0.17];
        mus=[8 8];
        L = l(j);
       % plot the sensitivity factor wrt the different reduction in the flow
       % values vs perticular L values. so @ fixed L > flow vs sensitivity
       % factor again with different L values.
        F=[4e-8 2e-9];
        beta1 = 0.45;
        beta2 = 0.1568;
        % F = [mean(aDb1(1,1:600)) mean(aDb1(2,1:600))];
        
        [G1(1,:),g2(1,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta1);
        [G1(2,:),g2(3,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta2);
        
        L = l(j);
        F=[(2e-8) 2e-9];
%         F=[(4e-8)*change(k) 2e-9];
        % F = [mean(aDb1(1,600:900)),mean(aDb1(2,600:900))];
        [G1(3,:),g2(2,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta1);
        [G1(4,:),g2(4,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta2);
    
        rsd=[rho rho2];
        mua = 0.17; %cm^-1 baseline absorption coefficient
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
%% 2 layered model
close all;
beta=0.5;
LayNo=2;
rho=1;
rho2=2.5;
% l=[0.4:0.1:3.5]; % First layer thicknes in cm
val = [1:0.05:10];

fact= zeros(1,length(val));


    for j=1:length(fact)
        mua=[0.17 0.17];
        mus=[10 10];
        L = 1.2;
       % plot the sensitivity factor wrt the different reduction in the flow
       % values vs perticular L values. so @ fixed L > flow vs sensitivity
       % factor again with different L values.
        F=[4e-8 2e-9];
        
        % F = [mean(aDb1(1,1:600)) mean(aDb1(2,1:600))];
        
        [G1(1,:),g2(1,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta);
        [G1(2,:),g2(3,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta);
        
%         L = l(j);
        F=[4e-8/val(j) 2e-9];
        % F = [mean(aDb1(1,600:900)),mean(aDb1(2,600:900))];
        [G1(3,:),g2(2,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta);
        [G1(4,:),g2(4,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta);
    
        rsd=[rho rho2];
        mua = 0.17; %cm^-1 baseline absorption coefficient
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
        fact(j) =sens_fact;
        
    
    
    end
    x = (100-(1./val)*100);
    plt(x,fact,"Flow reduction Vs Sensitivity Factor","Flow reduction (percent)","Sensitivity Factor","");
%     plot(val,fact)
%     legend(string(100/k))


%%

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
rsd=[rho rho2];
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
DeltaFc = DeltaOD_l - ((DeltaOD_lp./DeltaOD_sp).*(DeltaOD_s))
