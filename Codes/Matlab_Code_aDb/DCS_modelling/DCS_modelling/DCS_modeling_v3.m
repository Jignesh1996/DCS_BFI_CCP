clear all
close all
%%
lambda = 786e-7; % Wavelenght in cm
tau=load('data_tau.csv');
beta=0.4;

% tau=1e-6:10e-6:5e-3;

%% 2 layered model
beta=0.5;
LayNo=2;
rho=1;
rho2=2.5;
L=1.5; % First layer thicknes in cm
mua=[0.17 0.17];
mus=[8 8];
% F=[0.2e-7 1e-7];
F = [mean(aDb1(1,1:600)) mean(aDb1(2,1:600))];

[G1(1,:),g2(1,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta);
[G2(2,:),g2(3,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta);

% F=[0.1e-7 1e-7];
F = [mean(aDb1(1,600:1200)),mean(aDb1(2,600:1200))];
[G1(3,:),g2(2,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta);
[G2(4,:),g2(4,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta);


%%
figure()

semilogx(tau,g2(1,:),'r')
hold on
semilogx(tau,g2(2,:),'--r')
semilogx(tau,g2(3,:),'k')
semilogx(tau,g2(4,:),'--k')
ylabel('g_2')
xlabel('Tau (s)')
legend('1cm','1cm_{red}', '2.5cm','2.5cm_{red}','Location','SouthWest')

%%
rsd=[rho rho2];
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10; %cm^-1 baseline reduced scattering coefficient

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
%%

figure()

bar(1,100*(aDb(2)/aDb(1)-1),'r')
hold on
bar(2,100*(aDb(4)/aDb(3)-1),'k')
ylabel('\DeltaBFi (%)')
xlabel('Tau (s)')
% legend('1cm','1cm_{red}', '2.5cm','2.5cm_{red}','Location','SouthWest')


%%

