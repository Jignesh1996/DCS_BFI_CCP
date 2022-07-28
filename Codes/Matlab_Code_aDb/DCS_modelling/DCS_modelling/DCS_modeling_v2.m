clear all
close all
%%
lambda = 786e-7; % Wavelenght in cm
tau=load('data_tau.csv');
beta=0.5;


%% 2 layered model
beta=0.5;
LayNo=2;
rho=1;
rho2=2.5;
L=1.5; % First layer thicknes in cm
mua=[0.1 0.1];
mus=[10 10];
F=[1.8284e-08 4.9120e-09];

[G1(1,:),g2(1,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta);
[G2(3,:),g2(3,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta);

F=[1.1198e-08 3.3309e-09];
[G1(2,:),g2(2,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta);
[G2(4,:),g2(4,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta);



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
rsd=[1 2.5];

mua = 0.1; %cm^-1 baseline absorption coefficient
mus = 10; %cm^-1 baseline reduced scattering coefficient

for i=1:4
    if i<3
        rho=rsd(1);
    else
        rho=rsd(2);
    end    
        LB = [0.1e-12];
        UB = [10e-7];
        
        Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
        beta= g2(i,1); %0.1568;
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau,g2(i,:),mua,mus,rho,beta);
        aDb1_f(i) = FittedParams(1);
end
%%

figure()

bar(1,100*(aDb1(2)/aDb1(1)-1),'r')
hold on
bar(2,100*(aDb1(4)/aDb1(3)-1),'k')
ylabel('\DeltaBFi (%)')
xlabel('Tau (s)')
% legend('1cm','1cm_{red}', '2.5cm','2.5cm_{red}','Location','SouthWest')


%%

