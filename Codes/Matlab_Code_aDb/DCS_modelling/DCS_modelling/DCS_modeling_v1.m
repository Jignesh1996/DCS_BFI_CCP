clear all
close all
%%
lambda = 786e-7; % Wavelenght in cm
rho=3; % Source detetcor spearation
tau=load('data_tau.csv');
tau = tau;
beta=Data(1,1,1)-1;

%% 3 layered model
% LayNo=3;
% 
% L=[0.5 0.5]; % First and second layer thicknes in cm
% mua=[0.17 0.17 0.17];
% mus=[10 10 10];
% % F=[1e-10 1e-10 1e-10]; 
% F=[1e-10 1e-10 1e-10]; 
% 
% [G1] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta);

%%



%% 2 layered model
LayNo=2;

L=1; % First layer thicknes in cm
mua=[0.17 0.17];
mus=[10 10];
% F=[1e-10 1e-9];
F = [mean(aDb1(1,1:600)) mean(aDb1(2,1:600))];

[G1a] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta);


%% Figure
% semilogx(tau,G1)
% hold on
semilogx(tau,G1a,'r')

%% 2 layered model
LayNo=2;
rho=1;
L=1.5; % First layer thicknes in cm
mua=[0.17 0.17];
mus=[10 10];
F=[1e-8 1e-8];
% F = aDb1(1,1:40);
% F = [mean(aDb1(1,1:600)) mean(aDb1(2,1:600))];

[G1(1,:),g2] = G1fun2022(rho,tau,mua,mus,L,F, lambda,LayNo,beta);
g1(1,:)=G1(1,:)/G1(1,1);
g2(1,:)=beta*(g1.^2);


F=[0.5*1e-8 1e-8];
% F = [mean(aDb1(1,600:1200)),mean(aDb1(2,600:1200))];
[G1(2,:)] = G1fun2022(rho,tau,mua,mus,L,F, lambda,LayNo, beta);
g1(1,:)=G1(2,:)/G1(2,1);
g2(2,:)=beta*(g1.^2);


F=[1e-8 1e-8];
% F =aDb1(2,1:40);
% F = [mean(aDb1(1,1:600)) mean(aDb1(2,1:600))];
rho=3;
[G1a(1,:)] = G1fun2022(rho,tau,mua,mus,L,F, lambda,LayNo,beta);
g1a(1,:)=G1a(1,:)/G1a(1,1);
g2a(1,:)=beta*(g1a.^2);

% F = [mean(aDb1(1,600:1200)),mean(aDb1(2,600:1200))];
F=[0.5*1e-8 1e-8];
[G1a(2,:)] = G1fun2022(rho,tau,mua,mus,L,F, lambda,LayNo,beta);
g1a(1,:)=G1a(2,:)/G1a(2,1);
g2a(2,:)=beta*(g1a.^2);


semilogx(tau,g2(1,:),'r')
hold on
semilogx(tau,g2(2,:),'b')

%% 

aDb1 = standalone_dcs(Data,Data_tau);
%%
% beta=0.5;
% Data(:,1,:) = g2_new(1,:,:);
% Data(:,2,:) = g2_new(2,:,:);
% Data(:,3,:) = g2_new(3,:,:);
% Data(:,4,:) = g2_new(4,:,:);

LayNo=2;
rho=1;
rho2=2.5;
L=1.5; % First layer thicknes in cm
mua=[0.17 0.17];
mus=[10 10];
bsl_span = 1:500;
p_span = 650:1150;
% F=[0.2e-7 1e-7];
F = [mean(aDb1(1,bsl_span)) mean(aDb1(2,bsl_span))];
beta = Data(bsl_span(1),1,1)-1;
[G1(1,:),g2(1,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta);
beta = Data(bsl_span(1),2,1)-1;
[G2(2,:),g2(3,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta);

% F=[0.1e-7 1e-7];
F = [mean(aDb1(1,p_span)),mean(aDb1(2,p_span))];
beta = Data(p_span(1),1,1)-1;
    [G1(3,:),g2(2,:)] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta);
beta = Data(p_span(1),2,1)-1;
[G2(4,:),g2(4,:)] = G1fun2022(rho2,tau,mua,mus,L,F,lambda,LayNo,beta);


%%

DeltaOD_1=-log((g2(3,:)-1)./(g2(1,:)-1));
DeltaOD_3=-log((g2(4,:)-1)./(g2(2,:)-1));

DeltaFec(1,:)=DeltaOD_1/(0.3*10^-9);
DeltaFec(2,:)=DeltaOD_3/(0.3*10^-9);

ratio = DeltaFec(2,:)./DeltaFec(1,:);
semilogx(Data_tau,DeltaFec(2,:)./DeltaFec(1,:));
mean(ratio(1:31))

figure()

semilogx(tau,DeltaFec(1,:),'r')
hold on
semilogx(tau,DeltaFec(2,:),'b')

%% 
g2_n(:,1,:) = g2_new(1,:,:);
g2_n(:,2,:) = g2_new(2,:,:);
g2_n(:,3,:) = g2_new(3,:,:);
g2_n(:,4,:) = g2_new(4,:,:);

%% 
=