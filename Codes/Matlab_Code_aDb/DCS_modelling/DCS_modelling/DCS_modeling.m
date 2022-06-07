clear all
close all
%%
lambda = 786e-7; % Wavelenght in cm
rho=3; % Source detetcor spearation
tau=load('data_tau.csv');
beta=0.5;

%% 3 layered model
LayNo=3;

L=[0.5 0.5]; % First and second layer thicknes in cm
mua=[0.01 0.01 0.01];
mus=[10 10 10];
F=[1e-10 1e-10 1e-10]; 

[G1] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta);
%% 2 layered model
LayNo=2;

L=1; % First layer thicknes in cm
mua=[0.01 0.01];
mus=[10 10];
F=[1e-10 1e-9];

[G1a] = G1fun2022(rho,tau,mua,mus,L,F,lambda,LayNo,beta);


%% Figure
semilogx(tau,G1)
hold on
semilogx(tau,G1a,'r')

%% 2 layered model
LayNo=2;
rho=1;
L=1.5; % First layer thicknes in cm
mua=[0.01 0.01];
mus=[10 10];
F=[1e-9 1e-9];

[G1(1,:),g2] = G1fun2022(rho,tau,mua,mus,L,F, lambda,LayNo,beta);

%%
F=[0.5*1e-9 1e-9];
[G1(2,:)] = G1fun2022(rho,tau,mua,mus,L,F, lambda,LayNo);
g1(1,:)=G1(2,:)/G1(2,1);
g2(2,:)=beta*(g1.^2);


F=[1e-9 1e-9];
rho=3;
[G1a(1,:)] = G1fun2022(rho,tau,mua,mus,L,F, lambda,LayNo);
g1a(1,:)=G1a(1,:)/G1a(1,1);
g2a(1,:)=beta*(g1a.^2);

F=[0.5*1e-9 1e-9];
[G1a(2,:)] = G1fun2022(rho,tau,mua,mus,L,F, lambda,LayNo);
g1a(1,:)=G1a(2,:)/G1a(2,1);
g2a(2,:)=beta*(g1a.^2);


semilogx(tau,G1(1,:),'r')
hold on
semilogx(tau,G1(2,:))

%%

DeltaOD_1=-log(g2(2,:)./g2(1,:));
DeltaOD_3=-log(g2a(2,:)./g2a(1,:));

DeltaFec(1,:)=DeltaOD_1/(0.3*1e-9);
DeltaFec(2,:)=DeltaOD_3/(0.3*1e-9);


semilogx(tau,DeltaFec(1,:),'r')
hold on
semilogx(tau,DeltaFec(2,:),'b')