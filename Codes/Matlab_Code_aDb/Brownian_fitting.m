%% DCS fitting function
function sse = Brownian_fitting(params,tau,g2_1_raw,ua,us,rho,beta) %for one param add 'beta' after rho
aDb = params(1);
% beta = params(2); %remove this when calculating for one param


%% Define the parameters used in the correlation function
v=3;
zo = 1/(us); %first point source term in cm
% zo = 1/(ua+us); %first point source term in cm
D = 1/(3*(us)); %diffusion coefficient in cm
% D = 1/(3*(ua+us));%diffusion coefficient in cm
zb = 2*D*(1+0.493)/(1-0.493); %2nd point source term in cm
r1 = sqrt(rho^2+zo^2);
r2 = sqrt((zo+2*zb)^2+(rho^2));
n = 1.4;%refractive index of the medium
c=v/n;
% lamda = 780e-7; %786.5e-007;%wavelength of the light in cm
lamda = 850e-7; %786.5e-007;%wavelength of the light in cm
k = (2*pi*n)/lamda;%Wavenumber of light in the medium
k_D = sqrt((3*us*ua)+(6*us^2*k^2*aDb*tau));
Amp = (3*us)/(4*pi);
% beta = 0.15;
%% Field autocorrelation function
G1 = Amp*((exp(-k_D*r1)./r1)-(exp(-k_D*r2)./r2));
g1 = (G1)./max(G1);
g2_1_fit = beta*(g1.^2);
% plot(g2_1_fit);

%% Vector to minimize;
Error_Vector = g2_1_fit - g2_1_raw;
sse=sum(Error_Vector.^2);
end