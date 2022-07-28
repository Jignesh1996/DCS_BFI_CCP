function [G1,F1,F2,F3] = G1fun(rho,tau,ua1,us1,ua2,us2,ua3,us3,l1,l2,F1,F2,F3)
N = 1000;
dp = 20;
p = (0:N-1)/dp;
ds = dp/N;
s = pi*(0:N-1)*(ds); % [cm^-1] %radial spatial frequency, cycles/unit-time

%% Define the parameters used in the correlation function
n = 1.4; % refractive index of the medium
lambda = 786.5e-7; % [cm] wavelength of the light
Reff = 0.493; % The effective reflection coefficient

s0 = 1;
k0 = (2*pi*n)/lambda; % [cm^-1] Wavenumber of light in the medium %assume the same wave number for each layer
    
D1 = 1/(3*us1 + 3*ua1); % [cm] diffusion coefficient for 1st layer
D2 = 1/(3*us2 + 3*ua3); % [cm] diffusion coefficient for 2nd layer
D3 = 1/(3*us3 + 3*ua3); % [cm] diffusion coefficient for 3rd layer
    
z0 = 1/(ua1+us1); % z0 = 3*D1; % [cm] first point source term
zb = 2*D1*(1+Reff)/(1-Reff); % z' [cm] 2nd point source term


parfor i=1:length(tau)
    
    a1 = sqrt(3*ua1*us1 + 6*us1^2*k0^2*F1*tau(i) + s.^2); % [cm^-1]
    a2 = sqrt(3*ua2*us2 + 6*us2^2*k0^2*F2*tau(i) + s.^2); % [cm^-1]
    a3 = sqrt(3*ua3*us3 + 6*us3^2*k0^2*F3*tau(i) + s.^2); % [cm^-1]
    
    %% Exponential Functions
    
    A = exp(a1*(l1-zb));
    B = exp(-a1*(l1-zb));
    C = exp(a2*l2);
    D = exp(-a2*l2);
    E = exp(a1*l1);
    F = exp(-a1*l1);
    
    
    Exp1 = (A + B)/2;
    Exp2 = (C + D)/2;
    Exp3 = (C - D)/2;
    Exp4 = (A - B)/2;
    Exp7 = (E + F)/2;
    Exp8 = (E - F)/2;
    
    Part1 = a1.*D1.*Exp1.*(a2.*D2.*Exp2 + a3.*D3.*Exp3);
    Part2 = a2.*D2.*Exp4.*(a3.*D3.*Exp2 + a2.*D2.*Exp3);
    Part3 = a1.*Exp7.*(D1+a3*D3*z0)
    Part4 = Exp8.*(a3*D3 + a1.^2*D1*z0);
    Part5 = a1.*Exp7.*(a3*D1*D3 + a2.^2*D2^2*z0);
    Part6 = Exp8.*(a2.^2*D2^2 + a1.^2.*a3*D1*D3*z0);
    
    num = s0*z0*(Part1 + Part2);
    denom = a2.*D2.*Exp2.*(Part3+Part4) + Exp3.*(Part5+Part6);
    G1_f = num./denom;
    
    G1(i) = sum(ds*exp(-a1(:)*z0).*G1_f(:).*s(:).*besselj(0,rho.*s(:)))/((2*pi)^2); %% Is it the same weighting factor as previously?
end