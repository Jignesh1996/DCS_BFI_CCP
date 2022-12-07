function [G1,F1,F2] = G1_fun_2_layer(rho,tau,ua,us,l1,F)

% Inputting the data for this function is as follows
%us = [us(ec) us(c)]; us = [ua(ec) ua(c)]; F = [Fec Fc];
N = 5000;
dp = 20;
p = (0:N-1)/dp;
ds = dp/N;
s = pi*(0:N-1)*(ds); % [cm^-1] %radial spatial frequency, cycles/unit-time

n = 1.4; % refractive index of the medium
lambda = 786.5e-7; % [cm] wavelength of the light
Reff = 0.493; % The effective reflection coefficient

ua1 = ua(1);
ua2 = ua(2);
us1 = us(1);
us2 = us(2);
F1 = F(1);
F2 = F(2);

parfor i=1:length(tau)
    s0 = 1;
    k0 = (2*pi*n)/lambda; % [cm^-1] Wavenumber of light in the medium %assume the same wave number for each layer
    
    D1 = 1/(3*us1 + 3*ua1); % [cm] diffusion coefficient for 1st layer
    D2 = 1/(3*us2 + 3*ua2); % [cm] diffusion coefficient for 2nd layer
    
    z0 = 1/(ua1+us1); % z0 = 3*D1; % [cm] first point source term
    zb = 2*D1*(1+Reff)/(1-Reff); % z' [cm] 2nd point source term
    
    k1 = sqrt((D1.*s.^2 + ua1 +2*us1*(k0^2)*F1*tau(i))./D1);
    k2 = sqrt((D2.*s.^2 + ua2 +2*us2*(k0^2)*F2*tau(i))./D2);
    a1 = k1;
    a2 = k2;
%     a1 = sqrt(3*ua1*us1 + 6*us1^2*k0^2*F1*tau(i) + s.^2); % [cm^-1]
%     a2 = sqrt(3*ua2*us2 + 6*us2^2*k0^2*F2*tau(i) + s.^2); % [cm^-1]
    
    %% Hyperbolic Functions

    part1 = a1.*D1.*cosh(a1.*l1)+ a2.*D2.*sinh(a1.*l1);
    part2 = a1.*D1.*cosh(a1.*(l1+zb))+ a2.*D2.*sinh(a1.*(l1+zb));

    part3 = sinh(a1.*(zb+z0));
    part4 = D1.*a1;
    part5 = sinh(a1.*z0);
    part6 = D1.*a1;

    G1_f = ((part1.*part3)./(part2.*part4)) - (part5./part6);
    
    
    G1(i) = sum(ds.*G1_f(:).*s(:).*besselj(0,rho.*s(:)))/((2*pi)); 
end