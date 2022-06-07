function [beta,Bspt,Bept,spt,ept,g1] = Points(tau,g2,rho)
%% Truncate g2 and tau
spt_g2 = find(tau>1e-7,1,'first');
g2 = g2(spt_g2:end);
tau = tau(spt_g2:end);

% Determine beta
% Bspt = find(tau>1.5e-6,1,'first');
if rho<11
    Bspt = find(tau>1e-6,1,'first');
    Bept = Bspt + 5;
elseif (11<rho) && (rho<25)
    Bspt = find(tau>4e-7,1,'first');
    Bept = Bspt + 5;
else
    Bspt = find(tau>4e-7,1,'first');
    Bept = Bspt + 5;
end

beta = mean(g2(Bspt:Bept))-1;

% Determine spt and ept of fit
g1 = abs(sqrt((g2-1)/beta));

spt = find(tau>4e-7,1,'first');
ept = length(tau);
% ept = find(tau>1e-3,1,'first');
% ept = find(g1<0.2,1,'first'); % Baker

spt = spt + spt_g2 - 1;
ept = ept + spt_g2 - 1;
Bspt = Bspt + spt_g2 - 1;
Bept = Bept + spt_g2 - 1;