function [beta,Bspt,Bept,spt,ept,g1] = PointsSC2019(tau,g2,rho)
%% Truncate g2 and tau
spt_g2 = find(tau>1e-7,1,'first');
g2 = g2(1:end);
tau = tau(1:end);

% Determine beta

Bspt = tau(1);
Bept = tau(end);

beta = (g2(1));

% Determine spt and ept of fit
g1 = abs(sqrt((g2-1)/beta));

spt = 1;
ept = length(tau);
% ept = find(tau>1e-3,1,'first');
% ept = find(g1<0.2,1,'first'); % Baker

spt = spt;
ept = ept;
Bspt = Bspt;
Bept = Bept;