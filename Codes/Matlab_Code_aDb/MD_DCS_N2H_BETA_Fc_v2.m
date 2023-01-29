%% DCS two_layered fitting function
function sse = MD_DCS_N2H_beta_Fc_v2(params,tau1,tau2,g2_raw1,g2_raw2,ua1,us1,ua2,us2,ua3,us3,rho1,rho2,beta1,beta2,dFec)
F1 = params(1); % Flow of extracerebral tissue cm2/s
F2 = 0.01*F1; % Flow of extracerebral tissue cm2/s
F3 = params(3); % Flow of cerebral tissue cm2/s
l1 = params(4);
l2 = params(5);
% beta = params(4);

%%
rho1 = rho1/10;
rho2 = rho2/10;
% rho3= rho3/10;

%% First Detector
[G1_first,F1,F2,F3] = G1fun(rho1,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,F1,F2,F3);

% Field autocorrelation function
g1_first = (G1_first)./max(G1_first);
g2_1_first = beta1*(g1_first.^2);
g2_first = g2_1_first + 1;
% Vector to minimize;
Error_Vector1 = sum((g2_first - g2_raw1).^2);

% 
% figure(1)
% 
% semilogx(tau1, g2_raw1,'k')
% hold on
% semilogx(tau1, g2_first,'r')


%% Second Detector

[G1_second,~,~,~] = G1fun(rho2,tau2,ua1,us1,ua2,us2,ua3,us3,l1,l2,F1*dFec,F2,F3);
% [G1_second,~,~,~] = G1fun(rho2,tau2,ua1,us1,ua2,us2,ua3,us3,l1,l2,F1,F2,F3);
% Field autocorrelation function
g1_second = (G1_second)./max(G1_second);
g2_1_second = beta2*(g1_second.^2);
g2_second = g2_1_second + 1;
% Vector to minimize;
Error_Vector2 = sum((g2_second - g2_raw2).^2);
% 
% % 
% subplot(1,2,2)
% semilogx(tau1, g2_raw2,'k')
% hold on
% semilogx(tau1, g2_second,'r')

%% Third Detector

% [G1_third,F1,F2,F3] = G1fun(rho3,tau3,ua1,us1,ua2,us2,ua3,us3,l1,l2,F1,F2,F3);
% % Field autocorrelation function
% g1_third = (G1_third)./max(G1_third);
% g2_1_third = beta3*(g1_third.^2);
% g2_third = g2_1_third + 1;
% % Vector to minimize;
% Error_Vector3 = sum((g2_third(spt2:ept2)' - g2_raw3(spt2:ept2)).^2);

%% Vector to minimize;
Error_Vector = Error_Vector1 + Error_Vector2;
% Error_Vector = Error_Vector1;
% Error_Vector = Error_Vector1 + Error_Vector2 + Error_Vector3;
sse=sum(Error_Vector);

end