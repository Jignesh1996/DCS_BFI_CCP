%% DCS two_layered fitting function
function sse = MD_DCS_N2H_BETA(params,tau1,tau2,g2_raw1,g2_raw2,ua,us,rho1,rho2,spt1,spt2,ept1,ept2,l1,beta1,beta2)
F1 = params(1); % Flow of extracerebral tissue cm2/s
F2 = params(2); % Flow of cerebral tissue cm2/s


%%
rho1 = rho1/10;
rho2 = rho2/10;
F = [F1 F2];

%% First Detector
[G1_first,F1,F2] = G1_fun_2_layer(rho1,tau1,ua,us,l1,F);

% Field autocorrelation function
g1_first = (G1_first)./G1_first(1);
g2_1_first = beta1*(g1_first.^2);
g2_first = g2_1_first + 1;
% Vector to minimize;
Error_Vector1 = sum((g2_first(spt1:ept1)' - g2_raw1(spt1:ept1)).^2);

%% Second Detector
[G1_second,F1,F2] = G1_fun_2_layer(rho2,tau2,ua,us,l1,F);
% Field autocorrelation function
g1_second = (G1_second)./G1_second(1);
g2_1_second = beta2*(g1_second.^2);
g2_second = g2_1_second + 1;
% Vector to minimize;
Error_Vector2 = sum((g2_second(spt2:ept2)' - g2_raw2(spt2:ept2)).^2);


%% Vector to minimize;
Error_Vector = Error_Vector1 + Error_Vector2;
sse=sum(Error_Vector);

end