%% DCS two_layered fitting function
function sse = MD_DCS_N2H(params,tau1,tau2,g2_raw1,g2_raw2,ua1,us1,ua2,us2,ua3,us3,rho1,rho2,spt1,spt2,ept1,ept2,l1,l2)
F1 = params(1); % Flow of extracerebral tissue cm2/s
F2 = params(2); %Thickness of top layer
F3 = params(3); % Flow of cerebral tissue cm2/s
% l1 = params(4);
% l2 = params(5);
beta = params(4);

%%
rho1 = rho1/10;
rho2 = rho2/10;


%% First Detector
[G1_first,F1,F2,F3] = G1fun(rho1,tau1,ua1,us1,ua2,us2,ua3,us3,l1,l2,F1,F2,F3);
% Field autocorrelation function
g1_first = (G1_first)./max(G1_first);
g2_1_first = beta*(g1_first.^2);
g2_first = g2_1_first + 1;
% Vector to minimize;
Error_Vector1 = sum((g2_first(spt1:ept1)' - g2_raw1(spt1:ept1)).^2);

%% Second Detector
[G1_second,F1,F2,F3] = G1fun(rho2,tau2,ua1,us1,ua2,us2,ua3,us3,l1,l2,F1,F2,F3);
% Field autocorrelation function
g1_second = (G1_second)./max(G1_second);
g2_1_second = beta*(g1_second.^2);
g2_second = g2_1_second + 1;
% Vector to minimize;
Error_Vector2 = sum((g2_second(spt2:ept2)' - g2_raw2(spt2:ept2)).^2);

%% Vector to minimize;
Error_Vector = Error_Vector1 + Error_Vector2;
sse=sum(Error_Vector);

end