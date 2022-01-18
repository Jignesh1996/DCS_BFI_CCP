% Checking the upsampling factor for the interp and interp1 functions, if
% it is linear or not
y = interp(dcs_1,50);
x = 1:length(y);
x_d = (1:length(dcs_1));
[pks_u,locs_u] = findpeaks(y, 'MinPeakHeight', 0.5,'MinPeakDist',600,'MinPeakProminence',0.1);  %Determine peaks and Indices
figure()
[pks_d,locs_d] = findpeaks(dcs_1, 'MinPeakHeight', 0.5,'MinPeakProminence',10);
[pks_u1,locs_u1] = findpeaks(dcs_1up, 'MinPeakHeight', 0.5,'MinPeakDist',600,'MinPeakProminence',0.1);
a = locs_u./locs_d;
b = locs_u1./locs_d;
plot(x_d,dcs_1)
hold on
plot(x_d(locs_d),pks_d, '+r')
hold off

%% 

