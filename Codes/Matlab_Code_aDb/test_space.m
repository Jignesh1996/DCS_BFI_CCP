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
dcs_1res = resample(dcs_1,50,1)
plot(dcs_1res)
hold on
plot(dcs_1up)

%%
xu = (1:1:length(ecg1))/1000;
xd = (1:1:length(dcs_1))/20;

%%
ecg_res = resample(ecg1,1,50)
te = (1:1:length(ecg_res));
plot(xd,ecg_res,'b')
legend("Resampled","Original")
hold on
plot(xu, ecg1,'r')
%%
ini = dcs_1(locs_d(1):locs_d(2));
cyc =zeros(length(pks)-1,length(ini)); 
cyc(1,:)= ini;
x = (1:1:length(ini))/1000;
count = 0;
avg = ini;
for i=1:1:length(locs_d)-1
    count = count+1;
 
    sig = dcs_1(find(xd==round(xu(locs(i)),1)):find(xd==round(xu(locs(i+1)),1)));
    hold on;
    if length(avg) > length(sig)
        sig(length(sig):length(avg)) = 0;
    elseif length(avg)< length(sig)
        sig = sig(1:length(avg));
    end
    avg = avg+sig;
    plot((1:length(sig))/20,sig, 'DisplayName',""+i+"");
%     legend show
    cyc(count,:) = sig;
        
end
xlabel("Time(s)")
ylabel("aDb value")
title("Marker=ECG R peak, DCS 1cm Baseline Forearm")
hold off
avg = avg/count;

figure()
plot(x,avg)

%Plotting the ensemble average
ensavg = mean(cyc,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cyc,[],1)/sqrt(count);                             % Calculate 95% Confidence Intervals         
figure()
plot(x, ensavg, '-r', 'LineWidth',1)
hold on;
plot(x, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
grid
legend('Ensemble Average', '95% Confidence Intervals')

%%resampling

%%

figure()
hold on
plot(xd,dcs_1/max(dcs_1),'r')
plot(xd,ecg_res/max(ecg),'b')

%% Alligning the peaks of the DCS signal
% The idea is to coincide all the cycles by shifting the peaks to coincide
% with eachother to remove the temporal shift.
Extract=ones(size(pks_ECG_smooth,2)-1,min(diff(locs_ECG_smooth)));
Extract=Extract*NaN;

dcs_1_smooth2=circshift(dcs_1_smooth,shift)

for i=1:size(pks_ECG_smooth,2)-1
    locs_ECG_smooth(i)
    locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))
    signal = dcs_1_smooth2(1,locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1);
    base_sig = dcs_1_smooth2(1,locs_ECG_smooth(1):locs_ECG_smooth(1)+min(diff(locs_ECG_smooth))-1);
    sig = signal + (max(base_sig) - max(signal));
    [pks_b,locs_b] = findpeaks(base_sig, 'MinPeakHeight',0.2);
    [pks_s,locs_s] = findpeaks(sig, 'MinPeakHeight',0.2);
%     if locs_b>=locs_s
%         sig = shi%%%%%%%% Start from here. Subtract the difference betweent base signal location from the current signal
   
    if i==5
        break;
    end
    
%     Extract(i,:)=dcs_1_smooth2(1,locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1);

    Extract(i,:)=sig;
end

%%
bp = 1:100:length(dcs_1lp);
dcs_d = detrend(normalize(dcs_1lp),1,bp);
plot(normalize(dcs_1lp));
hold on;
plot(dcs_d);

%% Signal reconstruction using harmonics
step = 50;
final = zeros(1,step);
for i=1:step
    final = final + squeeze(Data(i+800,4,:));
end
final = final/step;
plot(final)

%%
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
signal = ecg1;
L = length(signal);             % Length of signal
t = (0:L-1)*T;  

Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% calculating the inverse fft to get the amplitude of the signal



f = Fs*(0:(L/2))/L;
figure()
plot(f,P1) 
hold on;
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

m = max(P1(100:length(P1)));


[p_hr,l_hr] = findpeaks(P1(100:length(P1)), "MinPeakHeight", m-0.01);
f_hr = Fs*(l_hr+98)/(L);
scatter(f_hr,p_hr);
hold off;

y = Y;
a = abs(ifft(y,40960));
figure();
plot(1:1:length(a),normalize(a))
% figure();
hold on;
plot(1:1:length(a),ecg1(1:length(a)))

%%
signal = bp_a;
L = length(signal);             % Length of signal
t = (0:L-1)*T;  

Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1) 
hold on;
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% m = max(P1(100:length(P1)));
% 
% [p_sig,l_sig] = findpeaks(P1(100:length(P1)), "MinPeakHeight", m-0.01);
% f_sig = Fs*(l_sig+98)/(L);
% scatter(f_sig,p_sig);


y = Y(1:l_hr+150);
% y = Y;
a = abs(ifft(y,4096))/30;
figure();
plot(a)



%%

t_bp = 0:(2*pi*l_bp)/794:2*pi*l_bp;
sin_bp = p_ecg*sin(t_bp);
plot(t_bp,sin_bp);
hold on;
t_sig = 0:(2*pi*l_sig)/794:2*pi*l_sig;
sin_sig = p_sig*sin(t_sig);
plot(t_sig,sin_sig);

figure();
plot(sin_sig,sin_bp)

%%
pBP =  (max(bp_a)+ 2*min(bp_a))/3;  % MAP pressure
psig = mean(signal);   % mean value of signal

CrCP = pBP - (psig/p_sig)*p_ecg;