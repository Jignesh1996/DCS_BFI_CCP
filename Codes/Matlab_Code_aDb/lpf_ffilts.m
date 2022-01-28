function sig_lpf = lpf_ffilts(signal,fc,fs)
Fs = fs;

fn = fc/fs;
sbf = fn/0.85;
fprintf("%d",sbf)
T = 1/Fs;             % Sampling period    

L = length(signal);             % Length of signal
t = (0:L-1)*T;  

d = designfilt('lowpassfir', ...
    'PassbandFrequency',fn,'StopbandFrequency',sbf, ...
    'PassbandRipple',1,'StopbandAttenuation',120, ...
    'DesignMethod','kaiserwin');
sig_lpf = filtfilt(d,signal);
% [h,t] = impz(d);
% figure();
% plot(t,h);
% freqz(d,1024);

end