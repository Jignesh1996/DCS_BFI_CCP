function sig_lpf = lpf_ffilts(signal,fc,fs)
Fs = fs;
fc = fc;
fn = fc/fs;
sbf = fn/0.85;
T = 1/Fs;             % Sampling period    
signal = signal;
L = length(signal);             % Length of signal
t = (0:L-1)*T;  

d = designfilt('lowpassfir', ...
    'PassbandFrequency',fn,'StopbandFrequency',sbf, ...
    'PassbandRipple',1,'StopbandAttenuation',60, ...
    'DesignMethod','equiripple');
sig_lpf = filtfilt(d,signal);

end