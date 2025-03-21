function sig_lpf = lpf(signal,fc,fs)
Fs = fs;
fc = fc;

T = 1/Fs;             % Sampling period    
signal = signal;
L = length(signal);             % Length of signal
t = (0:L-1)*T;  

[sig_lpf,d] = lowpass(signal,fc,Fs,'ImpulseResponse','fir','Steepness',0.85,'StopbandAttenuation',120);
% [h,t] = impz(d)
% figure();
% plot(t,h);
% freqz(d,1024);

end
