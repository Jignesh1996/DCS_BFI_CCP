close all;
clear all;

tic;
dcs_path = "C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\MPCM Study\DCS\";
param_path =  "C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\MPCM Study\20220310\";
file_list = ["20220218 - 3\"]; %Edit this list
file_name = ["MPCM005"];
body_param_files = ["MPCM005_tcd_20220218.mat"] ;

for j=1:length(file_list)

    % Read files;
    dcs_dir = strcat(dcs_path,file_list(j),"Data.mat") ;     % Directory for DCS data
    
    load(dcs_dir);
    data_mat_all(j,:,:,:) = Data;
    aDb1 = hybrid_dcs(Data,Data_tau);
    
    % time resultion - aqusition time used to aquire data
    figure();
    sgtitle(file_name(j))
    t_res=1; % seconds
    time=t_res*(1:1:size(aDb1,2));
    
    subplot(2,2,1)
    
    plot(time,aDb1(1,:))
    title('{\itr}_{SD}=1 cm')
    % set(gca,'xticklabel',{})
    
    subplot(2,2,2)
    
    plot(time,aDb1(2,:))
    title('{\itr}_{SD}=1.5 cm')
    xlabel('Time (s)')
    
    subplot(2,2,3)
    
    plot(time,aDb1(3,:))
    title('{\itr}_{SD}=2 cm')
    % set(gca,'xticklabel',{})
    
    subplot(2,2,4)
    
    plot(time,aDb1(4,:))
    title('{\itr}_{SD}=2.5 cm')
    xlabel('Time (s)')
    
    % Assigning the channels
    dcs_1cm = aDb1(1,:).*10^9;
    dcs_1lp = lpf_ffilts(dcs_1cm,10,20);
    dcs_15 = aDb1(2,:).*10^9;
    dcs_15lp = lpf_ffilts(dcs_15,10,20);
    dcs_2 = aDb1(3,:).*10^9;
    dcs_2lp = lpf_ffilts(dcs_2,10,20);
    dcs_25 = aDb1(4,:).*10^9;
    dcs_25lp = lpf_ffilts(dcs_25,7,20);
    
    % Recombine the DCS signal
    adb_lp = [dcs_1lp;dcs_15lp;dcs_2lp;dcs_25lp];

    %Loading the ECG,TCD and BP file
    param_file = strcat(param_path,body_param_files(j));
    load(param_file);

    ecg_a = data(datastart(1):dataend(1));
    bp_a = data(datastart(2):dataend(2));
    tcd_a = data(datastart(3):dataend(3));
    
    sig_roi = 60; % Length of signal under consideration
    
    ecg1 = ecg_a(1:sig_roi*1000);
    ecg1 = normalize(ecg1);
    ecg1 = lpf(ecg1,5,1000);
        
    tcd = tcd_a(1:length(ecg1));
    % tcd = normalize(tcd);
    tcd = lpf(tcd,15,1000);
    
    bp_a = bp_a(1:length(ecg1));
    bp_a = lpf(bp_a,3,1000);

    

    % Calculating the FFT for frequency spectrum
    sig = [ecg1];
    p_0 = zeros(1,size(sig,1));
    p_f = zeros(1,size(sig,1));
    for i=1:size(sig,1)
        signal = sig(i,:);
        Fs = 1000; 
        T = 1/Fs; 
        L = length(signal);             % Length of signal
        t = (0:L-1)*T;  
        
        Y = fft(signal);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        P1_ecg = P1;
        figure();
        f = Fs*(0:(L/2))/L;
        plot(f,P1) 
        hold on;
        title('Single-Sided Amplitude Spectrum of X(t)')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
        m = max(P1);
        [p_sig,l_sig] = findpeaks(P1(50:length(P1)), "MinPeakHeight", m-0.0001);
        p_sig = p_sig+P1(1);
        f_sig = Fs*(l_sig+49)/(L);
        scatter(f_sig,p_sig);
        l_hr = l_sig+49;
        
        y = Y(1:l_hr+150);
        % y = Y;
        a = abs(ifft(y(1:200),4096));
        figure();
        plot(a)
    
    end
    
    %Calculating the amplitude at f(hr) for DCS signal
    Fs = 20;            % Sampling frequency                    
    T = 1/Fs;             % Sampling period    
    signal = dcs_1lp(1:sig_roi*Fs);
    L = length(signal);             % Length of signal
    t = (0:L-1)*T;  
    
    Y = fft(signal);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    P1_DCS = P1;
    
    pDCS = P1(1);
    pDCS_f = max(P1((f_sig*L/Fs)-20:(f_sig*L/Fs)+20));
    
    %Calculating the amplitude at f(hr) for ABP signal
    signal = tcd;
    Fs = 1000;   
    T = 1/Fs;  
    L = length(signal);             % Length of signal
    t = (0:L-1)*T;  
    
    Y = fft(signal);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    P1_TCD = P1;
    pTCD = P1(1);
    pTCD_f = P1(l_hr);
    
    %Calculating the amplitude at f(hr) for ABP signal
    signal = bp_a;
    Fs = 1000;   
    T = 1/Fs;  
    L = length(signal);             % Length of signal
    t = (0:L-1)*T;  
    
    Y = fft(signal);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    P1_abp = P1;
    pABP = P1(1);
    pABP_f = P1(l_hr);
    
    %
    t = 0:0.01:1/f_sig;
    
    sin_dcs =pDCS+ pDCS_f*sin(2*pi*f_sig*t);
    figure()
    plot(t,sin_dcs)
    title("Fundamental DCS")
    
    sin_tcd =pTCD+ pTCD_f*sin(2*pi*f_sig*t);
    figure()
    plot(t,sin_tcd)
    title("Fundamental TCD")


    pBP =  (max(bp_a)+ 2*min(bp_a))/3;
    sin_bp =pABP+ pABP_f*sin(2*pi*f_sig*t);
    figure();
    plot(t,sin_bp)
    title("Fundamental BP")
    figure();
    scatter(sin_bp,sin_tcd)
    title("Scatter plot of TCD vs BP")
    
    
    
    p_tcd = robustfit(sin_bp,sin_tcd);
    x= -p_tcd(1)/p_tcd(2):1:max(max(sin_bp,[],2));
    y = p_tcd(1)+p_tcd(2)*x;
    figure();
    plot(x,y);
    hold on; 
    scatter(sin_bp,sin_tcd)
    xlabel("ABP (mmHg)")
    ylabel("CBFV (cm/s)")
    title("Regression TCD vs BP")
    hold off;
    ccp_tcd = -p_tcd(1)/p_tcd(2);   % simplifying the linear equation y=m*x + c for y=0 will lead to this equation 
    fprintf("\nCcCP TCD is = %d\n",ccp_tcd);
    
    p_dcs = robustfit(sin_bp,sin_dcs);
    x= -p_dcs(1)/p_dcs(2):1:max(max(sin_bp,[],2));
    y = p_dcs(1)+p_dcs(2)*x;
    figure();
    plot(x,y);
    xlabel("ABP (mmHg)")
    ylabel("cBFi (aDb)")
    hold on; 
    scatter(sin_bp,sin_dcs)
    title("Regression DCS vs BP")
    hold off;
    ccp_dcs = -p_dcs(1)/p_dcs(2);   % simplifying the linear equation y=m*x + c for y=0 will lead to this equation 
    fprintf("CrCP DCS is = %d\n",ccp_dcs);
    
end