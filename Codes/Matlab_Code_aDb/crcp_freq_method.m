function crcp_freq = crcp_freq_method(ecg,tcd,bp_a,dcs)

%%% This function calculates the value of critical closing pressure using
%%% the blood pressure and the DCS signal or TCD signal. It uses regression
%%% method to find the intercept of the dynamic pressure-flow graph between
%%% the DCS/TCD and ABP graph.
%%% This fuction uses a harmonics of the signal and then reconstructing the
%%% sinusoidal signal with the fundamental harmonics frequency and using
%%% that to calculate the critical closing pressure. 



%%%% Written by Jignesh Mistry %%%%


    
    

    ecg1 = ecg;

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
    signal = dcs;
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
    a = 0:0.05:1;
    plot(a,normalize(dcs(2:22)));
    hold on;
    plot(t,normalize(sin_dcs))
    title("Fundamental DCS")
    
    sin_tcd =pTCD+ pTCD_f*sin(2*pi*f_sig*t);
    figure()
    x = 0:0.001:1;
    plot(x,normalize(tcd(100:1100)));
    hold on;
    plot(t,normalize(sin_tcd))
    title("Fundamental TCD")
    
    pBP =  (max(bp_a)+ 2*min(bp_a))/3;  % MAP pressure 
    sin_bp =pBP+ pABP_f*sin(2*pi*f_sig*t);
    figure();
    plot(x,normalize(bp_a(100:1100)));
    hold on;
    plot(t,normalize(sin_bp))
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
    
    crcp_freq = [ccp_tcd;ccp_dcs];
end