function ccp = ccp_measure(varargin)

%%% This function calculates the value of critical closing pressure using
%%% the blood pressure and the DCS signal or TCD signal. It uses regression
%%% method to find the intercept of the dynamic pressure-flow graph between
%%% the DCS/TCD and ABP graph.
%%% This fuction uses only the tail of the cycle to calculate the CCP.



%%%% Written by Jignesh Mistry %%%%


% Mandatory arguments
ecg1 = varargin{1};
s = varargin{2};
bp = varargin{3};

% Non-mandatory arguments

if(nargin>3 && ~isempty(varargin{4}) )
    step_size=varargin{4};
else
    step_size = 10;
end

step_size = varargin{4};
bp = bp(1:length(ecg1));
%Upsampling the signal if the signal is dcs
if length(s)<length(ecg1)
%         signal = interp(singal,50);
    x = 1:1:length(s);
    uf = length(ecg1)/length(s);   % Upsampling factor
    xq = (1:(1/uf):length(s)+((uf-1)/uf)); 
    sig_up = interp1(x, s,xq,'makima');
    count  = 0;
else
    sig_up = s(1:length(ecg1));
    count = 1;    
end

% Stacking the DCS/TCD signal with BP to make a single matrix for the processing
d = [sig_up; bp];


for sig=1:length(d(:,1))
    signal = d(sig,:);
   
    time_signal=0.001*(1:1:size(signal,2));
    time_ECG=0.001*(1:1:size(ecg1,2));
    
    figure()
    hold on
    plot(time_signal,signal/max(signal),'r')
    plot(time_ECG,ecg1/max(ecg1),'b')
    
    % Finding the peaks of the signal
    [pks_sig,locs_sig]=findpeaks(signal/max(signal),'MinPeakHeight',0.35);
    [pks_ECG,locs_ECG]=findpeaks(ecg1/max(ecg1),'MinPeakHeight',0.65);
    
    % Filtering the signal
    
    windowsize=5; % how many points you want to use (it will depend on your resolution, we were using 10 so it was 3 s window)
    wages=ones(1,windowsize)/windowsize;
    % wages = window_1;
    sig_smooth=filtfilt(wages,1,signal); % Y is your time course you want to filter, Y_smoth is filtered data
    ecg1_smooth=filtfilt(wages,1,ecg1);
    
    [pks_ECG_smooth,locs_ECG_smooth]=findpeaks(ecg1_smooth./max(ecg1_smooth),'MinPeakHeight',0.65);
    [pks_sig_smooth,locs_sig_smooth]=findpeaks(sig_smooth./max(sig_smooth),'MinPeakHeight',0.35,'MinPeakDistance',500);
    
    % Plotting the signal
    fig1=figure('units','centimeters', 'Position',[2 2 35 13]); %18 width 15 heigh
    hold on
    plot(sig_smooth/max(sig_smooth),'r')
    plot(locs_sig_smooth, pks_sig_smooth,'*k')
    
    plot(ecg1_smooth/max(ecg1_smooth),'b')
    plot(locs_ECG_smooth, pks_ECG_smooth,'*k')
    
    if locs_ECG_smooth(1)>locs_sig_smooth(1)
        Difference=locs_ECG_smooth(1)-locs_sig_smooth(1);
        shift = floor(1.75*Difference);
    else
        shift =0;
    end

    % Extracting the cycle from the signal using R peak as a marker
%     Extract=ones(2,size(pks_ECG_smooth,2)-2,min(diff(locs_ECG_smooth)));
%     Extract=Extract*NaN;
    if sig==2
        dcs_1_smooth2=circshift(sig_smooth,-1075); % Personalizing the shift according to the signal
    elseif count==0
        dcs_1_smooth2=circshift(sig_smooth,775);
    else
        dcs_1_smooth2=circshift(sig_smooth,0);
    end
    
%     for i=1:size(pks_ECG_smooth,2)-2
%         locs_ECG_smooth(i);
%         locs_ECG_smooth(i)+min(diff(locs_ECG_smooth));
%         signal = dcs_1_smooth2(1,locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1);
% %         base_sig = dcs_1_smooth2(1,locs_ECG_smooth(2):locs_ECG_smooth(2)+min(diff(locs_ECG_smooth))-1);
% %         sig = signal + (max(base_sig) - max(signal));
%       
%         
%     %     Extract(i,:)=dcs_1_smooth2(1,locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1);
%     
%         Extract(sig,i,:)=signal;
%     end
    if sig ==1
        Extract=zeros(length(d(:,1)),ceil(size(pks_ECG_smooth,2)/step_size),min(diff(locs_ECG_smooth)));
        Extract=Extract*NaN;
    end

    for i=1:step_size:size(pks_ECG_smooth,2)-2
        sig_a = dcs_1_smooth2(1,locs_ECG_smooth(i):locs_ECG_smooth(i)+min(diff(locs_ECG_smooth))-1);
        cnt = 1;
        for j=1:1:step_size-1
            if i+j>size(pks_ECG_smooth,2)-2
                break;
            end
            locs_ECG_smooth(i+j);
            locs_ECG_smooth(i+j) + min(diff(locs_ECG_smooth));
            cycle = dcs_1_smooth2(1,locs_ECG_smooth(i+j):locs_ECG_smooth(i+j)+min(diff(locs_ECG_smooth))-1);
            sig_a = sig_a+cycle;
            cnt = cnt+1;
        end
        sig_a = sig_a./cnt;
          
    
        Extract(sig,floor(i/step_size)+1,:)=sig_a;
    end

    avg_sig = Extract(sig,:,:);
    avg_sig = reshape(avg_sig,[length(Extract(1,:,1)), length(Extract(1,1,:))]); % Reshaping the signal to make it a 2D matrix

    figure();
    x = (1:1:length(avg_sig'))/1000;
    plot(x,(avg_sig'))
    xlabel("Time(s)")
    ylabel("aDb value")
    title("Marker=ECG R peak, DCS 2.5cm")
    
    % Finding the ensemble average of the signal
    if sig==1
        avg = zeros(2,length(avg_sig(1,:)));
    end 

    if sig==1 && count==0
        ttle = 'DCS 2.5cm Ensemble Avg';
    elseif sig==1 && count==1
        ttle = 'TCD Ensemble Avg';
    else
        ttle = 'Brachial BP Ensemble Avg';
    end

    avg(sig,:) = ens_avg(avg_sig,ttle);
    

end



% Code the calculation of the CCP using the polyfit function
bp_stack = reshape(Extract(2,:,:),[length(Extract(1,:,1)), length(Extract(1,1,:))]);
sig_stack = reshape(Extract(1,:,:),[length(Extract(1,:,1)), length(Extract(1,1,:))]);

[pks_signal,locs_signal] = findpeaks(sig_stack);
% disp(bp_stack)
% plot(bp_stack)
ccp = zeros(1,length(bp_stack(:,1)));
for i=1:length(bp_stack(:,1))
    p = polyfit(bp_stack(i,:)',sig_stack(i,:)',1);
    x = -20:0.005:130;
    f = polyval(p,x);
    
    ccp(i) = mean(x(round(f,2)==0));
%     fprintf("%d",ccp(i));
end

end
