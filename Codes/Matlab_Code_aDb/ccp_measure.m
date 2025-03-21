function ccp = ccp_measure(varargin)

%%% This function calculates the value of critical closing pressure using
%%% the blood pressure and the DCS signal or TCD signal. It uses regression
%%% method to find the intercept of the dynamic pressure-flow graph between
%%% the DCS/TCD and ABP graph.
%%% This fuction uses a full cycle to calculate the CCP rather than only
%%% the tail.



%%%% Written by Jignesh Mistry %%%%


% Mandatory arguments
ecg = varargin{1};
s = varargin{2};   % Flow signal
bp = varargin{3};

% Non-mandatory arguments

if(nargin>3 && ~isempty(varargin{4}) )
    step_size=varargin{4};
else
    step_size = 10;
end

disp(length(ecg));
%Upsampling the signal if the signal is dcs
if length(s)<length(ecg)
%         signal = interp(singal,50);
    x = 1:1:length(s);
    uf = length(ecg)/length(s);   % Upsampling factor
    xq = (1:(1/uf):length(s)+((uf-1)/uf)); 
    sig_up = interp1(x, s,xq,'makima');  % Upsampling the signal
    count  = 0;
    ecg1 = ecg(1:length(s)*50);
    sig_up = circshift(sig_up,0);
%     disp(length(sig_up))
%     disp(length(ecg))
else
    sig_up = s;
    count = 1;    
    ecg1 = ecg(1:length(s));
end

% shifting the signal to time align both bp and flow signals
bp = bp(1:length(ecg1));
[a_bp,l_bp] = findpeaks(normalize(bp),"MinPeakHeight",1.25,'MinPeakDistance',700) ;
[a_sig,l_dcs] = findpeaks(normalize(sig_up),"MinPeakHeight",1.5,'MinPeakDistance',700);
shift = l_bp(1)-l_dcs(1);
bp = circshift(bp,0);
% disp(length(l_dcs))
% figure();
% plot(normalize(bp),'r')
% hold on;
% plot(l_bp, a_bp,'*k')
% hold on;   
% plot(normalize(sig_up),'b')
% plot(l_dcs, a_sig,'*k')
% Stacking the DCS/TCD signal with BP to make a single matrix for the processing
d = [sig_up; bp];

% disp(d)
figure();

% disp(shift)
plot(normalize(circshift(sig_up,0)));
hold on;
plot(normalize(bp));

for sig=1:length(d(:,1))
    signal = d(sig,:);
   
    time_signal=0.001*(1:1:size(signal,2));
    time_ECG=0.001*(1:1:size(ecg1,2));
    
%     figure()
%     hold on
%     plot(time_signal,signal/max(signal),'r')
%     plot(time_ECG,ecg1/max(ecg1),'b')
    
 
    % Filtering the signal
    
    windowsize=5; % how many points you want to use (it will depend on your resolution, we were using 10 so it was 3 s window)
    wages=ones(1,windowsize)/windowsize;
    % wages = window_1;
%     sig_smooth=filtfilt(wages,1,signal); % Y is your time course you want to filter, Y_smoth is filtered data
%     ecg1_smooth=filtfilt(wages,1,ecg1);
    sig_smooth = signal;
    ecg1_smooth = ecg1;
    
    [pks_ECG_smooth,locs_ECG_smooth]=findpeaks(normalize(ecg1_smooth),'MinPeakHeight',1.5,'MinPeakDistance',300);
    [pks_sig_smooth,locs_sig_smooth]=findpeaks(normalize(sig_smooth),'MinPeakHeight',0.35,'MinPeakDistance',300);

    if length(locs_ECG_smooth)<step_size
        step_size = length(locs_ECG_smooth)-1;
    end
    disp(length(pks_ECG_smooth))
    % Plotting the signal
    fig1=figure('units','centimeters', 'Position',[2 2 35 13]); %18 width 15 heigh
    hold on
    plot(normalize(sig_smooth),'r')
    plot(locs_sig_smooth, pks_sig_smooth,'*k')
    
    plot(normalize(ecg1_smooth),'b')
    plot(locs_ECG_smooth, pks_ECG_smooth,'*k')
    

    % Change this condition, this is not right, becuase irrespective the
    % difference may not be there in every case
    % Make sure that the flow signal and BP signal both are peak alligned
%     if locs_ECG_smooth(1)>locs_sig_smooth(1)
%         Difference=locs_ECG_smooth(1)-locs_sig_smooth(1);
%         shift = floor(1.75*Difference);
%     else
%         shift =0;
%     end

    % Extracting the cycle from the signal using R peak as a marker
%     Extract=ones(2,size(pks_ECG_smooth,2)-2,min(diff(locs_ECG_smooth)));
%     Extract=Extract*NaN;
    if sig==2
        dcs_1_smooth2=circshift(sig_smooth,0); % Personalizing the shift according to the signal
    elseif count==0
        dcs_1_smooth2=circshift(sig_smooth,0);
    else
        dcs_1_smooth2=circshift(sig_smooth,0);
    end
 
% Defining the size of the variable  
    if sig ==1
        Extract=zeros(length(d(:,1)),floor(size(pks_ECG_smooth,2)/step_size)-3,min(diff(locs_ECG_smooth)));
        Extract=Extract*NaN;
    end
    disp(floor(size(pks_ECG_smooth,2)/step_size)-2)
    ind = 1;

    % Taking the average of the signal with average of step_size number of
    % cycles to remove the possible occurance of ny
    for i=2:step_size:size(pks_ECG_smooth,2)-2
%         if size(pks_ECG_smooth,2)-i <step_size
%             break
%         end
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
          
    
        Extract(sig,ind,:)=sig_a;
        ind = ind+1;
    end
%     disp(Extract)
% Reshaping the signal to remove the extra dimension
    avg_sig = squeeze(Extract(sig,:,:));

    figure();
    x = (1:1:length(avg_sig'))/1000;
    plot(x,(avg_sig'))
    legend("1","2","3","4")
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

bp_stack = squeeze(Extract(2,:,:));
sig_stack = squeeze(Extract(1,:,:));
% disp(sig_stack);
% plot(bp_stack);
% disp(size(sig_stack))
stack = [bp_stack;sig_stack];
save("ccp_var_stack.mat","stack");
ccp = zeros(1,length(bp_stack(:,1)));
for i=1:length(bp_stack(:,1))
%     p = polyfit(bp_stack(i,:)',sig_stack(i,:)',1);
%     x = -20:0.005:130;
%     f = polyval(p,x);
%     
%     ccp(i) = mean(x(round(f,2)==0));
    p = robustfit(bp_stack(i,:)',sig_stack(i,:)');
    x= -10:1:max(max(bp_stack,[],2));
    y = p(1)+p(2)*x;
    figure();
    plot(x,y);
    hold on; 
    scatter(bp_stack(i,:)',sig_stack(i,:)')
    xlabel("ABP (mmHg)")
    ylabel("CBFi (cm/s^2)*10^9");
    title("Regression Plot(r_s_d")
    ylim([0 max((sig_stack(i,:)))+1])
   
%     switch i
%         case 1
%             title("ABP vs TCD regression")
%             ylabel("CBFV(cm/s)")
%         case 2
%             title("ABP vs CBFi(r_s_d=1cm) regresson")
%         case 3
%             title("ABP vs CBFi(r_s_d=1.5cm) regresson")
%         case 4
%             title("ABP vs CBFi(r_s_d=2cm) regresson")
%         case 5
%             title("ABP vs CBFi(r_s_d=2.5cm) regresson")
%         otherwise
%             title("Regression graph")
%     end
    hold off;
    ccp(i) = -p(1)/p(2);   % simplifying the linear equation y=m*x + c for y=0 will lead to this equation 
    fprintf("%d",ccp(i));
end

end
