%% Loading the file
% filename=strcat('D:\Jignesh\MSc Western Uni\Research MSc\Data\Matlab export_carotid occlusion_08062021_TCD DCS.mat');
% load(filename)

%% Upscaling the data by 3
data_u = interp(data,3);
tcd = data_u(1:14400);
dcs_1 = data_u(14501:21500);
dcs_3 = data_u(28801:43200);

%% Plotting the signal
% figure();
% x_fig = (1:1:600)/20
% plot(x_fig,dcs_3a(1:600))
% xlabel('Time (s)')
% ylabel('aDb *10^9 ')
% title('Original signal: DCS-3cm (Fs= 20Hz)')

%% Plotting the frequency spectrum
Fs = 20;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
signal = dcs_3a;
L = length(signal);             % Length of signal
t = (0:L-1)*T;  

Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
%% Finding the minima to find the starting of the signal

minima = islocalmin(dcs_1,'MinProminence',10);
x = 1:length(minima);
plot(x,dcs_1,x(minima),dcs_1(minima),'r*');

%% Plotting the data based on the minima
ini = dcs_1(1:50);
count = 1;
avg = ini;
for i=1:1:length(minima)
    if (minima(i)==1) && (i+49<=length(minima))
        fprintf("%d\n",i)
        plot(dcs_1(i:i+49))
        hold on
        avg = avg+dcs_1(i:i+49);
        count = count+1;
    end
end
hold off
avg = avg/count;
figure()
plot(avg)

%% Plotting the TCD signal

minima = islocalmin(tcd,'MinProminence',10);
x = 1:length(minima);
plot(x,tcd,x(minima),tcd(minima),'r*');

l = length(tcd(1:5000));
ini = tcd(5:54);

cyc =zeros(sum(minima==1)-1,50); 
cyc(1,:)= ini;
count = 1;
avg = ini;
for i=51:1:length(minima)
    if (minima(i)==1) && (i+49<=length(minima))
        count = count+1;
        plot(tcd(i:i+49))
        hold on
        avg = avg+tcd(i:i+49);
        cyc(count,:) = tcd(i:i+49);
    end
end
hold off
avg = avg/count;
x = (1:1:50)/60;
figure()
plot(avg)

%Plotting the ensemble average
ensavg = mean(cyc,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cyc,[],1)/sqrt(count);                             % Calculate 95% Confidence Intervals         
figure()
plot(x, ensavg, '-r', 'LineWidth',1)
hold on
plot(x, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
grid
legend('Ensemble Average', '95% Confidence Intervals')
xlabel('Time (s)')
ylabel('BFi')
title("Ensemble average of TCD upsampled signal")


%% Finding the signal using the find signal
% take the initial segment as a prototype signal to match it for the whole
% signal and then try to find the similar occurances in a singal. Probably
% this will not work cause the signal shape may not be completely identical
% to the prototype

% findsignal

%% Taking the ensemble average

l = length(data_u);
ini = data_u(5:54);
avg = ini;
for i=55:50:l-10
    avg = (avg+data_u(i:i+49));
%     hold on;
    plot(avg);
end
avg = avg/((l-10)/50);
plot(avg);

%% Processing using the actual data
dcs_1a = data(4850:7110);
%% Finding the minima to find the starting of the signal

minima = islocalmin(dcs_1a,'MinProminence',10);
x = 1:length(minima);
figure();
plot(x,dcs_1a,x(minima),dcs_1a(minima),'r*');

%% Plotting the data based on the minima
ini = dcs_1a(1:16);
count = 1;
avg = ini;
for i=1:1:length(minima)
    if (minima(i)==1) && (i+15<=length(minima))
        count = count+1;
        plot(dcs_1a(i:i+15))
        hold on
        avg = avg+dcs_1a(i:i+15);
        
    end
end
hold off
avg = avg/count;
figure()
plot(avg)

%% Finding the minima to find the starting of the signal
sg_lp_30 = sg_lp_03(1:604);
sg_lp_30(109:112) = [];
% sg_lp_30(393:441) = [];
% sg_lp_30(439:443) = [];
minima = islocalmin(sg_lp_30,'MinSeparation',14,'MinProminence',5);
x = 1:length(minima);
figure();
plot(x,sg_lp_30,x(minima),sg_lp_30(minima),'r*');

%% Plotting the data based on the minima
ini = sg_lp_30(1:16);
cyc =zeros(sum(minima==1)-1,16); 
cyc(1,:)= ini;
count = 1;
avg = ini;
for i=17:1:length(minima)
    if (minima(i)==1) && (i+15<=length(minima))
        count = count+1;
        plot(sg_lp_30(i:i+15))
        hold on
        avg = avg+sg_lp_30(i:i+15);
         cyc(count,:) = sg_lp_30(i:i+15);
        
    end
end
hold off
avg = avg/count;
x = (1:1:16)/20;
figure()
plot(avg)

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
%% Finding the minima to find the starting of the upsampled signal 

sg_lp_30_u = interp(sg_lp_30,3);
minima_u = islocalmin(sg_lp_30_u,'MinSeparation',42,'MinProminence',5);
x_u = 1:length(minima_u);
figure();
plot(x_u,sg_lp_30_u,x_u(minima_u),sg_lp_30_u(minima_u),'r*');

%% Plotting the data based on the minima for the upsampled signal
ini = sg_lp_30_u(1:48);
cycu =zeros(sum(minima_u==1)-1,48); 
cycu(1,:)= ini;
count_u = 1;
avg = ini;
for i=49:1:length(minima_u)
    if (minima_u(i)==1) && (i+47<=length(minima_u))
        count_u = count_u+1;
        fprintf("%d\n",i)
        plot(sg_lp_30_u(i:i+47))
        hold on
        avg = avg+sg_lp_30_u(i:i+47);
        cycu(count_u,:) = sg_lp_30_u(i:i+47);
    end
end
hold off
x_u = (1:1:48)/60;
avg = avg/count_u;
figure()
plot(x_u,avg)

%Plotting  the ensemble average
ensavg = mean(cycu,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cycu,[],1)/sqrt(count_u);                             % Calculate 95% Confidence Intervals
figure()
plot(x_u, ensavg, '-r', 'LineWidth',1)
hold on
plot(x_u, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x_u, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
grid
legend('Ensemble Average', '95% Confidence Intervals')
xlabel('Time (s)')
ylabel('aDb *10^9')
title("DCS-1.5cm Upsampled signal")

%% Subplotting the original singal and upscaled signal

subplot(1,2,1)
ensavg = mean(cyc,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cyc,[],1)/sqrt(count);                             % Calculate 95% Confidence Intervals                                    
plot(x, ensavg, '-r', 'LineWidth',1)
hold on
plot(x, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
xlabel('Time (s)')
ylabel('aDb *10^9 ')
title("Ensemble average of original signal")
grid
legend('Ensemble Average', '95% Confidence Intervals')

subplot(1,2,2)
ensavg = mean(cycu,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cycu,[],1)/sqrt(count_u);                             % Calculate 95% Confidence Intervals  
plot(x_u, ensavg, '-r', 'LineWidth',1)
hold on
plot(x_u, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x_u, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
xlabel('Time (s)')
ylabel('aDb *10^9')
title("Ensemble average of Upsampled signal")
grid
legend('Ensemble Average', '95% Confidence Intervals')


%% Similar processing for the DCS 3CM data
dcs_3a = data(9650:12000);
sg_lp_30 = dcs3a_sg_lp(1:620);
% sg_lp_30 = interp(sg_lp_30, 3)
% sg_lp_30(393:441) = [];
% sg_lp_30(439:443) = [];
minima3 = islocalmin(sg_lp_30,'MinSeparation',14,'MinProminence',3);
x3 = 1:length(minima3);
figure();
plot(x3,sg_lp_30,x3(minima3),sg_lp_30(minima3),'r*');


ini = sg_lp_30(1:17);
cyc =zeros(sum(minima3==1)-1,17); 
cyc(1,:)= ini;
count = 2;
avg = ini;
for i=18:1:length(minima3)
    if (minima3(i)==1) && (i+16<=length(minima3))
        fprintf("%d\n",i)
        plot(sg_lp_30(i:i+16))
        hold on
        avg = avg+sg_lp_30(i:i+16);
        cyc(count,:) = sg_lp_30(i:i+16);
        count = count+1;
    end
end
hold off
avg = avg/count;
x3 = (1:1:17)/20;
figure()
plot(avg)

%Plotting the ensemble average
ensavg = mean(cyc,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cyc,[],1)/sqrt(count);                             % Calculate 95% Confidence Intervals         
figure()
subplot(1,2,1)
plot(x3, ensavg, '-r', 'LineWidth',1)
hold on
plot(x3, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x3, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
grid
xlabel('Time (s)')
ylabel('aDb *10^9 ')
title('Original Signal')
legend('Ensemble Average', '95% Confidence Intervals')


sg_lp_30_u = interp(sg_lp_30,3);

minima3u = islocalmin(sg_lp_30_u,'MinSeparation',14,'MinProminence',3);
% x3u = 1:length(minima3u);
% figure();
% plot(x3u,sg_lp_30_u,x3u(minima3u),sg_lp_30_u(minima3),'r*');

ini = sg_lp_30_u(1:51);
cycu =zeros(sum(minima3u==1)-1,51); 
cycu(1,:)= ini;
count = 2;
avg = ini;
for i=52:1:length(minima3u)
    if (minima3u(i)==1) && (i+51<=length(minima3u))
        fprintf("%d\n",i)
%         plot(sg_lp_30_u(i:i+50))
%         hold on
        avg = avg+sg_lp_30_u(i:i+50);
        cycu(count,:) = sg_lp_30_u(i:i+50);
        count = count+1;
    end
end
hold off
avg = avg/count;
x3u = (1:1:51)/60;

%Plotting the ensemble average
ensavg = mean(cycu,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(cycu,[],1)/sqrt(count);                             % Calculate 95% Confidence Intervals         
subplot(1,2,2)
plot(x3u, ensavg, '-r', 'LineWidth',1)
hold on
plot(x3u, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x3u, ensavg-ci95, ':g', 'LineWidth',1.5)
hold off
grid
title('DCS - 2.7cm Upsampled Signal')
xlabel('Time (s)')
ylabel('aDb *10^9 ')
legend('Ensemble Average', '95% Confidence Intervals')