function [avg,ensemble_curve] = ens_avg(signals,fs,ttl)

% The issue with the ensemble averaging with NaNs is that the extra signal
% at the end iduces discontinuities and the signal goes down till the
% maximum time of the cycle. To compensate for that, I'm cutting the signal
% till the minimum point. 

x = (1:1:length(signals(1,:)))/fs;
count = length(signals(:,1));

%Plotting the ensemble average
ensavg = mean(signals,1,"omitnan");                                                   % Calculate Ensemble Average
ci95 = 1.96*std(signals,[],1,"omitnan")/sqrt(count);                             % Calculate 95% Confidence Intervals         
figure()
plot(x, ensavg, '-r', 'LineWidth',1)
hold on;
plot(x, ensavg+ci95, ':g', 'LineWidth',1.5)
plot(x, ensavg-ci95, ':g', 'LineWidth',1.5)
xlabel("Time (s)");
hold off
grid
legend('Ensemble Average', '95% Confidence Intervals')
avg = ensavg;
ensemble_curve = [avg+ci95;avg;avg-ci95];
title(ttl)

end