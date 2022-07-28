function avg = ens_avg(signals,ttl)

% The issue with the ensemble averaging with NaNs is that the extra signal
% at the end iduces discontinuities and the signal goes down till the
% maximum time of the cycle. To compensate for that, I'm cutting the signal
% till the minimum point. 

x = (1:1:length(signals(1,:)))/1000;
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
title(ttl)

end