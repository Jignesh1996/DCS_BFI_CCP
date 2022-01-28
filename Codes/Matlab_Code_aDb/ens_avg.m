function avg = ens_avg(signals,ttl)
x = (1:1:length(signals(1,:)))/1000;
count = length(signals(:,1));

%Plotting the ensemble average
ensavg = mean(signals,1);                                                   % Calculate Ensemble Average
ci95 = 1.96*std(signals,[],1)/sqrt(count);                             % Calculate 95% Confidence Intervals         
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