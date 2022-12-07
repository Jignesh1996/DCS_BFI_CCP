function aDb = sfit_stdaln_dcs(Data_avg,Data_tau)
%Calculates the aDb values for each channels

    g2(1,:)=(Data_avg(1,:)-1); %g2-1 curve generation
    g2_2_temp=(Data_avg(2,:)-1); %g2-1 curve generation
    g2_3_temp=(Data_avg(3,:)-1); %g2-1 curve generation
    g2_4_temp=(Data_avg(4,:)-1); %g2-1 curve generation
    
    
    % average g2 curve for large source detector separation
%     for i=1:size(g2,2)
%         g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
%     end
    g2(2,:) = (g2_2_temp+g2_3_temp+g2_4_temp)/3;


% aDb calculation
rho = [1 2.5]; %source detector separations in cm 
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10; %cm^-1 baseline reduced scattering coefficient

tau_values=Data_tau;

for chan=1:size(g2,1)
        rsd=rho(chan);
        g2_temp(chan,:)=(g2(chan,:));
        LB = [0];
        UB = [inf];
        Starting = [1e-10]; %[aDb, Beta; cm^2/s, a.u.]
        beta= squeeze(g2(chan,1)); %0.1568;
%         beta= squeeze(mean(g2(chan,1:500,1))); %0.1568;
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(chan,:),mua,mus,rsd,beta);
        aDb(chan) = FittedParams(1);
end

end

