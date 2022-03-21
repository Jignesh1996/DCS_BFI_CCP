function aDb = hybrid_dcs(Data,Data_tau)
%Calculates the aDb values for each channels
g2(1,:,:)=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2(2,:,:)=squeeze(Data(:,2,:)-1); %g2-1 curve generation
g2(3,:,:)=squeeze(Data(:,3,:)-1); %g2-1 curve generation
g2(4,:,:)=squeeze(Data(:,4,:)-1); %g2-1 curve generation

% aDb calculations

rho = [1 1.5 2 2.5]; %source detector separations in cm 
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10; %cm^-1 baseline reduced scattering coefficient

tau_values=Data_tau;

for chan=1:size(g2,1)
    for i=1:size(g2,2)
        rsd=rho(chan);
        g2_temp(i,:)=squeeze(g2(chan,i,:));
        LB = [0];
        UB = [inf]; 
        Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
        beta= squeeze(g2(chan,i,1)); %0.1568;
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rsd,beta);
        aDb(chan,i) = FittedParams(1);
 
    end
%     semilogx(Data_tau,squeeze(aDb(chan,1),Data_tau,g2_temp(1,:)),Color='b',LineWidth=0.9);
%     it is wrong to plot that, cause we need 50 tau points for adb  as
%     well
end
end
