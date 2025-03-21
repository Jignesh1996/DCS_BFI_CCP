function aDb = standalone_tr_dcs(Data_avg,Data_tau)
%Calculates the aDb values for each channels
if size(Data_avg,2)==8
    g2(1,:,:)=squeeze(Data_avg(:,1,:)-1); %g2-1 curve generation
    g2_2_temp=squeeze(Data_avg(:,2,:)-1); %g2-1 curve generation
    g2_3_temp=squeeze(Data_avg(:,3,:)-1); %g2-1 curve generation
    g2_4_temp=squeeze(Data_avg(:,4,:)-1); %g2-1 curve generation
    g2_5_temp=squeeze(Data_avg(:,5,:)-1); %g2-1 curve generation
    g2_6_temp=squeeze(Data_avg(:,6,:)-1); %g2-1 curve generation
    g2_7_temp=squeeze(Data_avg(:,7,:)-1); %g2-1 curve generation
    g2_8_temp=squeeze(Data_avg(:,8,:)-1); %g2-1 curve generation
    
    % average g2 curve for large source detector separation
    for i=1:size(g2,2)
        g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
        g2(3,i,:) = (g2_5_temp(i,:)+g2_6_temp(i,:)+g2_7_temp(i,:)+g2_8_temp(i,:))/4;
    end
else
    g2(1,:,:) = Data_avg(:,1,:);
    g2(2,:,:) = Data_avg(:,2,:);
    g2(3,:,:) = Data_avg(:,3,:);
end

% aDb calculation
rho = [1 2.5 2.5]; %source detector separations in cm 
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
        beta= squeeze(mean(g2(chan,1:500,1))); %0.1568;
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rsd,beta);
        aDb(chan,i) = FittedParams(1);
    end
end

 %Plotting the fit
Channel=1;
Curve_no=1;
rho = [1 2.5 2.5];

beta=g2(Channel, Curve_no,1);
aDb1=aDb(Channel,Curve_no);

g2_fit=gen_DCS_fit(Data_tau,mua,mus,rho(Channel),beta,aDb1)

semilogx(Data_tau,squeeze(g2(1,1,:)),'k')
hold on
semilogx(Data_tau,g2_fit,'r')


end

