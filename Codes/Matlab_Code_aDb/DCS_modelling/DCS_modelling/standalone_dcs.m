function aDb = standalone_dcs(Data_avg,Data_tau)
%Calculates the aDb values for each channels
tau_range = (1:50);
if size(Data_avg,2)==4
    g2(1,:,tau_range)=squeeze(Data_avg(:,1,tau_range)-1); %g2-1 curve generation
    g2_2_temp=squeeze(Data_avg(:,2,tau_range)-1); %g2-1 curve generation
    g2_3_temp=squeeze(Data_avg(:,3,tau_range)-1); %g2-1 curve generation
    g2_4_temp=squeeze(Data_avg(:,4,tau_range)-1); %g2-1 curve generation
    
    % average g2 curve for large source detector separation
    for i=1:size(g2,2)
        g2(2,i,:)=(g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
    end
else
    g2(1,:,:) = Data_avg(:,1,tau_range);
    g2(2,:,:) = Data_avg(:,2,tau_range);
end

% aDb calculation
rho = [1 2.5]; %source detector separations in cm 
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 8; %cm^-1 baseline reduced scattering coefficient

tau_values=Data_tau(tau_range);
lc = size(g2,1);
ld = size(g2,2);
disp(size(g2))
parfor chan=1:lc
    g2_temp = zeros(size(g2,2),length(tau_range));
    for i=1:ld
        rsd=rho(chan);
        g2_temp(i,:)=squeeze(g2(chan,i,:));
        LB = [0];
        UB = [inf];
        Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
%         beta= squeeze(g2(chan,i,1)); %0.1568;
        beta= squeeze(mean(g2(chan,1:500,1))); %0.1568;
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rsd,beta);
        aDb(chan,i) = FittedParams(1);
    end
end

 %Plotting the fit
Channel=1;
Curve_no=1;
rho = [1 2.5];

beta=g2(Channel, Curve_no,1);
aDb1=aDb(Channel,Curve_no);

g2_fit=gen_DCS_fit(tau_values,mua,mus,rho(Channel),beta,aDb1);

semilogx(tau_values,squeeze(g2(1,1,:)),'k')
hold on
semilogx(tau_values,g2_fit,'r')


end

