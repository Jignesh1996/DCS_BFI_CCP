% close all;
% clear all;


% To insert the file, Either use the import data function and import the
% Data.mat file from the folder or you can define the file_path and then
% use the load function to get the data.


%The following 5 lines of code reassigns the different channels of the raw
%data to new variable and then recombine the later 3 channels into 1.

g2(1,:,:)=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2_1_temp=squeeze(Data(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(Data(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(Data(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(Data(:,4,:)-1); %g2-1 curve generation

for i=1:size(g2,2)
    g2(2,i,:)=( g2_2_temp(i,:)+g2_3_temp(i,:)+g2_4_temp(i,:))/3;
end


rho = [1 2.7]; %source detector separations in cm 
mua = 0.17; %cm^-1 baseline absorption coefficient
mus = 10; %cm^-1 baseline reduced scattering coefficient

tau_values=Data_tau;

%Following code fits the g2 curve to get the BFi values for both short and
%long source detector separation.


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
        aDb1(chan,i) = FittedParams(1);
    end
end

plot(aDb1');

