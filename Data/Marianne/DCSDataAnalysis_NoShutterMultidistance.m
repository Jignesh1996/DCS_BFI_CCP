
clear all
close all
%--------------------------------------------------------------------------

% % Folder name  - provide folder name for which you want to convert the data
Folder='20211116-7';

%--------------------------------------------------------------------------

filename=strcat(pwd, '\',Folder,'\','data.mat');
load(filename)

%
temp_res=0.25 %in s

%% marker for time resolved
TR=1;



for i=1:size(Data,1)
        Data_cut(i,:,:)=Data(i,:,:);
        Data_time(1,i)=i*temp_res;
    
end



%%
g2_1_temp=squeeze(Data_cut(:,1,:)-1); %g2-1 curve generation
g2_2_temp=squeeze(Data_cut(:,2,:)-1); %g2-1 curve generation
g2_3_temp=squeeze(Data_cut(:,3,:)-1); %g2-1 curve generation
g2_4_temp=squeeze(Data_cut(:,4,:)-1); %g2-1 curve generation

% average g2 curve for 4 channels
if TR==0
    g2(1,:,:)=(g2_1_temp);
    g2(2,:,:)=(g2_1_temp+g2_2_temp+g2_4_temp)/3;
else
    g2(1,:,:)=(g2_1_temp);
    g2(2,:,:)=(g2_2_temp);
    g2(3,:,:)=(g2_3_temp);
    g2(4,:,:)=(g2_4_temp);
end
 
%   g2=(g2_2_temp)/1;
 
% for stand alone

if TR==0
semilogx(Data_tau,squeeze((g2(1,1,:))),'r')
hold on
semilogx(Data_tau,squeeze((g2(2,1,:))),'b')
ylabel('g_2')
legend('1.5 cm','2.7 cm')

else
% for TR
semilogx(Data_tau,squeeze((g2(1,1,:))),'r')
hold on
semilogx(Data_tau,squeeze((g2(2,1,:))),'g')
semilogx(Data_tau,squeeze((g2(3,1,:))),'b')
semilogx(Data_tau,squeeze((g2(4,1,:))),'magenta')
ylabel('g_2')
legend('1 cm','1.5 cm','2 cm','2.5 cm')
end
 %%
%  j=1;
%  for i=1:moving_window:size(g2,1)
%      g2_new(j,:)=mean(g2(i:i+moving_window-1,:));
%      j=j+1;
%  end
%      

%% aDb calculations
%stand alone

if TR==0
    rho = [1.5 2.7]; %source detector separations in cm 
else 
    rho = [1 1.5 2 2.5]; %source detector separations in cm 
end

mua = 0.1; %cm^-1 baseline absorption coefficient
mus = 10; %cm^-1 baseline reduced scattering coefficient

tau_values=Data_tau;
for r=1:size(rho,2)
    for i=1:size(g2,2)
        g2_temp(i,:)=squeeze(g2(r,i,:))';
        LB = [0.1e-9];
        UB = [10e-8];
        Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
        beta= squeeze((g2(r,i,1:1))); %0.1568;
        options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
        [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rho(1,r),beta);
        aDb1(r,i) = FittedParams(1);
    end
end


% %% Data plotting
% 
% for i=1:size(aDb1,2)
%     aDb1_avg(i)=mean(aDb1);
%     time_avg(i)=mean(Data_time);
% end
% line([60 60],[min(aDb1_avg) max(aDb1_avg)],'LineWidth',2,'Color','red')
% line([75 75],[min(aDb1_avg) max(aDb1_avg)],'LineWidth',2,'Color','green')
% line([195 195],[min(aDb1_avg) max(aDb1_avg)],'LineWidth',2,'Color','green')
% line([210 210],[min(aDb1_avg) max(aDb1_avg)],'LineWidth',2,'Color','red')
% line([330 330],[min(aDb1_avg) max(aDb1_avg)],'LineWidth',2,'Color','green')
% line([345 345],[min(aDb1_avg) max(aDb1_avg)],'LineWidth',2,'Color','red')
% % line([360 360],[min(aDb1_avg) max(aDb1_avg)],'LineWidth',2,'Color','green')
% % line([375 375],[min(aDb1_avg) max(aDb1_avg)],'LineWidth',2,'Color','red')
% % line([390 390],[min(aDb1_avg) max(aDb1_avg)],'LineWidth',2,'Color','green')
% % line([405 405],[min(aDb1_avg) max(aDb1_avg)],'LineWidth',2,'Color','red')
% 
% hold on
% 
% plot(time_avg,aDb1_avg,':k','LineWidth',2)
% 
% title('rCBF')
% % set(gca,'xticklabel',{})
% xlabel('Time (s)')


%%

if TR==0
windowsize=10;
wagi=ones(1,windowsize)/windowsize;
aDb1_plot(1,:)=filtfilt(wagi,1,aDb1(1,:));
aDb1_plot(2,:)=filtfilt(wagi,1,aDb1(2,:));

figure(2)
subplot(2,1,1)
% hold on
% rectangle('Position',[30,-79,30,180],'FaceColor',[0.9 .9 .9],'EdgeColor','none',...
%     'LineWidth',3)
plot(Data_time,aDb1(1,:),'r','LineWidth',1,'MarkerSize',2)
title('rCBF')
legend('1 cm')
hold on
subplot(2,1,2)
plot(Data_time,aDb1(2,:),'k','LineWidth',1,'MarkerSize',2)
legend('2.7cm')
hold on

xlabel('Time (s)')

else
    
windowsize=10;
wagi=ones(1,windowsize)/windowsize;
aDb1_plot(1,:)=filtfilt(wagi,1,aDb1(1,:));
aDb1_plot(4,:)=filtfilt(wagi,1,aDb1(4,:));

figure(2)
subplot(2,1,1)
plot(Data_time,aDb1(1,:),'r','LineWidth',1,'MarkerSize',2)
title('rCBF')
legend('1cm')
hold on
subplot(2,1,2)
plot(Data_time,aDb1(2,:),'k','LineWidth',1,'MarkerSize',2)
legend('2.5 cm')
hold on

xlabel('Time (s)')

end
    

%%

