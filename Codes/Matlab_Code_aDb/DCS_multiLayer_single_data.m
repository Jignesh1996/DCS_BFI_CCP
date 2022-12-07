fclose('all'); % Close all open files
clear all; 
close all;
clc;

%%

% load('data_tau.csv')
% load('BFi_changes.mat')

data_cut= [510 570; 30 90; 30 90; 30 90; 270 330; 270 330; 30 90; 510 570; 30 90];
data_cut2=[270 330; 90 150; 90 150; 90 150; 90 150; 90 150; 90 150; 90 150; 90 150];

for subject=2:9

filename=strcat(pwd,'\2 layer model\00',num2str(subject),'_TCD_g2.mat')    
load(filename);

g2_1_temp=g2(:,data_cut(subject,1)/0.05+1:data_cut(subject,2)/0.05,:);

n = 60; % average every n values
for chan=1:4
    for tau=1:50
        a=squeeze(g2_1_temp(chan,:,tau));
        g2_1(chan,tau,:) = arrayfun(@(i) mean(a(i:i+n-1)),1:n:length(a)-n+1)'; % the averaged vecto
    end
end
clear g2

%

rsd=[1 1.5 2 2.5];

% mua = 0.1; %cm^-1 baseline absorption coefficient
% mus = 10; %cm^-1 baseline reduced scattering coefficient

mua = OpticalProp(subject,1); %cm^-1 baseline absorption coefficient
mus = OpticalProp(subject,2); %cm^-1 baseline reduced scattering coefficient

temp_res=10;

baseline_t=10; % time in seconds
baseline=baseline_t/temp_res;
tau_values=data_tau;

for r=1:4
        for i=1:size(g2_1,3)
            g2_temp(i,:)=squeeze(g2_1(r,:,i));
            LB = [0.1e-9];
            UB = [10e-8];
            rho=rsd(r);
            Starting = 1e-9; %[aDb, Beta; cm^2/s, a.u.]
            beta= g2_1(r,1,i); %0.1568;
            % beta= squeeze(g2_new(1,1)); %0.1568;
            options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
            [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rho,beta);
            aDb1(r,i) = FittedParams(1);
        end
end
            

filename=strcat(pwd,'\2 layer model\00',num2str(subject),'_CUFF_g2.mat')    
load(filename);

g2_1_temp=g2(:,data_cut2(subject,1)/0.05+1:data_cut2(subject,2)/0.05,:);

n = 60; % average every n values
for chan=1:4
    for tau=1:50
        a=squeeze(g2_1_temp(chan,:,tau));
        g2_1(chan,tau,:) = arrayfun(@(i) mean(a(i:i+n-1)),1:n:length(a)-n+1)'; % the averaged vecto
    end
end
clear g2

rsd=[1 1.5 2 2.5];

% mua = 0.1; %cm^-1 baseline absorption coefficient
% mus = 10; %cm^-1 baseline reduced scattering coefficient

mua = OpticalProp(subject,1); %cm^-1 baseline absorption coefficient
mus = OpticalProp(subject,2); %cm^-1 baseline reduced scattering coefficient

temp_res=10;

baseline_t=10; % time in seconds
baseline=baseline_t/temp_res;
tau_values=data_tau;

for r=1:4
        for i=1:size(g2_1,3)
            g2_temp(i,:)=squeeze(g2_1(r,:,i));
            LB = [0.1e-9];
            UB = [10e-8];
            rho=rsd(r);
            Starting = [1e-9]; %[aDb, Beta; cm^2/s, a.u.]
            beta= g2_1(r,1,i); %0.1568;
            % beta= squeeze(g2_new(1,1)); %0.1568;
            options = optimset('Display','final','TolX',1e-30,'MaxIter',2000000, 'MaxFunEvals', 200000);
            [FittedParams] = fminsearchbnd(@Brownian_fitting,Starting,LB,UB,options,tau_values,g2_temp(i,:),mua,mus,rho,beta);
            aDb2(r,i) = FittedParams(1);
        end
end
    
BFi_all(subject,:,:)=aDb1;
BFi_all2(subject,:,:)=aDb2;

end
%% Analyze INDIVIDUAL
% Determine the Blood Flow Index for each Channel
tau=data_tau;
time=(1:1:20)*3;


fig1=figure('units','centimeters', 'Position',[2 2 35 15]) %18 width 15 heigh
for i=1:4
    subaxis(2,4,i,'SpacingVert',0.07,'SpacingHoriz',0.05,'MR',0.05, 'ML',0.05,'MT',0.08,'MB',0.14)
    plot(time, squeeze(BFi_all(:,i,:))');
    set(gca,'xticklabel',{})
    if i==1
        ylabel('aDb_{1cm}')
    end
    
    subaxis(2,4,i+4)
    plot(time, squeeze(BFi_all2(:,i,:))');
    xlabel('Time (s)')
    
    if i==1
        ylabel('aDb_{2cm}')
    end
    
end


%% Analyze INDIVIDUAL
% Determine the Blood Flow Index for each Channel
tau=data_tau;

for subject=1:9
    for chan=1:4
        X=squeeze(BFi_all(subject,chan,:))';
        BFi_all_scaled(subject,chan,:)=X/mean(X(1,6:10));
        X=squeeze(BFi_all2(subject,chan,:))';
        BFi_all_scaled2(subject,chan,:)=X/mean(X(1,6:10));
    end
end



fig1=figure('units','centimeters', 'Position',[2 2 35 15]) %18 width 15 heigh
for i=1:4
    subaxis(2,4,i,'SpacingVert',0.07,'SpacingHoriz',0.05,'MR',0.05, 'ML',0.06,'MT',0.08,'MB',0.14)
    plot(time, squeeze(BFi_all_scaled(:,i,:))');
    set(gca,'xticklabel',{})
    if i==1
        ylabel('aDb_{no cuff}')
    end
    rsd=[1,1.5,2, 2.5]
    tit=strcat('r_{SD}=',num2str(rsd(i)),'cm');
    title(tit)
    
    subaxis(2,4,i+4)
    plot(time, squeeze(BFi_all_scaled2(:,i,:))');
    xlabel('Time (s)')
    
    if i==1
        ylabel('aDb_{cuff}')
    end
    
end

%%
%% Analyze INDIVIDUAL
% Determine the Blood Flow Index for each Channel
tau=data_tau;

for subject=1:9
    for chan=1:4
        X=squeeze(BFi_all(subject,chan,:))';
        BFi_all_scaled(subject,chan,:)=100*(X/mean(X(1,6:10))-1);
        X=squeeze(BFi_all2(subject,chan,:))';
        BFi_all_scaled2(subject,chan,:)=100*(X/mean(X(1,6:10))-1);
    end
end



fig1=figure('units','centimeters', 'Position',[2 2 35 15]) %18 width 15 heigh
for i=1:4
    subaxis(2,4,i,'SpacingVert',0.07,'SpacingHoriz',0.05,'MR',0.05, 'ML',0.08,'MT',0.08,'MB',0.14)
    plot(time, squeeze(BFi_all_scaled(:,i,:))');
    set(gca,'xticklabel',{},'ylim',[-100 100])
    if i==1
        ylabel('aDb_{no cuff}')
    end
        rsd=[1,1.5,2, 2.5]
    tit=strcat('r_{SD}=',num2str(rsd(i)),'cm');
    title(tit)
    
    subaxis(2,4,i+4)
    plot(time, squeeze(BFi_all_scaled2(:,i,:))');
    xlabel('Time (s)')
    set(gca,'ylim',[-100 100])
    if i==1
        ylabel('aDb_{cuff}')
    end
    
end

%% Analyze INDIVIDUAL
% Determine the Blood Flow Index for each Channel
tau=data_tau;

for subject=1:9
    for chan=1:4
        X=squeeze(BFi_all(subject,chan,:))';
        BFi_all_scaled(subject,chan,:)=100*(X/mean(X(1,6:10))-1);
        X=squeeze(BFi_all2(subject,chan,:))';
        BFi_all_scaled2(subject,chan,:)=100*(X/mean(X(1,6:10))-1);
    end
end




fig1=figure('units','centimeters', 'Position',[2 2 40 12]) %18 width 15 heigh
for i=1:4
    
    mean_BFi=squeeze(mean(BFi_all_scaled(:,i,:),1));
    std_BFi=std(squeeze(BFi_all_scaled(:,i,:)));

    mean_BFi2=squeeze(mean(BFi_all_scaled2(:,i,:),1));
    std_BFi2=std(squeeze(BFi_all_scaled2(:,i,:)));
    
    
    subaxis(1,4,i,'SpacingVert',0.07,'SpacingHoriz',0.02,'MR',0.02, 'ML',0.06,'MT',0.085,'MB',0.15)
    rectangle('Position', [30 -61 15  122], 'Facecolor', [0.9 0.9 0.9]);
    lineProps.col{1} = 'k';
    lineProps.edgestyle = ':';
    mseb(time,mean_BFi,std_BFi,lineProps,1)
    lineProps.col{1} = 'r';
    lineProps.edgestyle = ':';
    mseb(time,mean_BFi2,std_BFi2,lineProps,1) 
    
    rsd=[1,1.5,2, 2.5]
    tit=strcat('r_{SD}=',num2str(rsd(i)),'cm');
    title(tit)
    
    set(gca,'ylim',[-60 60],'xlim',[10 60])
    if i==1
        ylabel('\DeltaBFi (%)')
        
    else
        set(gca,'yticklabel',{})
    end
    xlabel('Time (s)')   
   
    if i==4
        legend('Cuff OFF','Cuff ON')
    end
    
end
%%

BFI=squeeze(BFi_all_scaled(:,1,:));
BFI(11:19,:)=squeeze(BFi_all_scaled(:,4,:));

BFI2=squeeze(BFi_all_scaled2(:,1,:));
BFI2(11:19,:)=squeeze(BFi_all_scaled2(:,4,:));

%%
tau=data_tau;

for subject=1:9
    for chan=1:4
        X=squeeze(BFi_all(subject,chan,:))';
        BFi_all_scaled(subject,chan,:)=100*(X/mean(X(1,6:10))-1);
        X=squeeze(BFi_all2(subject,chan,:))';
        BFi_all_scaled2(subject,chan,:)=100*(X/mean(X(1,6:10))-1);
    end
end

fig1=figure('units','centimeters', 'Position',[2 2 45 18]) %18 width 15 heigh

for i=1:4
    for subject=1:9
    subaxis(4,9,subject+9*(i-1),'SpacingVert',0.01,'SpacingHoriz',0.02,'MR',0.02, 'ML',0.06,'MT',0.04,'MB',0.15)
    rectangle('Position', [30 -71 15  142], 'Facecolor', [0.9 0.9 0.9]);
    hold on
    plot(time,squeeze(BFi_all_scaled(subject,i,:)),'k','LineWidth',2)
    plot(time,squeeze(BFi_all_scaled2(subject,i,:)),'r','LineWidth',2)
    rsd=[1,1.5,2, 2.5]

    set(gca,'ylim',[-70 70],'xlim',[20 60])
    if subject==1&&i==1
        ylabel('\DeltaBFi_{1cm} (%)')
    elseif subject==1&&i==2
        ylabel('\DeltaBFi_{1.5cm} (%)')
    elseif subject==1&&i==3
        ylabel('\DeltaBFi_{2cm} (%)')
    elseif subject==1&&i==4
        ylabel('\DeltaBFi_{2.5cm} (%)')
%         set(gca,'yticklabel',{})
    end
    if subject>1
        set(gca,'yticklabel',{})
    end
%         tit=strcat('r_{SD}=',num2str(rsd(i)),'cm');
%     title(tit)
    if i<4
       set(gca,'xticklabel',{})
    else
        xlabel('Time (s)')
    end
    

       legend('OFF','ON')

    
    end
    
end

