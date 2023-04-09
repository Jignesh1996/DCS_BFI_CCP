function G1= g1_curve(Data)
Data_avg = Data;
tau_range = (1:50);
if size(Data_avg,2)==4
    g2(1,:,1:length(tau_range))=squeeze(Data_avg(:,1,tau_range)-1); %g2-1 curve generation
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

G1 = abs(sqrt((g2)./g2(:,:,1)));

end