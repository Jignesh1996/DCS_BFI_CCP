%% Load the data
time = 0.1:0.1:150;
DCS_data = DCSandBPcombined{:,2:8};
fing_p = DCS_data(:,1);
dcs_1cm = DCS_data(:,6);
dcs_3cm = DCS_data(:,7);
sys = DCS_data(:,2);
MAP = DCS_data(:,3);
dias = DCS_data(:,4);
hr = DCS_data(:,5);

%% Plotting 


