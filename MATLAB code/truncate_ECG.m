% == Truncate ECG data to start at the same time as radar meas ==
% Since the ECG meas is started before the radar meas the data is not
% matched in time, this function will give a new start and end index for 
% the ECG data which can be used to truncate the beginning and end of the 
% ECG data to sync it with the radar meas
%
% Inputs:
% radar_start_time = start time of radar meas as datetime object
% ECG_start_time = start time of ECG meas as datetime object
% ECG_data = data to be truncated
% t_meas = length of radar measurement
% add_t_diff = additional time difference, a positive value would cut off
% more from the beginning of the ECG data
%
% Outputs:
% ECG_trunc = the truncated ECG data

function ECG_trunc = truncate_ECG(radar_start_time, ECG_start_time, ECG_data, T_ECG, t_meas, add_t_diff)
    % Time difference between data sets
    t_diff = floor(seconds(radar_start_time - ECG_start_time)) + add_t_diff;
    
    start_i = floor(t_diff/T_ECG); 
    % The index of the ECG data which corresponds to the start time of the
    % radar meas. Data points before this index will be truncated.
    
    end_i = start_i + t_meas/T_ECG; 
    % The index of the ECG data which corresponds to the end time of the
    % radar meas. Data points after this index will be truncated.
    
    ECG_trunc = ECG_data;
    % Truncate the data 
    if end_i <= numel(ECG_trunc)
        ECG_trunc = ECG_data(start_i:end_i);
    else
        ECG_trunc = ECG_data(start_i:end);
    end
end