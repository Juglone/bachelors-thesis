% == Truncate ECG data to start at the same time as radar meas ==
% Since the ECG meas is started before the radar meas the data is not
% matched in time, this function will give a new start and end index for 
% the ECG data which can be used to truncate the beginning and end of the 
% ECG data

function [qrs_i_trunc, qrs_amp_trunc, ECG_data_trunc] = truncate_ECG(t_diff, t_meas, qrs_i, qrs_amp, ECG_data, T_ECG)
    start_i = 1 + floor(t_diff/T_ECG);
    % The index of the ECG data which corresponds to the start time of the
    % radar meas, if the index of the R-peaks are lower than this they will
    % be truncated
    
    end_i = start_i + t_meas/T_ECG;  
    % The index of the ECG data which corresponds to the end time of the
    % radar meas. Data points after this index will be truncated
    
    start_counter = 1;  % 1 + num elements to be truncated from start
    end_counter = 0;    % Num elements to be truncated from end
    
    for i = (1:numel(qrs_i))
        if qrs_i(i) < start_i
            start_counter = start_counter + 1;
        elseif qrs_i(i) > end_i
            end_counter = end_counter + 1;
        end
    end 
    % Truncate all the data
    qrs_i_trunc = qrs_i(start_counter:end - end_counter) - start_i + 1;
    qrs_amp_trunc = qrs_amp(start_counter:end - end_counter);
    ECG_data_trunc = ECG_data(start_i:end_i);
end