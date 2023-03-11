% == Truncate ECG data to start at the same time as radar meas ==
function [index_shift, offset] = truncate_ECG(qrs_i_shifted, ECG_start_time, radar_start_time, T_ECG)
    meas_delay = radar_start_time - ECG_start_time;
    time_diff = seconds(meas_delay);

    index_start = 1 + time_diff/T_ECG; 

    counter = 0;
    offset = 0;
    for i = (1:numel(qrs_i_shifted))
        if qrs_i_shifted(i) < index_start
        counter = counter +1;
        offset = qrs_i_shifted(i);
        end
    end
    
    index_shift = counter+1;

end