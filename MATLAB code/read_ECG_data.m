% == Load and read ECG data ==
% This function will load and read the data of the ECG measurement, it
% gives both the measurement data and the start time of each measurement

function [ECG_data, ECG_start_time] = read_ECG_data(date, meas_nr)
    % The .txt file you want to read (contains measurement data and start time)
    ECG_filename = strcat(date, "_ekg", meas_nr, ".txt");
    
    % The base path of the files
    base_path = strcat("C:\Users\tess-\OneDrive\Skrivbord\Kandidatarbete\data\", date, "\");
    
    ECG_data_path = strcat(base_path, ECG_filename);
    ECG_data_raw = transpose(readlines(ECG_data_path));
    
    ECG_data_split = zeros(numel(ECG_data_raw)-4,4);
    
    for i = (4:numel(ECG_data_raw)-1)
        ECG_data_split (i-3,:) = strsplit(ECG_data_raw(i));
    end
    
    ECG_data = ECG_data_split(:,3);
    
    ECG_json_raw = ECG_data_raw(2);
    ECG_json_raw_trunc = extractAfter(ECG_json_raw, "# ");
    
    ECG_json = jsondecode(ECG_json_raw_trunc);
    ECG_time = getfield(getfield(ECG_json, 'x00_07_80_4B_15_F7'), 'time');
    
    ECG_start_time = datetime(ECG_time, "InputFormat", "HH:mm:ss.SSS");  
end