% == Load and read radar data ==
% This function will load and read the data of the radar measurement, it
% gives both the measurement data and the start time of each measurement
% It also reshapes the data so that the data of each chirp is stored in
% separate rows.

function [radar_data, radar_start_time] = read_radar_data(date, meas_nr, num_samples, num_frames)
    % The .bin file you want to read (contains measurement data)
    % The files are named on the form "2023_03_02_radar8.bin"
    radarDataFilename = strcat(date, "_radar", meas_nr, ".bin");
    
    % The .csv log file for the same measurement (contains start time of meas)
    % The files are named on the form "2023_03_02_radar8_Raw_LogFile.csv"
    radarLogFilename = strcat(date, "_radar", meas_nr, "_Raw_LogFile.csv");
    
    % The base path of the files
    radarBasePath = strcat("C:\Users\tess-\OneDrive\Skrivbord\Kandidatarbete\data\", date, "\");

    % Read start time of measurement
    radarLogPath = strcat(radarBasePath, radarLogFilename);
    % Read log file
    log_data = readlines(radarLogPath);
    startTimeSplit = split(log_data(17), "- ");
    % Note, may want to use "H" instead of "HH" or "s" instead of "ss".
    radar_start_time = datetime(startTimeSplit(2), "InputFormat", "eee MMM dd HH:mm:ss y");

    % Read radar data
    radarDataPath = strcat(radarBasePath, radarDataFilename);
    radarDataRaw = readDCA1000(radarDataPath, num_samples);
    radarDataOneRow = radarDataRaw(1:2:end); % Here the data is all stored in one row

    % Reshape the data so each chirp is stored in seperate rows
    radar_data = transpose(reshape(radarDataOneRow,[num_samples, num_frames]));
end