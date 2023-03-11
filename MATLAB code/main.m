% == Load radar parameters ==
radar_parameters_B;     % write name of file with parameters used in meas

% This will give access to the following variables, write over them if
% needed but remember that they are derived from each other

num_frames;             % number of frames = number of chirps
num_samples;            % amount of ADC samples per chirp               
T_f;                    % period of frames = period of chirps
T_s;                    % period of chirp samples

% == Load and read data and start time of meas ==
date = "2023_03_08";    % date of meas
meas_nr = "4";          % meas number

[radar_data, radar_start_time] = read_radar_data(date, meas_nr, num_samples, num_frames);

[ECG_data, ECG_start_time] = read_ECG_data(date, meas_nr);

% == Perform range FFT on radar data ==
% Performs an FFT on each row (on each chirp) in our radar data
range_fft = fft(radar_data,[],2);

fft_magnitude = abs(range_fft); % The magnitudes of each FFT
fft_phase = angle(range_fft);   % The angles of each FFT

% == Find and select range bin(s) == 
% Find the range bins of each chirp
[bin_magnitudes, bin_indices] = find_range_bins(num_frames, fft_magnitude);
bin_index = mode(bin_indices); % Gives the most frequent range bin
possible_bins = unique(bin_indices); % Gives a list of possible bins

% If the unwrapped phase does not look good for the most frequent bin you
% can switch the selected index to one of the other possible bins and see
% if the signal looks better. In most cases it is just 1-3 different bins.

% == Extract and unwrap the phase signal ==
% Extract the phase signal based on the most frequent range bin
% Currently we only use one bin
wrapped_phase = fft_phase(:,11);     % Wrapped phase signal
unwrapped_phase = unwrap(wrapped_phase);    % The unwrapped phase signal

% == Plot the unwrapped phase signal ==
t_radar = T_f*(0:num_frames-1);   % Time vector for chirps/frames

plot(t_radar, unwrapped_phase)
title('Unwrapped phase signal')
xlabel('Time [s]');
ylabel('Amplitude');

% == Perform Pan Tompkins on ECG data and plot R-peaks ==
%unique(isfinite(ECG_data)) % Use to check if ECG_data contains any NaNs
[qrs_amp_raw,qrs_i_raw,delay] = pan_tompkin(ECG_data, 1000,0);

qrs_i_shifted = qrs_i_raw - delay; % Shift due to delay in Pan Tompkin

T_ECG = 1/1000;                         % Period of ECG samples
t_ECG = T_ECG*(0:numel(ECG_data)-1);    % Time vector for ECG meas

% Truncate the ECG meas so it starts closerto radar meas
[index_shift, offset] = truncate_ECG(qrs_i_shifted, ECG_start_time, radar_start_time, T_ECG);

qrs_i_trunc = qrs_i_shifted(index_shift:end) - offset;
qrs_amp_trunc = qrs_amp_raw(index_shift:end);
ECG_data_trunc = ECG_data(offset+1:end); % Truncating original data needs to be fixed
t_ECG_trunc = t_ECG(1:numel(ECG_data)-offset);

hold on
plot(t_ECG_trunc, rescale(ECG_data_trunc))      % Plot the original ECG meas
plot(t_ECG(qrs_i_trunc), qrs_amp_trunc, "o")    % Only plot the R-peaks
%plot(t_radar, unwrapped_phase)




%% Commented this out because it is not really needed 
% == Plot the range FFTs and detected peaks corresponding to range bins ==
plot_FFT_and_bins(num_frames, T_f, num_samples, fft_magnitude, bin_indices)
