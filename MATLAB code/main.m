% == Load radar parameters ==
radar_parameters_A;     % write name of file with parameters used in meas

% This will give access to the following variables, write over them if
% needed but remember that they are derived from each other

num_frames;             % number of frames = number of chirps
num_samples;            % amount of ADC samples per chirp               
T_f;                    % period of frames = period of chirps
T_s;                    % period of chirp samples


% == Load and read radar data and start time of meas ==
date = "2023_03_02";    % date of meas
meas_nr = "2";          % meas number

[radar_data, radar_start_time] = read_radar_data(date, meas_nr, num_samples, num_frames);


% == Perform range FFT on radar data ==
% Performs an FFT on each row (on each chirp) in our radar data
range_fft = fft(radar_data,[],2);

fft_magnitude = abs(range_fft); % The magnitudes of each FFT
fft_phase = angle(range_fft);   % The angles of each FFT


% == Find and select range bin(s) == 
% Find the range bins of each chirp
[bin_magnitudes, bin_indices] = find_range_bins(num_frames, fft_magnitude);
possible_bins = unique(bin_indices); % Gives a list of possible bins
bin_index = mode(bin_indices); % Gives the most frequent range bin

% If the unwrapped phase does not look good for the most frequent bin you
% can switch the selected index to one of the other possible bins and see
% if the signal looks better. In most cases it is just 1-3 different bins.


% == Extract and unwrap the phase signal ==
% Extract the phase signal based on the most frequent range bin
% Currently we only use one bin
wrapped_phase = fft_phase(:,bin_index);     % Wrapped phase signal
unwrapped_phase = unwrap(wrapped_phase);    % The unwrapped phase signal


% == Plot the unwrapped phase signal ==
t = T_f*(0:num_frames-1);   % Time vector for chirps/frames

plot(t, unwrapped_phase)
title('Unwrapped phase signal')
xlabel('Time [s]');
ylabel('Amplitude');

%%
% == Plot the range FFTs and detected peaks corresponding to range bins ==
plot_FFT_and_bins(num_frames, T_f, num_samples, fft_magnitude, bin_indices)
