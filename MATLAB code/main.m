% == Specificy the meas to be analysed ==
date = "2023_03_08";    % date of meas
meas_nr = "2";          % meas number

% == Load radar parameters ==
radar_parameters_B;     % write name of file with parameters used in meas

% This will give access to the following variables, write over them if
% needed but remember that they are derived from each other

num_frames;             % number of frames = number of chirps
num_samples;            % amount of ADC samples per chirp               
T_f;                    % period of frames = period of chirps
T_s;                    % period of samples per chirp

% == Read meas data and start time of meas ==
% Radar
[radar_data, radar_start_time] = read_radar_data(date, meas_nr, num_samples, num_frames);

% ECG
[ECG_data, ECG_start_time] = read_ECG_data(date, meas_nr);


% == Perform Pan Tompkins on ECG data to detect R-peaks ==
%unique(isfinite(ECG_data)) % Use to check if ECG_data contains any NaNs
[qrs_amp, qrs_i_raw, delay] = pan_tompkin(ECG_data, 1000,0);
%qrs_i = qrs_i_raw - delay; % Shift due to delay in Pan Tompkin
qrs_i = qrs_i_raw; % Doesn't seem to be a delay

% Cut off start and end of ECG data to match radar meas
t_diff = floor(seconds(radar_start_time - ECG_start_time));
% The time between the start of both meas
t_meas = 20;                            % Length of the radar meas
T_ECG = 1/1000;                         % Period of ECG samples
% Truncate the data
[qrs_i_trunc, qrs_amp_trunc, ECG_data_trunc] = truncate_ECG(t_diff, t_meas, qrs_i, qrs_amp, ECG_data, T_ECG);


% == Perform range FFT on radar data ==
% Performs an FFT on each row (each chirp) in our radar data
range_fft = fft(radar_data,[],2);

fft_magnitude = abs(range_fft); % The magnitudes of each FFT
fft_phase = angle(range_fft);   % The angles of each FFT


% == Find and select range bin(s) == 
% To plot the detected peaks, switch the '0' to '1'
[possible_bins, selected_bin] = find_range_bins(num_frames, num_samples, T_f, fft_magnitude, 0);

% If the unwrapped phase does not look good for the most frequent bin you
% can switch the selected index to one of the other possible bins and see
% if the signal looks better. In most cases the no. of possible bins < 3.


% == Extract and unwrap the phase signal ==
% Extract the phase signal based on the most frequent range bin
% Currently we only use one bin, might want to sum over multiple bins
wrapped_phase = fft_phase(:,selected_bin);  % Wrapped phase signal
unwrapped_phase = unwrap(wrapped_phase);    % The unwrapped phase signal


% == Define filters for the unwrapped phase signal ==
fs_radar = 1/T_f;           % Sampling frequency of chirps

% Highpass filter to filter out breathing below 1 Hz
f_hp = [0.3 1];             % Band edges of filter
amp_hp = [0 1];             % Band amplitude, highpass 
dev_hp = [0.0005 0.01];     % Maximum deviation

[n_hp,Wn_hp,beta_hp,ftype_hp] = kaiserord(f_hp,amp_hp,dev_hp,fs_radar);
hp_filter = fir1(n_hp,Wn_hp,ftype_hp,kaiser(n_hp+1,beta_hp),'noscale');

% Lowpass filter to filter out high frequencies above 10 Hz
f_lp = [10 15];             % Band edges of filter
amp_lp = [1 0];             % Band amplitude, highpass 
dev_lp = [0.05 0.01];       % Maximum deviation

[n_lp,Wn_lp,beta_lp,ftype_lp] = kaiserord(f_lp,amp_lp,dev_lp,fs_radar);
lp_filter = fir1(n_lp,Wn_lp,ftype_lp,kaiser(n_lp+1,beta_lp),'noscale');


% == Filter the unwrapped phase signal ==
% Apply highpass filter
hp_filtered_phase = filtfilt(hp_filter, 1, unwrapped_phase);

% Apply the lowpass filter
filtered_phase = filtfilt(lp_filter, 1, hp_filtered_phase);

% Smooth the filtered signal to get rid of noise by using moving average
smooth_phase = smoothdata(filtered_phase, 'movmean', 5);


% == Locate the peaks in the phase signal ==
[pks, locs] = findpeaks(smooth_phase);

% == Plot data ==
% Run this section first to load in the time vectors
t_ECG = T_ECG*(0:numel(ECG_data_trunc)-1);  % Time vector for trunc ECG data
t_radar = T_f*(0:num_frames-1);             % Time vector for chirps/frames


%% Plot ECG data
hold on
plot(t_ECG, rescale(ECG_data_trunc))        % Plot truncated ECG data
plot(t_ECG(qrs_i_trunc), qrs_amp_trunc, "o")% Plot detected R-peaks
title('ECG measurement and located R-peaks')
xlabel('Time [s]');
ylabel('Voltage [arb. unit]');
hold off

%% Plot unwrapped phase signal
plot(t_radar, unwrapped_phase)
title('Unwrapped phase signal')
xlabel('Time [s]');
ylabel('Amplitude [arb. unit]');

%% Plot signal before and after filter
hold on
plot(t_radar, unwrapped_phase)      % Unfiltered phase
plot(t_radar, hp_filtered_phase)    % HP filtered phase
plot(t_radar, filtered_phase)       % HP and LP filtered phase
plot(t_radar, smooth_phase)         % Filtered and smoothed phase
legend('Unfiltered', 'HP filtered', 'HP and LP filtered', 'Filtered and smoothed')
title('Phase signal')
xlabel('Time [s]');
ylabel('Amplitude [arb. unit]');
hold off

%% Plot filtered and smoothed signal with detected peaks
hold on
plot(t_radar, smooth_phase)     % Filtered and smoothed phase
plot(t_radar(locs), pks, "o")   % Located peaks
title('Phase signal with detected peaks')
xlabel('Time [s]');
ylabel('Amplitude [arb. unit]');
hold off

%% Plot signal and detected peaks compared to R-peaks
hold on
plot(t_radar, smooth_phase)     % Filtered and smoothed phase
plot(t_radar(locs), pks, "o")   % Located peaks
plot(t_ECG(qrs_i_trunc), 0, "ok") % R-peaks
title('Phase signal with detected peaks')
xlabel('Time [s]');
ylabel('Amplitude [arb. unit]');
hold off

%% Plot inverted signal and detected peaks compared to R-peaks
inv_phase = -1*smooth_phase;
[inv_pks, inv_locs] = findpeaks(inv_phase);
hold on
plot(t_radar, inv_phase)     % Filtered and smoothed phase
plot(t_radar(inv_locs), inv_pks, "or")   % Located peaks
plot(t_ECG(qrs_i_trunc), 0, "ok") % R-peaks
title('Inverted phase signal with detected peaks')
xlabel('Time [s]');
ylabel('Amplitude [arb. unit]');
hold off
