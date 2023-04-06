% == Specificy the meas to be analysed and load parameters ==
date = "2023_03_14";    % date of meas
meas_nr = "20";         % meas number
both_lanes = 0;         % was data collected from both lanes, 1 = yes, 0 = no

radar_parameters_C;     % name of file with parameters used in meas

% This will give access to the following variables, write over them if
% needed but remember that they are derived from each other

num_frames;     % number of frames = number of chirps
num_samples;    % number of ADC samples per chirp               
T_f;            % period of frames = period of chirps
T_s;            % period of samples per chirp
t_meas;         % length of radar meas in seconds

T_ECG = 1/1000; % period of ECG samples


% == Read meas data and start time of meas ==
% Comment: the functions that read the data are adapted to my computer,
% they need to be adjusted to work for someone elses data paths
% Comment: the start time is no longer used to sync ECG and radar

% Radar
[radar_data, radar_start_time] = read_radar_data(date, meas_nr, num_samples, num_frames, both_lanes);

% ECG
[ECG_data, ECG_start_time] = read_ECG_data(date, meas_nr);


% == Truncate ECG meas to sync with radar meas ==             
% Use the following line if using start times
%ECG_trunc = truncate_ECG(radar_start_time, ECG_start_time, ECG_data, T_ECG, t_meas, 0);

% Use the following lines if both meas were stopped at the same time
num_i_ECG = floor(t_meas/T_ECG);            % num of i corresponding to meas length
ECG_trunc = ECG_data(end-num_i_ECG+1:end);  % truncated ECG data


% == Perform Pan Tompkins on ECG data to detect R-peaks ==
%unique(isfinite(ECG_data)) % use to check if ECG_data contains any NaNs
[qrs_amp, qrs_i_raw, delay] = pan_tompkin(ECG_trunc, 1000,0);
qrs_i = qrs_i_raw - delay;  % shift indexes due to delay in Pan Tompkins
%qrs_i = qrs_i_raw;         % without shift due to delay


% == Perform range FFT on radar data ==
% Performs an FFT on each row (each chirp) in our radar data
range_fft = fft(radar_data,[],2);

fft_magnitude = abs(range_fft); % The magnitudes of each FFT
fft_phase = angle(range_fft);   % The angles of each FFT


% == Find and select range bin(s) == 
% To plot the detected peaks, switch the '0' to '1'
[bin_possible, bin_selected] = find_range_bins(num_frames, num_samples, T_f, fft_magnitude, 0);

% If the unwrapped phase does not look good for the most frequent bin you
% can switch the selected index to one of the other possible bins and see
% if the signal looks better.
% Comment: we want to develop the method of picking the range bin so that 
% we choose the range bin that contains the most information about the 
% heart and less information about random body movement


% == Extract and unwrap the phase signal ==
% Comment: currently we extract the phase signal based on the most frequent 
% range bin, the method of selecting a range bin will be further developed
wrapped_phase = fft_phase(:,bin_selected);  % wrapped phase signal
unwrapped_phase = unwrap(wrapped_phase);    % the unwrapped phase signal 


% == Create filters ==
fs_radar = 1/T_f;           % sampling frequency
fn = 0.5/T_f;               % nyquist frequency

% Highpass filter 8-50 Hz
f_08_50 = [7.6 8];          % band edges of filter
amp_08_50 = [0 1];          % band amplitude
dev_08_50 = [0.01 0.01];    % maximum deviation
[n_08_50,Wn_08_50,beta_08_50,ftype_test] = kaiserord(f_08_50,amp_08_50,dev_08_50,fs_radar);
filter_08_50 = fir1(n_08_50,Wn_08_50,ftype_test,kaiser(n_08_50+1,beta_08_50),'noscale');

% Additional filter to play around with
f_test = [9.6 10];         % band edges of filter
amp_test = [0 1 ];         % band amplitude, highpass 
dev_test = [0.01 0.01];    % maximum deviation
[n_test,Wn_test,beta_test,ftype_test] = kaiserord(f_test,amp_test,dev_test,fs_radar);
filter_test = fir1(n_test,Wn_test,ftype_test,kaiser(n_test+1,beta_test),'noscale');


% == Filter the unwrapped phase signal ==
% Apply highpass filter
filt_phase_08_50 = filtfilt(filter_08_50, 1, unwrapped_phase);
filt_phase_test = filtfilt(filter_test, 1, unwrapped_phase);

% Smooth the filtered signal by using moving average
smooth_phase_08_50 = smoothdata(filt_phase_08_50, 'movmean', 5);
smooth_phase_test = smoothdata(filt_phase_test, 'movmean', 5);


% == Plot data ==
% Run this to load in time vectors and then run a section below to plot
t_ECG = T_ECG*(0:numel(ECG_data)-1);    % time vector for trunc ECG data
t_radar = T_f*(0:num_frames-1);         % time vector for chirps/frames

%% Wavelet analysis to filter out baseline drift and random body movement
% TO BE DEVELOPED

%% Plot the possible range bins or other bins
% bin_possible = []     % add other bins to plot if needed
for i = bin_possible
    hold on
    plot(t_radar, unwrap(fft_phase(:,i)))
end

%% Plot unwrapped phase signal compared to filtered signal
% Unwrapped phase signal
t = tiledlayout(2,1);
ax1 = nexttile;
plot(t_radar, unwrapped_phase)      % unfiltered phase
title('Unwrapped phase signal, BB*, meas 031420')
xlabel('Time [s]');
ylabel('Phase [rad]');

% Filtered signal
ax2 = nexttile;
hold on
%plot(t_radar, filt_phase)          % filtered phase
plot(t_radar, filt_phase_test)      % filtered and smoothed phase

% plot R-peaks
r_peaks = t_ECG(abs(qrs_i));
for i = r_peaks
    xline(i,'--k')
end

legend('Phase signal', 'R-peaks')
title('Filtered phase signal (10-50Hz), BB*, meas 031420')
xlabel('Time [s]');
ylabel('Phase [rad]');
hold off

linkaxes([ax1 ax2],'x')

%% Plot and compare filters
% Filter one
t = tiledlayout(2,1);
ax1 = nexttile;
hold on
plot(t_radar, filt_phase_08_50)

% Plot R-peaks
r_peaks = t_ECG(abs(qrs_i));
for i = r_peaks
    xline(i,'--k')
end

title('Filtered phase signal (8-50Hz), BB*, meas 031420')
xlabel('Time [s]');
ylabel('Phase [rad]');

% Filter two
ax2 = nexttile;
hold on
plot(t_radar, filt_phase_test)   

% plot R-peaks
r_peaks = t_ECG(abs(qrs_i));
for i = r_peaks
    xline(i,'--k')
end

title('Filtered phase signal (10-50Hz), BB*, meas 031420')
xlabel('Time [s]');
ylabel('Phase [rad]');
hold off

linkaxes([ax1, ax2], 'xy')

%% Plot R-peaks from ECG
r_peaks = t_ECG(abs(qrs_i));
for i = r_peaks
    xline(i,'--k')
end

%% Plot FFT of signal
signal = fft(filt_phase_08_50); % signal to do FFT on
Fs = 100;
L = num_frames;

P2 = abs(signal/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1)
title('FFT of HP filtered phase signal, FB*, meas 031414')
xlabel('Frequency [Hz]');
ylabel('Amplitude');

%% Plot detrended data based on found "change points"
ipt = findchangepts(unwrapped_phase,'Statistic','linear','MinThreshold',100);
%findchangepts(unwrapped_phase,'Statistic','linear','MinThreshold',100)

t = t_radar;
x = unwrapped_phase;
bp = t_radar(ipt);
y = detrend(x,1,bp,'SamplePoints',t,'Continuous',true);
plot(t,x,t,y,t,x-y,':k')
legend('Input Data','Detrended Data','Trend','Location','northwest') 

%% Plot data where deviating data points have been clipped out
% Should be altered and developed further for the spcified signal

plot(t_radar, filloutliers(filt_phase_08_50, 'clip', 'mean'));
%plot(t_radar,filloutliers(filt_phase, 'spline', 'mean'));

