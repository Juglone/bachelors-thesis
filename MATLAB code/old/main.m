% == Specificy the meas to be analysed ==
date = "2023_03_14";    % date of meas
meas_nr = "6";          % meas number
both_lanes = 0;         % 1 = yes, 0 = no

% == Load radar parameters ==
radar_parameters_C;     % write name of file with parameters used in meas

% This will give access to the following variables, write over them if
% needed but remember that they are derived from each other

num_frames;             % number of frames = number of chirps
num_samples;            % number of ADC samples per chirp               
T_f;                    % period of frames = period of chirps
T_s;                    % period of samples per chirp

% == Read meas data and start time of meas ==
% Radar
[radar_data, radar_start_time] = read_radar_data(date, meas_nr, num_samples, num_frames, both_lanes);

% ECG
[ECG_data, ECG_start_time] = read_ECG_data(date, meas_nr);

% == Truncate ECG meas to start at the same time as radar meas ==
t_meas = num_frames*T_f;    % Length of radar meas in seconds
T_ECG = 1/1000;             % Periodicity of ECG samples             

% If using the start times of measurements
ECG_trunc = truncate_ECG(radar_start_time, ECG_start_time, ECG_data, T_ECG, t_meas, 0);

% If the measurements have been ended at the same time
%num_i_ECG = floor(t_meas/T_ECG);    % Num of i that the meas length corresponds to
%ECG_trunc = ECG_data(end-num_i_ECG+1:end);


% == Perform Pan Tompkins on ECG data to detect R-peaks ==
%unique(isfinite(ECG_data)) % Use to check if ECG_data contains any NaNs
[qrs_amp, qrs_i_raw, delay] = pan_tompkin(ECG_trunc, 1000,0);
qrs_i = qrs_i_raw - delay;  % Shift due to delay in Pan Tompkin
%qrs_i = qrs_i_raw;          % Without shift due to delay


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

% Filter
%filter_08_50 = fir1(1000, [8/fn,50/fn]); % Use this filter
Fs = 100;                                                   % Sampling Frequency (Hz)
fcuts = [14.5 15 16 16.5 19.5 20 21 21.5 26 26.5 27.5 28];           % Frequencies
mags = [0 1 0 1 0 1 0];                                    % Passbands & Stopbands
devs = [0.05 0.01 0.05 0.01 0.05 0.01 0.05];                % Tolerances
[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fs_radar);          % Kaiser Window FIR Specification
n = n + rem(n,2);
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');           % Filter Realisation

filt_phase = filtfilt(hh, 1, unwrapped_phase);

% == Filter the unwrapped phase signal ==
% Apply highpass filter
%hp_filt_phase = filtfilt(hp_filter, 1, unwrapped_phase);
%filt_phase = filtfilt(lp_filter, 1, hp_filt_phase);

% Smooth the filtered signal to get rid of noise by using moving average
smooth_phase = smoothdata(filt_phase, 'movmean', 5);

% == Plot data ==
% Run this section first to load in the time vectors
t_ECG = T_ECG*(0:numel(ECG_data)-1);  % Time vector for trunc ECG data
t_radar = T_f*(0:num_frames-1);             % Time vector for chirps/frames

%% == Locate the peaks in the phase signal ==
[pks, locs] = findpeaks(smooth_phase);

% Inverted signal
inv_phase = -1*smooth_phase;
[inv_pks, inv_locs] = findpeaks(inv_phase);

%%
for i = possible_bins
    hold on
    plot(t_radar, unwrap(fft_phase(:,i)))
end

%% Plot phase signal before and after filter

subplot(2,1,1);
plot(t_radar, unwrapped_phase)      % Unfiltered phase
title('Unwrapped phase signal, FB')
xlabel('Time [s]');
ylabel('Phase [rad]');

subplot(2,1,2);
hold on
%plot(t_radar, filt_phase)            % HP filtered phase
plot(t_radar, smooth_phase)        % Filtered and smoothed phase
%plot(t_radar,filloutliers(filt_phase, 'clip', 'mean'));
r_peaks = t_ECG(abs(qrs_i));
for i = r_peaks
    xline(i,'--k')
end

% Plot R-peaks
%legend('Unfiltered', 'HP filtered', 'HP and LP filtered', 'Filtered and smoothed')
%legend('Phase signal', 'R-peaks')
title('Filtered phase signal (8-50Hz), FB')
xlabel('Time [s]');
ylabel('Phase [rad]');
hold off

%% Plot signal before and after filter
hold on
%plot(t_radar, unwrapped_phase)      % Unfiltered phase
%plot(t_radar, hp_filt_phase)    % HP filtered phase
%plot(t_radar, filt_phase)       % HP and LP filtered phase
plot(t_radar, smooth_phase)         % Filtered and smoothed phase
r_peaks = t_ECG(abs(qrs_i));
for i = r_peaks
    xline(i,'--k')
end
%legend('Unfiltered', 'HP filtered', 'HP and LP filtered', 'Filtered and smoothed')
%legend('Phase signal', 'R-peaks')
title('Filtered phase signal (1.2-1.6Hz), BB')
xlabel('Time [s]');
ylabel('Amplitude [arb. unit]');
hold off

%% Plot smoothed IBI
hold on
plot(t_ECG(qrs_i_trunc(1:end-1)), IBI_ECG, "k")
plot(t_radar(locs(1:end-1)),smoothdata(IBI_phase, 'movmean', 5), "r")
plot(t_radar(inv_locs(1:end-1)),smoothdata(IBI_inv_phase, 'movmean', 2), "b")
legend('ECG', 'Phase signal', 'Inverted phase signal')
title('Time between peaks, movmean 4 , FB')
xlabel('Time [s]');
ylabel('IBI [s]');

%% Plot IBI
hold on
plot(t_ECG(qrs_i_trunc(1:end-1)), IBI_ECG, "k")
plot(t_radar(locs(1:end-1)),IBI_phase, "r")
plot(t_radar(inv_locs(1:end-1)), IBI_inv_phase, "b")
legend('ECG', 'Phase signal', 'Inverted phase signal')
title('Time between peaks, FB')
xlabel('Time [s]');
ylabel('IBI [s]');


%%
hold on
plot(t_ECG(qrs_i_trunc(1:end-1)), 0, "ok")
plot(t_radar(locs(1:end-1)),0, "or")
plot(t_radar(inv_locs(1:end-1)),0, "ob")

%%
hold on
plot(t_ECG(qrs_i_trunc(1:end-1)), t_ECG(qrs_i_trunc(1:end-1)), "ok")
plot(t_radar(locs(1:end-1)),t_radar(locs(1:end-1)), "or")
plot(t_radar(inv_locs(1:end-1)),t_radar(inv_locs(1:end-1)), "ob")

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
ylabel('Phase displacement [rad]');

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
legend('Phase signal', 'Detected peaks','R-peaks', 'location', 'southeast')
title('Filtered phase signal (1,2 - 1,5Hz), FB')
xlabel('Time [s]');
ylabel('Amplitude [arb. unit]');
hold off

%% Plot inverted signal and detected peaks compared to R-peaks
plot(t_radar, inv_phase)     % Filtered and smoothed phase
hold on
plot(t_radar(inv_locs), 0, "or")   % Located peaks
plot(t_ECG(qrs_i_trunc), 0, "ok") % R-peaks
title('Inverted phase signal with detected peaks 1.2-2 Hz')
xlabel('Time [s]');
ylabel('Amplitude [arb. unit]');
hold off

%% == FFT on phase signal
Y = fft(filt_phase);
Fs = 100;
L = num_frames;

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1)
title('FFT of unwrapped phase signal, FB')
xlabel('Frequency [Hz]');
ylabel('Amplitude [arb. unit]');

%% STFT
hold on
%plot(t_radar, filt_I)
%plot(t_radar, filt_Q)
%plot(t_radar, s)
plot(t_radar, henv)
%plot(t_radar, normalize(smooth_phase))
%plot(t_radar, unwrapped_I)
%plot(t_radar, unwrapped_Q)
%plot(t_radar, normalize(unwrapped_phase))

%% STFT
hold on
x = henv;
fs = 100;
[a,b,c] = stft(x,fs,"Window", gausswin(10),"OverlapLength",2, "FFTLength", 800);
stft(x,fs,"Window", gausswin(10),"OverlapLength",2, "FFTLength", 40);
plot(t_ECG(qrs_i_trunc), 0, "ok") % R-peaks
%% Envelope
[up,lo] = envelope(smooth_phase);
hold on
plot(t_radar,up)
plot(t_radar, smooth_phase)   
%plot(t,filtered_phase,t,up,t,lo,'linewidth',1.5)
%legend('q','up','lo')
hold off

%% PSD
fs = 100;
t = t_radar;
x = unwrapped_phase;

N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(x):fs/2;

plot(freq,pow2db(psdx))
grid on
title("Periodogram Using FFT")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")
%% Average
IBI_average = zeros(numel(IBI_phase),1);
if locs(1) < inv_locs(1)
    IBI_average(1) = (IBI_phase(1)+IBI_inv_phase(1)/2)/1.5;
    for i = (2:numel(IBI_phase))
        IBI_average(i) = (IBI_phase(i)+(IBI_inv_phase(i-1)+ IBI_inv_phase(i))/2)/2;
    end
elseif locs(1) > inv_locs(1)
    IBI_average(1) = (IBI_inv_phase(1)+IBI_phase(1)/2)/1.5;
    for i = (2:numel(IBI_inv_phase))
        IBI_average(i) = (IBI_inv_phase(i)+(IBI_phase(i-1)+ IBI_phase(i))/2)/2;
    end
end
