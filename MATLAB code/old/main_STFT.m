% == Specificy the meas to be analysed ==
date = "2023_03_08";    % date of meas
meas_nr = "5";          % meas number
both_lanes = 1;         % 1 = yes, 0 = no

% == Load radar parameters ==
radar_parameters_D;     % write name of file with parameters used in meas

% This will give access to the following variables, write over them if
% needed but remember that they are derived from each other

num_frames;             % number of frames = number of chirps
num_samples;            % amount of ADC samples per chirp               
T_f;                    % period of frames = period of chirps
T_s;                    % period of samples per chirp

% == Read meas data and start time of meas ==
% Radar
[radar_data, radar_start_time] = read_radar_data(date, meas_nr, num_samples, num_frames, both_lanes);

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
t_meas = num_frames*T_f;                % Length of the radar meas
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

I = real(radar_data);
Q = imag(radar_data);

unwrapped_I = I(:,selected_bin);
unwrapped_Q = Q(:,selected_bin);

unwrapped_s = unwrapped_I + unwrapped_Q*sqrt(-1);

% == Define filters for the unwrapped phase signal ==
fs_radar = 1/T_f;           % Sampling frequency of chirps

% Highpass filter to filter out breathing below 1 Hz
f_hp = [7.6 8];             % Band edges of filter
amp_hp = [0 1];             % Band amplitude, highpass 
dev_hp = [0.01 0.01];     % Maximum deviation

[n_hp,Wn_hp,beta_hp,ftype_hp] = kaiserord(f_hp,amp_hp,dev_hp,fs_radar);
hp_filter = fir1(n_hp,Wn_hp,ftype_hp,kaiser(n_hp+1,beta_hp),'noscale');

% Lowpass filter to filter out high frequencies above 10 Hz
%f_lp = [40 40.4];             % Band edges of filter
%amp_lp = [1 0];             % Band amplitude, highpass 
%dev_lp = [0.01 0.01];       % Maximum deviation
%
%[n_lp,Wn_lp,beta_lp,ftype_lp] = kaiserord(f_lp,amp_lp,dev_lp,fs_radar);
%lp_filter = fir1(n_lp,Wn_lp,ftype_lp,kaiser(n_lp+1,beta_lp),'noscale');

% == Filter the unwrapped phase signal ==
% Apply highpass filter

%filt_phase = filtfilt(lp_filter, 1, unwrapped_phase);
%filt_I = filtfilt(lp_filter, 1, unwrapped_I);
%filt_Q = filtfilt(lp_filter, 1, unwrapped_Q);

filt_phase = filtfilt(hp_filter, 1, unwrapped_phase);
filt_I = filtfilt(hp_filter, 1, unwrapped_phase);
filt_Q = filtfilt(hp_filter, 1, unwrapped_phase);

s_filt = filt_I + sqrt(-1)*filt_Q;
sig = sqrt(filt_I.^2 + filt_Q.^2);

f_s = [20 20.4];             % Band edges of filter
amp_s = [1 0];             % Band amplitude, highpass 
dev_s = [0.01 0.01];       % Maximum deviation

[n_s,Wn_s,beta_s,ftype_s] = kaiserord(f_s,amp_s,dev_s,fs_radar);
s_filter = fir1(n_s,Wn_s,ftype_s,kaiser(n_s+1,beta_s),'noscale');

env = filt_phase;

filt_s = filtfilt(s_filter, 1, log(env));
henv = exp(filt_s);

% == Plot data ==
% Run this section first to load in the time vectors
t_ECG = T_ECG*(0:numel(ECG_data_trunc)-1);  % Time vector for trunc ECG data
t_radar = T_f*(0:num_frames-1);             % Time vector for chirps/frames


%% STFT new, run this to perform and plot STFT
hold on
x = sig;
fs = 200;
[a,b,c] = stft(x,fs,"Window", hamming(20),"OverlapLength",4, "FFTLength", 100);
stft(x,fs,"Window", hamming(20),"OverlapLength",4, "FFTLength", 100)
%plot(t_ECG(qrs_i_trunc), 0, "ok") % R-peaks

%% Integrate, filter and plot IBI from STFT
stft_pos = a((58:end-1),:);
stft_neg = flip(a((1:42),:));

mean_pos = zeros(1,size(stft_pos,2));
mean_neg = zeros(1,size(stft_pos,2));

for i = 1:size(stft_pos,2)
    pos_weights = 0;
    pos_sum = 0;
    neg_weights = 0;
    neg_sum = 0;
    for j = 1:size(stft_pos,1)
        pos_weights = pos_weights + stft_neg(j,i);
        pos_sum = j*stft_neg(j,i);
        neg_weights = pos_weights + stft_neg(j,i);
        neg_sum = j*stft_neg(j,i);
    end
    mean_pos(i) = pos_sum/pos_weights;
    mean_neg(i) = neg_sum/neg_weights;
end

mean_sum = mean_pos + mean_neg;
mean_abs = abs(mean_sum);
mean_smooth = smoothdata(mean_abs, 'gaussian', 5);
hold on
%plot(c, mean_pos)
%plot(c, mean_neg)
plot(c, mean_smooth)
r_peaks = t_ECG(qrs_i_trunc);
for i = r_peaks
    xline(i,'--k')
end

% == Locate the peaks in the phase signal ==
[pks, locs] = findpeaks(mean_smooth);

T_mean = T_f*num_frames/numel(c);
% == Analyse IBI ==
% Analyse IBI from R-waves 
IBI_ECG = zeros(numel(qrs_i_trunc)-1,1);
for i = 2:numel(qrs_i_trunc)
    i_diff = qrs_i_trunc(i)-qrs_i_trunc(i-1);
    IBI_ECG(i-1) = i_diff*T_ECG;
end

% Analyse IBI from ordinary phase signal
IBI_phase = zeros(numel(pks)-1,1);
for i = 2:numel(pks)
    i_diff = locs(i)-locs(i-1);
    IBI_phase(i-1) = i_diff*T_mean;
end

%% Plot signal and detected peaks compared to R-peaks
hold on
plot(c, mean_smooth)     % Filtered and smoothed phase
plot(c(locs), 0, "or")   % Located peaks
plot(t_ECG(qrs_i_trunc), 0, "ok") % R-peaks
legend('Phase signal', 'Detected peaks','R-peaks', 'location', 'southeast')
title('Filtered phase signal (1,2 - 1,5Hz), FB')
xlabel('Time [s]');
ylabel('Amplitude [arb. unit]');
hold off

%% Plot IBI
hold on
plot(t_ECG(qrs_i_trunc(1:end-1)), IBI_ECG, "k")
plot(c(locs(1:end-1)),IBI_phase, "r")
legend('ECG', 'Phase signal')
title('Time between peaks, FB')
xlabel('Time [s]');
ylabel('IBI [s]');


%% Plot smoothed IBI
hold on
plot(t_ECG(qrs_i_trunc(1:end-1)), IBI_ECG, "k")
plot(c(locs(1:end-1)),smoothdata(IBI_phase, 'movmean', 5), "r")
legend('ECG', 'Phase signal')
title('Time between peaks, movmean 4 , FB')
xlabel('Time [s]');
ylabel('IBI [s]');


%% STFT
hold on
x = s_filt;
fs = 100;
[a,b,c] = stft(x,fs,"Window", gausswin(10),"OverlapLength",2, "FFTLength", 400);
stft(x,fs,"Window", gausswin(20),"OverlapLength",10, "FFTLength", 100);
plot(t_ECG(qrs_i_trunc), 0, "ok") % R-peaks

b(1:231,:) = [];
a(1:231,:) = []; 

stft_mean = mean(a);

%%
hold on
plot(c, abs(stft_mean))
plot(t_ECG(qrs_i_trunc), 0, "ok") % R-peaks

%%
hold on
plot(t_radar, normalize(filt_I))
plot(t_ECG(qrs_i_trunc), 0, "ok") % R-peaks

%%
hold on
%plot(t_radar, unwrapped_I)
%plot(t_radar, unwrapped_Q)
plot(t_radar, unwrapped_s)

%% Plot smoothed IBI
hold on
plot(t_ECG(qrs_i_trunc(1:end-1)), IBI_ECG, "k")
plot(t_radar(locs(1:end-1)),smoothdata(IBI_phase, 'movmean', 2), "r")
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
ylabel('Amplitude [arb. unit]');


%% Plot signal before and after filter
hold on
%plot(t_radar, unwrapped_phase)      % Unfiltered phase
%plot(t_radar, hp_filtered_phase)    % HP filtered phase
plot(t_radar, filt_phase)       % HP and LP filtered phase
%plot(t_radar, smooth_phase)         % Filtered and smoothed phase
%plot(t_ECG(qrs_i_trunc), 0, "ok") % R-peaks
%legend('Unfiltered', 'HP filtered', 'HP and LP filtered', 'Filtered and smoothed')
%legend('Phase signal', 'R-peaks')
title('Filtered phase signal (1.2-1.6Hz), BB')
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
plot(t_ECG(qrs_i_trunc), 0, "|k") % R-peaks
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
Y = fft(log(env));
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
plot(t_radar, normalize(filt_Q))
plot(t_ECG(qrs_i_trunc), 0, "ok") % R-peaks
%plot(t_radar, filt_Q)
%plot(t_radar, s)
%plot(t_radar, henv)
%plot(t_radar, normalize(smooth_phase))
%plot(t_radar, unwrapped_I)
%plot(t_radar, unwrapped_Q)
%plot(t_radar, normalize(unwrapped_phase))
%plot3(unwrapped_I, unwrapped_Q,t_radar)
%plot(t_radar, log(filt_phase))

%% STFT
hold on
x = filt_phase;
fs = 100;
[a,b,c] = stft(x,fs,"Window", gausswin(10),"OverlapLength",2, "FFTLength", 400);
stft(x,fs,"Window", gausswin(16),"OverlapLength",2, "FFTLength", 400);
plot(t_ECG(qrs_i_trunc), 0, "ok") % R-peaks
%% Envelope
[up,lo] = envelope(filt_phase);
hold on
plot(t_radar,smoothdata(up))
%plot(t_radar, smooth_phase)   
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
