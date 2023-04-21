% == Specificy the meas to be analysed and load parameters ==
date = "2023_04_12";    % date of meas
meas_nr = "15";         % meas number
both_lanes = 1;         % was data collected from both lanes, 1 = yes, 0 = no

radar_parameters_E;     % name of file with parameters used in meas

% This should give access to the following variables
num_frames;             % number of frames = number of chirps
num_samples;            % number of ADC samples per chirp              
T_f;                    % period of frames = period of chirps
T_s;                    % period of samples per chirp
t_meas;                 % length of radar meas in seconds

T_ECG = 1/1000; % period of ECG samples

% == Read meas data and start time of meas ==
% Comment: the functions that read the data need to be adjusted to work for 
% someone elses data paths
% Comment: the start time is no longer used to sync ECG and radar

% Radar
[radar_data, radar_start_time] = read_radar_data(date, meas_nr, num_samples, num_frames, both_lanes);

% ECG
[ECG_data, ECG_start_time] = read_ECG_data(date, meas_nr);

% == Truncate ECG meas to sync with radar meas ==             
% Use the following line if using start times (old)
%ECG_trunc = truncate_ECG(radar_start_time, ECG_start_time, ECG_data, T_ECG, t_meas, 1.15);

% Use the following lines if both meas were stopped at the same time
num_i_ECG = floor(t_meas/T_ECG);            % num of i corresponding to meas length
del_i_ECG = floor(0/T_ECG);          % additional time to take away from end of ECG meas
ECG_trunc = ECG_data(end-num_i_ECG-del_i_ECG+1:end-del_i_ECG);  % truncated ECG data

t_ECG = T_ECG*(0:numel(ECG_data)-1);    % time vector for trunc ECG data
t_radar = T_f*(0:num_frames-1);         % time vector for chirps/frames


% == Perform Pan Tompkins on ECG data to detect R-peaks ==
%unique(isfinite(ECG_data)) % use to check if ECG_data contains any NaNs
[qrs_amp, qrs_i_raw, delay] = pan_tompkin(ECG_trunc, 1000,0);
qrs_i = qrs_i_raw - delay;  % shift indexes due to delay in Pan Tompkins
%qrs_i = qrs_i_raw;         % without shift due to delay

% == Perform range FFT on radar data ==
% Performs an FFT on each row (each chirp) in our radar data
range_fft = fft(radar_data,[],2);

fft_magnitude = abs(range_fft); % magnitudes of each FFT
fft_phase = angle(range_fft);   % angles of each FFT


% == Find and select range bin(s) == 
% Comment: we want to develop this method
% To plot detected peaks, switch the '0' to '1'
[bin_possible, bin_selected] = find_range_bins(num_frames, num_samples, T_f, fft_magnitude, 0);

% == Extract and unwrap the phase signal ==
% Comment: currently we extract the phase signal based on the most frequent 
% range bin
wrapped_phase = fft_phase(:,bin_selected);  % wrapped phase signal
unwrapped_phase = unwrap(wrapped_phase);    % the unwrapped phase signal 
unwrapped_test = unwrap(fft_phase(:,bin_selected+1));  % wrapped phase signal

% == Create filters ==
fs_radar = 1/T_f;           % sampling frequency
fn = 0.5/T_f;               % nyquist frequency

% Filter one
f_1 = [7.6 8];          % band edges of filter
amp_1 = [0 1];          % band amplitude
dev_1 = [0.01 0.01];    % maximum deviation
[n_1,Wn_1,beta_1,ftype_1] = kaiserord(f_1,amp_1,dev_1,fs_radar);
filter_1 = fir1(n_1,Wn_1,ftype_1,kaiser(n_1+1,beta_1),'noscale');

% Filter two - additional filter to play around with
f_2 = [9.6 10 30 30.4];         % band edges of filter
amp_2 = [0 1 0];         % band amplitude, highpass 
dev_2 = [0.01 0.01 0.01];    % maximum deviation
[n_2,Wn_2,beta_2,ftype_2] = kaiserord(f_2,amp_2,dev_2,fs_radar);
filter_2 = fir1(n_2,Wn_2,ftype_2,kaiser(n_2+1,beta_2),'noscale');

signal_to_filt = unwrapped_phase;

% == Filter the unwrapped phase signal ==
% Apply highpass filter
filt_phase_08_50 = filtfilt(filter_1, 1, signal_to_filt);
filt_phase_test = filtfilt(filter_2, 1, signal_to_filt);

% Smooth the filtered signal by using moving average
smooth_phase_08_50 = smoothdata(filt_phase_08_50, 'movmean', 5);
smooth_phase_test = smoothdata(filt_phase_test, 'movmean', 5);

%% Find most interesting range bin by using spatial phase coherency (SPC)

% Unwrap all phase signals
p_unwrapped = zeros(num_frames, num_samples);
for i = 1:25
    p_unwrapped(:,i) = unwrap(fft_phase(:,i));
end

% Use wavelets to decompose signal 2 levels, save d1 and d2 seperately
p_d1 = zeros(num_frames, num_samples);
p_d2 = zeros(num_frames, num_samples);
for i = 1:num_samples
    [p_a1(:,i), p_d1(:,i)] = wavelet_leveldepend(p_unwrapped(:,i),'db6',1,0,10,0);
    [p_a2(:,i), p_d2(:,i)] = wavelet_leveldepend(p_unwrapped(:,i),'db6',2,0,20,0);
end

% Normalize the signals
std_d1 = std(p_d1);
p_d1_norm = p_d1./std_d1;

std_d2 = std(p_d2);
p_d2_norm = p_d2./std_d2;


% Calculate SPC for d1 and d2 of all bins
SPC_d1 = zeros(25);
SPC_d2 = zeros(25);
for i = 1:25
    for j = 1:25
        if i == j
            SPC_d1(i,j) = 0;
            SPC_d2(i,j) = 0;
        else
            SPC_d1(i,j) = sum(p_d1_norm(:,i).*p_d1_norm(:,j));
            SPC_d2(i,j) = sum(p_d2_norm(:,i).*p_d2_norm(:,j));
        end
    end
end
%% Plot SPC of filtered signal vs filtered and despiked signal 
[X,Y] = meshgrid(1:25,1:25);

figure
surface(X,Y,SPC_d1)
colorbar
title('SPC of d1')

figure
surface(X,Y,SPC_d2)
colorbar
title('SPC of d2')

%% Select main bin of interest
% Automatically pick bin of interest by computing the sum of the SPC for
% each bin and choosing the bin of highest peak
% d2 seems to show better snr so we use the SCP for d2, guessing that since 
% d1 contains higher frequencies it might be comparable to noise
[foo, b_max] = max(sum(SPC_d2));

% Manually pick main bin of interest
%b_max = 

hold on
plot(sum(SPC_d1))
plot(sum(SPC_d2)) % we find the highest peak in this plot

%% Select other bins of interest 
% Specify the SPC for the main bin of interest
SPC_max = SPC_d2(b_max,:);

% Automatically choose k bins of highest peaks
[foo, b_interest] = maxk(SPC_max,4);

% Add b_max to the bins of interest
b_interest(end+1) = b_max;

% Sort the bins from lowest to highest
b_interest = sort(b_interest);

% Manually select other bins of interest by looking at the plot
%b_interest = [b_max,4,10,11];

hold on
plot(SPC_max) % we fin the k highest peaks in this plot

%% Plot the decomposed signals of the bins of interest
for i = 1:numel(b_interest)
    subplot(numel(b_interest), 1, i)
    hold on
    plot(t_radar, p_d1_norm(:,b_interest(i)))
    plot(t_radar, p_d2_norm(:,b_interest(i)))
    r_peaks = t_ECG(abs(qrs_i));
    for j = r_peaks
        xline(j,'--k')
    end
    hold off
    legend('d1','d2')
    title(['Decomposed signal from bin ' num2str(b_interest(i))])
end

%% Plot absolute values of d1 and d2 
for i = 1:numel(b_interest)
    subplot(numel(b_interest), 1, i)
    hold on
    plot(t_radar,movmean(abs(p_d1_norm(:,b_interest(i))),1))
    plot(t_radar,movmean(abs(p_d2_norm(:,b_interest(i))),1))
    r_peaks = t_ECG(abs(qrs_i));
    for j = r_peaks
        xline(j,'--k')
    end
    hold off
    legend('abs(d1)','abs(d2)')
    title(['Absolute value of decomposed signal from bin ' num2str(b_interest(i))])
end


%% Compute the product between the absolute values of d1 and d2
p_d1_abs = movmean(abs(p_d1_norm),1);
p_d2_abs = movmean(abs(p_d2_norm),1);
SPC_d1d2 = zeros(num_frames,numel(b_interest));
for i = 1:numel(b_interest)
    product = movmean(p_d1_abs(:,b_interest(i)).*p_d2_abs(:,b_interest(i)),1);
    SPC_d1d2(:,i) = product;
end

hold on
plot(t_radar, movmean(sum(SPC_d1d2,2),10))
%plot(t_radar, movmean(SPC_d1d2,10))
r_peaks = t_ECG(abs(qrs_i));
for i = r_peaks
    xline(i,'--k')
end
hold off
%%
signal = movmean(SPC_d1d2,10);
pks = cell(numel(b_interest), 1);
locs = cell(numel(b_interest), 1);
for i = 1:numel(b_interest)
    [pks{i}, locs{i}] = findpeaks(signal(:,i),'MinPeakDistance',50);
end

for i = 1:numel(b_interest)
    subplot(numel(b_interest), 1, i)
    hold on
    plot(t_radar, signal(:,i))
    plot(t_radar(locs{i}), pks{i},'o')
    r_peaks = t_ECG(abs(qrs_i));
    for j = r_peaks
        xline(j,'--k')
    end
    hold off
    title(['Found peaks from bin ' num2str(b_interest(i))])
end



%% Find peaks of signal

% Moving average of window with 10 samples seems to be pretty good
signal = movmean(sum(SPC_d1d2(:,4),2),10);

[pks, locs] = findpeaks(signal,'MinPeakDistance',50);

% Remove certain peaks if necessary 
remove = [];
pks(remove) = [];
locs(remove) = [];

% Plot the found peaks compared to ECG
hold on
plot(t_radar, signal)
plot(t_radar(locs), pks,'o')

r_peaks = t_ECG(abs(qrs_i));
for i = r_peaks
    xline(i,'--k')
end
hold off

%% == Analyse IBI ==
% Analyse IBI from R-waves 
IBI_ECG = zeros(numel(qrs_i)-1,1);
for i = 2:numel(qrs_i)
    i_diff = qrs_i(i)-qrs_i(i-1);
    IBI_ECG(i-1) = i_diff*T_ECG;
end

% Analyse IBI from phase signal
IBI_phase = zeros(numel(pks)-1,1);
for i = 2:numel(pks)
    i_diff = locs(i)-locs(i-1);
    IBI_phase(i-1) = i_diff*T_f;
end

figure
hold on
plot(t_ECG(abs(qrs_i(1:end-1))), IBI_ECG, "k")
plot(t_radar(locs(1:end-1)),movmean(IBI_phase,2), "r")

%% Plot FFT of signal
signal = fft(signal); % signal to do FFT on
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

%% Envelope

%p_new_1 = spk.*abs(spk);
%p_new_1 = p_max.*p_max;
%p_new_1 = p_sum_despiked.*p_sum_despiked;
%p_new_1 = wavelet_filt.*wavelet_filt;
p_new_1 = p_spikes.*abs(p_spikes);

p_new_2 = p_new_1.*abs(p_new_1);
p_new_3 = p_new_2.*abs(p_new_2);
hold on
%plot(p_new_2)
plot(t_radar,p_new_1)
%plot(p_sum_despiked)

r_peaks = t_ECG(abs(qrs_i));
for i = r_peaks
    xline(i,'--k')
end
hold off

%% Wavelet 
[despk, spk] = wavelet_despike_movrms(p_spikes,'haar',5,2,50,1);

%%
[p_env_up, p_env_low] = envelope(p_spikes_sum,1,'peak');
movrms = sqrt(movmean(p_spikes .^ 2, 100));

%hold on
%envelope(p_spikes,5,'peak');
hold on
plot(p_env_up)
plot(new_sig)

%%

hold on

plot(t_radar, p_env_up)
plot(t_radar, -p_env_low)
r_peaks = t_ECG(abs(qrs_i));
for i = r_peaks
    xline(i,'--k')
end
hold off

%% Locate peaks from envelope

[pks, locs] = findpeaks(p_spikes_sum);
hold on
plot(t_radar(locs), 0,'o')

r_peaks = t_ECG(abs(qrs_i));
for i = r_peaks
    xline(i,'--k')
end
hold off

%% Test to get rid of high peaks by subtracting overlap of different bins
bin1 = 4;
bin2 = 6;
test0 = zeros(num_frames,1);
int_sig = (p_filtered(:,bin1)/std_filtered(bin1));
comb_sig = (p_filtered(:,bin1)/std_filtered(bin1)+p_filtered(:,bin2)/std_filtered(bin2))./2;
err_sig = p_filtered(:,bin2)/std_filtered(bin2);
SPC_terms = (p_filtered(:,bin1).*p_filtered(:,bin2))./(std_filtered(bin1)*std_filtered(bin2));
for i = 1:num_frames
    if abs(SPC_terms(i)) >= 2
        test0(i) = err_sig(i);
    end
end
hold on 
%plot(SPC_terms)
%plot(test0)
%plot(err_sig)
%plot(comb_sig)
plot(int_sig-test0)
%plot(comb_sig-test0-20)

%% Wavelet analysis to filter out baseline drift

ssds = zeros(3,1);

cur_lp = filt_phase_08_50;

iterations = 0;
counter = 0; 

name = 'haar';
%name = 'sym8';

while 1
    % Decompose one level
    [lp,hp] = dwt(cur_lp,name);
    
    % Shift and calculate energy of detail/high pass coefficient
    ssds = vertcat(sum(hp.^2), ssds);
    
    if (ssds(3) > ssds(2)) && (ssds(2) < ssds(1))
        break
    end
    cur_lp = lp;
    iterations = iterations + 1; 
    counter = counter + 1; 
        
end
baseline = cur_lp;

while 1
    baseline = idwt(baseline, zeros(numel(baseline),1), name);
    if iterations == 1
        break
    end
    iterations = iterations -1;
end

new_baseline = resample(baseline, numel(unwrapped_phase), numel(baseline));

%%
hold on
plot(unwrapped_phase)
plot(new_baseline)
plot(unwrapped_phase-new_baseline)

%% Detrend data based on found "change points"
chpts = findchangepts(unwrapped_phase,'Statistic','linear','MinThreshold',10);
findchangepts(unwrapped_phase,'Statistic','linear','MinThreshold',10)
%%
t = t_radar;
x = unwrapped_phase;
bp = t_radar(chpts);
y = detrend(x,1,bp,'SamplePoints',t,'Continuous',true);
plot(t,x,t,y,t,x-y,':k')
legend('Input Data','Detrended Data','Trend','Location','southwest') 

detrended_phase = y;


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
hold on
plot(t_radar, p_detrended(:,17))      % unfiltered phase
title('Unwrapped phase signal, BB*, meas 031420')
xlabel('Time [s]');
ylabel('Phase [rad]');

% Filtered signal
ax2 = nexttile;
hold on
%plot(t_radar, filt_phase)          % filtered phase
plot(t_radar, filt_phase_08_50)      % filtered and smoothed phase

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
plot(t_radar, phase_4)

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
plot(t_radar, phase_4-phase_delete_3)   

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
chpts = findchangepts(filt_phase_test,'Statistic','rms');
%findchangepts(unwrapped_phase,'Statistic','linear','MinThreshold',100)

t = t_radar;
x = filt_phase_test;
bp = t_radar(chpts);
y = detrend(x,1,bp,'SamplePoints',t,'Continuous',true);
plot(t,x,t,y,t,x-y,':k')
legend('Input Data','Detrended Data','Trend','Location','northwest') 

detrended_phase_test = y;

%%
hold on
plot(filt_phase_test+1000)
plot(detrended_phase_test)

%% Plot data where deviating data points have been clipped out
% Should be altered and developed further for the spcified signal
hold on
plot(t_radar, filt_phase_08_50)
plot(t_radar, filloutliers(filt_phase_08_50, 'linear','percentiles', [1 99]))
%%
hold on
plot(t_radar, filt_phase_test-500)
plot(t_radar, filloutliers(filt_phase_test, 'clip', 'median')+500);
plot(t_radar, filloutliers(filt_phase_test, 'spline', 'movmean',100,'Samplepoints', t_radar));
%plot(t_radar,filloutliers(filt_phase, 'spline', 'mean'));

%% Plotta en veckad fassignal vs uppveckad fassignal

subplot(1,2,1)
hold on
plot(wrapped_phase)
yline(pi,'--k')
yline(-pi,'--k')
title('Veckad fassignal')

subplot(1,2,2)
hold on
plot(unwrapped_phase)
yline(pi,'--k')
legend('fassignal','\pm\pi')
title('Uppveckad fassignal')
%plot(-pi*ones(num_frames),'--k')
