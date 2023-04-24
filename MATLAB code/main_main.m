% == Specificy the meas to be analysed and load parameters ==
date = "2023_03_08";    % date of meas
meas_nr = "3";          % meas number
both_lanes = 1;         % was data collected from both lanes, 1 = yes, 0 = no
use_start_times = 0;    % 1 = yes, 0 = no, for new meas this should be 0

radar_parameters_B;     % name of file with parameters used in meas

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
[radar_data, radar_start_time] = read_radar_data(date, meas_nr, num_samples, num_frames, both_lanes, use_start_times);

% ECG
[ECG_data, ECG_start_time] = read_ECG_data(date, meas_nr, use_start_times);


% == Truncate ECG meas to sync with radar meas == 
t_lag = 0; 
% Specify additional lag in order to sync the measurents. This will cut off
% the end of the ECG meas. 
% In the case of using start times this will instead take away more from 
% the beginning of the ECG meas   
if use_start_times
    ECG_trunc = truncate_ECG(radar_start_time, ECG_start_time, ECG_data, T_ECG, t_meas, t_lag);
else
    num_i_ECG = floor(t_meas/T_ECG);        % num of i corresponding to meas length
    del_i_ECG = floor(t_lag/T_ECG);         % additional time to take away from end of ECG meas
    ECG_trunc = ECG_data(end-num_i_ECG-del_i_ECG+1:end-del_i_ECG);  % truncated ECG data
end

t_ECG = T_ECG*(0:numel(ECG_data)-1);    % time vector for trunc ECG data
t_radar = T_f*(0:num_frames-1);         % time vector for chirps/frames


% == Perform Pan Tompkins on ECG data to detect R-peaks ==
[qrs_amp, qrs_i_raw, delay] = pan_tompkin(ECG_trunc, 1000,0);
qrs_i = qrs_i_raw - delay;  % shift indexes due to delay in Pan Tompkins


% == Perform range FFT on radar data ==
% Performs an FFT on each row (each chirp) in our radar data
range_fft = fft(radar_data,[],2);

fft_magnitude = abs(range_fft); % magnitudes of each FFT
fft_phase = angle(range_fft);   % angles of each FFT


% == Unwrap all phase signals == 
% Choose how many bins to include in the analysis
%num_bins = num_samples;
num_bins = 25;

p_unwrapped = zeros(num_frames, num_bins);
for i = 1:num_bins
    p_unwrapped(:,i) = unwrap(fft_phase(:,i));
end


% == Filter signals for SPC ==
% Filter one
fs_radar = 100;
f_1 = [24.5 25];          % band edges of filter
amp_1 = [0 1];          % band amplitude
dev_1 = [0.01 0.01];    % maximum deviation
[n_1,Wn_1,beta_1,ftype_1] = kaiserord(f_1,amp_1,dev_1,fs_radar);
filter_1 = fir1(n_1,Wn_1,ftype_1,kaiser(n_1+1,beta_1),'noscale');


p_filtered = zeros(num_frames, num_bins);
for i = 1:num_bins
    p_filtered(:,i) = filtfilt(filter_1, 1, p_unwrapped(:,i));
end

% == Perform wavelet decomposition == 
% Use wavelets to decompose signal 2 levels, save d1 and d2 seperately
p_d1 = zeros(num_frames, num_bins);
p_d2 = zeros(num_frames, num_bins);
for i = 1:num_bins
    [p_a1(:,i), p_d1(:,i)] = wavelet_leveldepend(p_unwrapped(:,i),'db6',1,0,10,0);
    [p_a2(:,i), p_d2(:,i)] = wavelet_leveldepend(p_unwrapped(:,i),'db6',2,0,10,0);
end


% == Zero the endpoints of the signals ==
for i = 1:num_bins
    num_zeros = 10;     % choose how many frames to zero at the ends of signal
    p_d1([1:num_zeros, end-num_zeros:end], i)= 0;
    p_d2([1:num_zeros, end-num_zeros:end], i)= 0;
    p_filtered([1:num_zeros, end-num_zeros:end], i)= 0;
end


% == Normalize the signals (details) using std ==
% detail 1
std_d1 = std(p_d1);
p_d1_norm = p_d1./std_d1;

% detail 2
std_d2 = std(p_d2);
p_d2_norm = p_d2./std_d2;

% filtered
std_filt = std(p_filtered);
p_filt_norm = p_filtered./std_filt;

% == Calculate SPC for d1 and d2 of all bins ==
SPC_d1 = zeros(num_bins);
SPC_d2 = zeros(num_bins);
SPC_filtered = zeros(num_bins);
for i = 1:num_bins
    for j = 1:num_bins
        if i == j
            SPC_d1(i,j) = 0;
            SPC_d2(i,j) = 0;
            SPC_filtered(i,j) = 0;
        else
            SPC_d1(i,j) = sum(p_d1_norm(:,i).*p_d1_norm(:,j));
            SPC_d2(i,j) = sum(p_d2_norm(:,i).*p_d2_norm(:,j));
            SPC_filtered(i,j) = sum(p_filt_norm(:,i).*p_filt_norm(:,j));
        end
    end
end

SPC_tot = SPC_d1 + SPC_d2;

figure
surface(X,Y,SPC_filtered)
colorbar
title('SPC of filtered signals')

%% Plot SPC of filtered signal vs filtered and despiked signal 
% This step is not necessary for the analysis
[X,Y] = meshgrid(1:num_bins,1:num_bins);

figure
surface(X,Y,SPC_d1)
colorbar
title('SPC of d1')

figure
surface(X,Y,SPC_d2)
colorbar
title('SPC of d2')

figure
surface(X,Y,SPC_tot)
colorbar
title('Sum of SPC d1 and d2')

%% Select main bin of interest
% Automatically pick bin of interest by computing the sum of the SPC for
% each bin and choosing the bin of highest peak
% d2 seems to show better snr so we use the SCP for d2, guessing that since 
% d1 contains higher frequencies it might be comparable to noise
[foo, b_max] = max(sum(SPC_filtered));

% Manually pick main bin of interest based on plots
%b_max = 

hold on
plot(sum(SPC_filtered))
plot(sum(SPC_d1))
plot(sum(SPC_d2))
plot(sum(SPC_tot)) % we find the highest peak in this plot
legend('filtered','sum of SPC d1', 'sum of SPC d2', 'sum of SPC d1 and d2')

%% Select other bins of interest 
% Specify the SPC for the main bin of interest
SPC_max = SPC_tot(b_max,:);

% Automatically choose k bins of highest peaks
[foo, b_interest] = maxk(SPC_max,2);

% Add b_max to the bins of interest
b_interest(end+1) = b_max;

% Sort the bins from lowest to highest
b_interest = sort(b_interest);

% Manually select other bins of interest by looking at the plot
%b_interest = [b_max,4,10,11];

hold on
plot(SPC_max) % we find the k highest peaks in this plot

%% == Plot d1 and d2 of the bins of interest ==
% This step is not necessary for the analysis
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

%% == Plot absolute values of d1 and d2 ==
% This step is not necessary for the analysis
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


%% == Compute the product between the absolute values of d1 and d2 ==
p_d1_abs = movmean(abs(p_d1_norm),1);
p_d2_abs = movmean(abs(p_d2_norm),1);

SPC_d1d2 = zeros(num_frames,numel(b_interest));
for i = 1:numel(b_interest)
    product = movmean(p_d1_abs(:,b_interest(i)).*p_d2_abs(:,b_interest(i)),10);
    SPC_d1d2(:,i) = product;
end


hold on
plot(t_radar, movmean(SPC_d1d2,1))
title('Product between abs values of d1 and d2') 

r_peaks = t_ECG(abs(qrs_i));
for i = r_peaks
    xline(i,'--k')
end
hold off

%% == Find and plot peaks of individual bins ==
% This step is not necessary for the analysis
signal = movmean(SPC_d1d2,20);
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

%% == Find and plot peaks of sum of all bins ==
% Moving average of window with 10 samples seems to be pretty good
signal = movmean(sum(SPC_d1d2,2),20);

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
plot(t_radar(locs(1:end-1)),movmean(IBI_phase,4), "r")


%% == Interpolate IBI from ECG and calculate accuracy ==
t_IBI_ECG = qrs_i(1:end-1).*T_ECG;  % time vector of the IBI from ECG
t_IBI_radar = locs(1:end-1).*T_f;   % time vector of the IBI from radar

IBI_ECG_interp = interp1(t_IBI_ECG,IBI_ECG,t_IBI_radar,'linear','extrap');

% choose window length for movmean
IBI_radar = movmean(IBI_phase,4);

nog = 0;
for i = 1:numel(IBI_radar)
    nog(i) = (IBI_ECG_interp(i)-abs(IBI_ECG_interp(i)-IBI_radar(i)))./IBI_ECG_interp(i);
end

figure
hold on
plot(IBI_ECG_interp)
plot(IBI_radar)
title('IBI of radar vs interpolated IBI of ECG')
legend('ECG','Radar')

medelv = mean(nog)
media = median(nog)

%% == Plot FFT of signal ==
% This step is not necessary for the analysis
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

