%%
% Run the code from this file. Note that this code may contain errors
%%

clc
clear all
close all

N_sample=96; % number of samples
N_frame=6000;  % number of frames = number of chirps
delta_t=40*1e-3; % time between two chirps
fsample=1/delta_t; % sampling frequency
fn=fsample/2; % nyquist frequency

name_file='Meas1_R_11_NP_JB_2';
path_file=strcat('C:\Users\julia\OneDrive\Dokument\Chalmers\År 3\Kandidatarbete\Data\', name_file,'.bin');
ADCdata=readDCA1000(path_file,N_sample);

%data=zeros(N_frame,1);
mag=zeros(N_frame,N_sample);
%range=zeros(N_frame,N_sample);
phase=zeros(N_frame,N_sample);
%realsignal=zeros(N_frame,N_sample);
%pos=zeros(N_frame,1);

for p=1:N_frame
    chirp=ADCdata(2*(p-1)*N_sample+1:2:p*2*N_sample);  % the data of each chirp

    spec=fft(chirp);  % fft of one chirp
    %range(p,:)= spec;

    mag(p,:)=abs(spec);  % magnitude of fft
    phase(p,:)=angle(spec);  % phase of fft

end

t=0:delta_t:(N_frame-1)*delta_t; % time array

% select range bin(s) and get an unwraped phase signal.
[phase_unwraped, finalRangeBins, range, m, bin, L] = selectRangeBin_4(mag, phase, N_frame);


%% ECG measurments
name_file_ecg='Meas1_R_11_NP_JB_2_EKG';
path_file_ecg=strcat('C:\Users\julia\OneDrive\Dokument\Chalmers\År 3\Kandidatarbete\Data\EKG\', name_file_ecg,'.txt');

[sample_nr, voltageECG_with_offset] = importECG(path_file_ecg); % import the ECG data
f_sample_ECG = 1000; % sampling frequency of ECG, to be modified if wanted different
delta_tECG = 1/f_sample_ECG; % time difference between two samples

% remove DC-offset:
voltageECG = voltageECG_with_offset-mean(voltageECG_with_offset);

%synchronize ECG data with data fram radar
% this can only be used if ECG measurement ended at the same time as the
% radar measurement ended and if the ECG measurement was started before the
% radar measurement
[ECG_signal, sample_nr_signal] = synchronize_ECG_Radar(voltageECG,sample_nr,delta_tECG,delta_t ,N_frame);
tECG = delta_tECG.*sample_nr_signal; % end time for the ECG signal

% Plot ECG signal
figure(1)
plot(tECG,ECG_signal)

% pan tompkin algorithm to find IBI (IBI = inter beat interval)
[qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(ECG_signal,f_sample_ECG,0);

% calculate heart rate. BPM_ECG = beats per minute. IBI_ECG = inter beat
% interval
[BPM_ECG, IBI_ECG] = Calculate_rate_ECG(qrs_i_raw, delta_tECG);

% use windowing function. (The reason to why we use windowing function for the ECG signal is that
% we do it for the radar signal to avoid false detections of peaks. When we calculate the accuracy and plot the heart rate as a
% function of time for both the radar and the ECG it is more clear if a windowing function is applied on both signals and not just
% the radar signal)
[evolution_ECG] = WindowingECG_HR(qrs_i_raw,delta_tECG,length(ECG_signal));

t_ECG=linspace(0,tECG(end),length(BPM_ECG)); % time array for the ECG signal
t_ECG_ev = linspace(7.5,tECG(end)-7.5,length(evolution_ECG)); % time array for evolution signal

% Plot BPM and IBI for ECG as functions of time
figure(2)
plot(t_ECG,BPM_ECG, 'LineWidth',0.9)
ylabel('Hjärtfrekvens [BPM]', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTex')
xlabel('Tid [s]', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTeX')

figure(3)
plot(t_ECG,IBI_ECG)
ylabel('IBI [s]', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTex')
xlabel('Time [s]', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTeX')

figure(4)
plot(t_ECG_ev,evolution_ECG)

%% filter out the heart rate, broad bandpassfilter
t=0:delta_t:(N_frame-1)*delta_t; % time array
filter_06_35 = fir1(1500, [0.6/fn,3.5/fn], 'bandpass'); % Broad bandpass filter. The order of the filter and the frequencies may be changed if wanted.
% Note: the order of the filter must be even.
HR_06_35 = movmean(filtfilt(filter_06_35,1,phase_unwraped),5); % filtered HR

%% filter out the heart rate, narrow bandpass
b = fir1(1500, [1/fn,1.5/fn]); % narrow bandpass filter.The order of the filter and the frequencies may be changed if wanted.
% Note: the order of the filter must be even.
HR_narrow = movmean((filtfilt(b,1, phase_unwraped)),7); % filtered HR

%% filter out the heart rate by looking at the HR from the ECG signal, apply different bandpass filters on different parts of the signal

number_of_samples_filter = 400; % length of each interval is 400 samples = 16 seconds, to be modified if wanted different
time_filter = number_of_samples_filter*delta_t; % length of an interval in seconds
ECG_frequency = 1./IBI_ECG; % frequancy of heart beats (in Hz)
HR_manually_filtered = zeros(length(phase_unwraped),1);

end_time = time_filter; % the time of the end point of each interval. For the first interval end_time = 16,
% for the second interval end_time = 32.... (if time_filter = 16)
start_time = 0; % the time of the start point of the interval. For the first interval start_time =0,
% for the second interval start_time =16....
stop = 0;
j=1;

while  stop==0
    % find max and min frequency for each interval:
    if end_time<t_ECG(end)
        index_of_first = find(t_ECG >= start_time, 1, 'first');
        index_of_last  = find(t_ECG <= end_time, 1, 'last');
        max_value = max(ECG_frequency(index_of_first:index_of_last));
        min_value = -max(-ECG_frequency(index_of_first:index_of_last));
    else
        index_of_first = find(t_ECG >= start_time, 1, 'first');
        max_value = max(ECG_frequency(index_of_first:end));
        min_value = -max(-ECG_frequency(index_of_first:end));
        stop = 1;
    end

    %Update start_time and end_time:
    start_time = end_time;
    end_time=end_time + time_filter;

% choose limits for bandpass filter, to be modified if wanted different
    filter_low_limit_ECG = min_value-0.5;
    filter_high_limit_ECG = max_value+0.5;
     if filter_low_limit_ECG < 0.6
         filter_low_limit_ECG = 0.6;
     end
%      % use a threashold if wanted:
%     if filter_high_limit_ECG < 1.3
%         filter_high_limit_ECG = 1.3;
%     end

% apply bandpass filter (the order may be changed if wanted):
    filter_HR_ECG = fir1(120, [filter_low_limit_ECG/fn,filter_high_limit_ECG/fn], 'bandpass');
    HR_manually_filtered(j:j+number_of_samples_filter-1) = filtfilt(filter_HR_ECG,1, phase_unwraped(j:j+number_of_samples_filter-1));
    %HR_manually_filtered(j:j+number_of_samples_filter-1) = filtfilt(filter_HR_ECG,1, unwrap(phase(j:j+number_of_samples_filter-1,10)));

    j = j+number_of_samples_filter;
end


%% Calculation of evolutions and plot heart rate as a function of time (manuel filter)

[evolutionHR_f, index_of_peaks, values_peaks]=WindowingHR(movmean(HR_manually_filtered,5),delta_t,N_frame); % calculate HR evolution from radar signal
radar_peak_time = index_of_peaks.*delta_t; % the time (in seconds) that corresponds to the detected peaks in the unwraped phase signal

t_HR=linspace(7.5,t(end)-7.5,length(evolutionHR_f)); % time array corresponding to heart rate evolution from radar (if each window is 15 seconds)
t_ECG=linspace(7.5,tECG(end)-7.5,length(BPM_ECG)); % time array corresponding to heart rate evolution from ECG (if each window is 15 seconds)

% Plot BPM and IBI for ECG as functions of time
figure(5)
plot(t_ECG,BPM_ECG)
hold off
figure(6)
plot(t_ECG,IBI_ECG)

% calculcation of accurancy
[accuracy_mean, accuracy_median]=accuracy(evolutionHR_f,evolution_ECG)

% plots
figure(9);
plot(t_HR,evolutionHR_f, 'r-x')
hold on
%plot(t_HR,evolutionHR_b, 'r-x')
hold on
%plot(t_ECG,BPM_ECG, 'k--')
hold on
plot(t_ECG_ev,evolution_ECG, 'k')
hold on
%plot(t_HR,evolutionHR_f, 'rx')
hold on
%plot(t_HR,evolutionHR_6, 'm')
hold off
ylabel('Hjärtfrekvens [BPM]', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTex')
xlabel('Tid [s]', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTeX')
%title('\fontsize{18} Hög puls för försöksperson A')

lgd = legend({'Radar','Referens'},'FontSize',16,'TextColor','black');
lgd.NumColumns = 1;


%%
% the code below was used to filter the phase signal automatically,
% without any information from the reference signal. However, many
% assumptions were made (see the bachelor thesis of 2022 to read about this more in detail) and it is not
% sure that this code will work for all cases. It is also recommended to
% use a function to do this and not have the entire code in the main script
% :)


L = 2500; % length of the zero-padded signal (change if you want longer och shorter signal)
number_of_samples = 400; % number of samples of the phase signal that we want to analyse in frequency domain
number_of_zeros = L-400; % number of zeros that have to be added
f = fsample.*(0:(L/2))/L;
q = fir1(1000, [0.1/fn,3.5/fn], 'bandpass'); % bandpass filter
f_vec = zeros(6000/number_of_samples,1); % vector that will contain the estimated frecuency from each interval
m=1; % counts iterations
pm_filter = 0.11; % set bandstop filter limits for removing heartbeat harmonics
pm_filter_b = 0.11; % set bandstop filter limits for removing breathing harmonics
HR = zeros(length(phase_unwraped),1); % vector that will contain the filtered HR-signal (after estimation of frecuency)
hb = filtfilt(q,1, phase_unwraped); % HR-signal 0,1 Hz to 3,5 Hz

% for-loop iterates through the phase signal, 400 samples each time:
for i = 1:number_of_samples:6000
heart = movmean(hb(i:i+number_of_samples-1),5); % moving average to remove high frecuency components from noise and random movements
hw = hann(length(heart)); % hanning window
r = [hw.*heart;zeros(number_of_zeros,1)]; % zero padded signal

Y = fft(r); % applt fft
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1); % one-sided amplitude spectrum
%figure(20)
%plot(f,P1)

f_breathing = f(17:75); % frequency interval for breathing
amp_breathing = P1(17:75); % amplitude of the frequency interval for breathing
[breathing_values, breathing_peaks] = findpeaks(amp_breathing); % find possible breathing peaks

% if we only have one possible breathing peak then that frequency will be
% considered the breathing frequency, f_b:
if length(breathing_peaks) == 1
[b_max_value, b_max_index] = maxk(breathing_values,1);
f_b = f_breathing(breathing_peaks(b_max_index(1))); % breathing frequency
else

% If we have more than one possible breathing peaks:
[b_max_value, b_max_index] = maxk(breathing_values,2);
f_b_max_1 = f_breathing(breathing_peaks(b_max_index(1))); % frequency with highest amplitude
f_b_max_2 = f_breathing(breathing_peaks(b_max_index(2))); % frequency with next highest amplitude
possible_harmonic_b = 0;

% if f_b_max_2 is the second harmonic of f_b_max_1 possible_harmonic_b = 1, if not possible_harmonic_b = 0
if (2*f_b_max_1 < f_b_max_2 + 0.07) && (2*f_b_max_1 > f_b_max_2 - 0.07)
    possible_harmonic_b = 1;
end

% check if we want to look at the peak with highest amplitude or the peak
% with next highest amplitude.
if f_b_max_1 < 0.25 && 0.65*b_max_value(1) < b_max_value(2) && possible_harmonic_b == 0
   f_b = f_b_max_2;
   breathing_max_value = b_max_value(1);
else
   f_b = f_b_max_1;
   breathing_max_value = b_max_value(2);
end
end

% chose frequncy interval for heart rate:
% Note: if you change the number of sample that you apply fft on,
% you should also change the values of h_first
if f_b < 0.35
    h_first = 76;
elseif f_b < 0.39
    h_first = 90;
elseif f_b < 0.45
    h_first = 96;
elseif f_b < 0.55
    h_first = 106;
else
    h_first = 113;
end

% find the highest peak in the frequency interval for the heart rate
[heart_peaks_value, heart_peaks_index] = findpeaks(P1(h_first:501));
[heart_max_peak, heart_max_peak_index] = maxk(heart_peaks_value,6); % heart_max_peak = "amplitude of highest peak in the frequency interval for the heart rate"
f_max_heart = f(h_first+heart_peaks_index(heart_max_peak_index)-1); % frequency of the 6 peaks with highest amplitude
possible_harmonic_heart = 0;

% check if one of the 6 peaks with highest amplitudes is a harmonic of the
% frequency of the highest peak
for k = 2:length(f_max_heart)
    if (2*f_max_heart(1) < f_max_heart(k) + 0.07) && (2*f_max_heart(1) > f_max_heart(k) - 0.07)
        possible_harmonic_heart = 1;
        break
    end
end

% remove the second harmonic from breathing if the heart rate is not
% assmued to have the same frequency
if (2*f_b < f_max_heart(1) + pm_filter_b+0.01) && (2*f_b > (f_max_heart(1) - pm_filter_b-0.01)) && (heart_max_peak(1) > 0.8*breathing_max_value) && possible_harmonic_heart == 1
    hb_1 = heart;
else
    low_limit = 2*f_b-0.11;
    high_limit = 2*f_b+0.11;
    hb_1 = bandstop(heart,[low_limit, high_limit],25);
end


% remove the third harmonic from breathing if the heart rate is not
% assmued to have the same frequency
if  (3*f_b < f_max_heart(1) + pm_filter_b+0.01) && (3*f_b > (f_max_heart(1) - pm_filter_b-0.01)) && (heart_max_peak(1) > 0.8*breathing_max_value) && possible_harmonic_heart == 1
    hb_2 = hb_1; % does not remove third harmonic
elseif f_b > 0.45
    f_b_harmonic_3 = 3*f_b;
    low_limit_3 = f_b_harmonic_3-pm_filter_b;
    high_limit_3 = f_b_harmonic_3+pm_filter_b;
    hb_2 = bandstop(hb_1,[low_limit_3, high_limit_3],25); % removes third harmonic
    yy=2;
else
    hb_2 = hb_1;
    yy=3;
end

% apply fft on the filtered signal:
hw_2 = hann(length(hb_2));
r_2 = [hw_2.*hb_2;zeros(number_of_zeros,1)];
Y2 = fft(r_2);
P2_2 = abs(Y2/L);
P1_2 = P2_2(1:L/2+1);
P1_2(2:end-1) = 2*P1_2(2:end-1); % new single sided spectrum
%figure(22)
%plot(f,P1_2)

% to remove harmonics we first divide the heart frequency interval in 3 parts.
% Below is the first part of the interval:
if f_b < 0.35
    f_limit = 0.75;
    f_h1 = f(76:151);
    a_h1 = P1_2(76:151);
    first = 76;
elseif f_b < 0.39
    f_limit = 0.9;
    f_h1 = f(90:151);
    a_h1 = P1_2(90:151);
    first = 90;
elseif f_b < 0.45
    f_limit = 0.95;
    f_h1 = f(96:151);
    a_h1 = P1_2(96:151);
    first = 96;
elseif f_b < 0.55
    f_limit = 1.05;
    f_h1 = f(106:151);
    a_h1 = P1_2(106:151);
    first = 106;
else
    f_limit = 1.12;
    f_h1 = f(113:151);
    a_h1 = P1_2(113:151);
    first = 113;
end


[max_heart, index_max_heart] = maxk([a_h1;P1_2(152:501)],1); % Find the highest peak in the entire frequency interval for the heart rate. Note: since we have applied
% a bandpass filter there is no need to look at the the last samples,
% therefore we only consider at the first 501 samples. (You can look at even fewer samples if wanted)
[peaks_value_heart_1, peak_index_heart_1] = findpeaks(a_h1); % find the peaks in the first part of the interval
[max_heart_1, index_heart_1] = max(peaks_value_heart_1); % highest peak in the first part of the interval

% if the highest peak in the first part of the heart rate interval is
% higher than 50% of the amplitude of the highest peak then that frequency
% is considered a potential heart frequency and the second harmonic is
% removed
if max_heart_1 > 0.5*max_heart
    f_max_1 = f_h1(peak_index_heart_1(index_heart_1));
    hb_3 = bandstop(hb_2, [(2*f_max_1-pm_filter),(2*f_max_1+pm_filter)],25);
    h=1;
else
    hb_3 = hb_2;
    h=0;
end

% apply fft on the filtered signal:
hw_3 = hann(length(hb_3));
r_3 = [hw_3.*hb_3;zeros(number_of_zeros,1)];
Y3 = fft(r_3);
P2_3 = abs(Y3/L);
P1_3 = P2_3(1:L/2+1);
P1_3(2:end-1) = 2*P1_3(2:end-1); % new single sided spectrum
%figure(23)
%plot(f,P1_3)

% the second part of the heart rate interval
f_h2 = f(150:300);
a_h2 = P1_3(150:300);

[max_heart_2, index_max_heart_2] = maxk(a_h2,1);

% remove second harmonic as above
if max_heart_2 > 0.5*max_heart
    f_max_2 = f_h2(index_max_heart_2);
    hb_4 = bandstop(hb_3, [(2*f_max_2-pm_filter),(2*f_max_2+pm_filter)],25);
else
    hb_4 = hb_3;
end

% apply fft on the filtered signal:
hw_4 = hann(length(hb_4));
r_4 = [hw_4.*hb_4;zeros(number_of_zeros,1)];
Y4 = fft(r_4);
P2_4 = abs(Y4/L);
P1_4 = P2_4(1:L/2+1);
P1_4(2:end-1) = 2*P1_4(2:end-1); % new single sided spectrum
%figure(24)
%plot(f,P1_4)

% choose the highest (remaining) peak as the heart frequency:
[peak_value,peak_index] = findpeaks(P1_4(first:501));
[max_value, index_value] = maxk(peak_value,1);
freq = f(first:501);
f_vec(m) = freq(peak_index(index_value));
m=m+1;
end

b=1;
frequency = zeros(6000/number_of_samples,1);

% apply bandpass filters:

% lower limit of filter is 0,3 Hz lower than the estimated frequency for
% each interval, higher limit of filter is 0,3 Hz higher than the estimated
% frequency for each interval. These limits may be modified if wanted.
for c = 1:length(f_vec)
    if c == 1 || c == length(f_vec)
        frequency(c) = f_vec(c);
        filter_low_limit = f_vec(c)-0.3;
        filter_high_limit = f_vec(c)+0.3;
        if filter_low_limit < 0.6
            filter_low_limit = 0.6;
        end
        if filter_high_limit < 1.2
            filter_high_limit = 1.2;
        end
        filter_HR = fir1(120, [filter_low_limit/fn,filter_high_limit/fn], 'bandpass');
        HR(b:b+number_of_samples-1) = filtfilt(filter_HR,1, phase_unwraped(b:b+number_of_samples-1));

    else
        if 60*f_vec(c-1) + 40 < 60*f_vec(c) && 60*f_vec(c+1)+40 < 60*f_vec(c)
            frequency(c) = mean([f_vec(c-1),f_vec(c+1)]);
        elseif 60*f_vec(c-1) - 40 > 60*f_vec(c) && 60*f_vec(c+1)-40 > 60*f_vec(c)
            frequency(c) = mean([f_vec(c-1),f_vec(c+1)]);
        else
            frequency(c) = f_vec(c);
        end

        filter_low_limit = frequency(c)-0.3;
        filter_high_limit = frequency(c)+0.3;
        if filter_low_limit < 0.6
            filter_low_limit = 0.6;
        end
        if filter_high_limit < 1.2
            filter_high_limit = 1.2;
        end
        filter_HR = fir1(120, [filter_low_limit/fn,filter_high_limit/fn], 'bandpass');
        HR(b:b+number_of_samples-1) = filtfilt(filter_HR,1, phase_unwraped(b:b+number_of_samples-1));
    end
   b = b+400;
end


%% calculate evolutions and plot heart rate as a function of time (automatic filter)

[evolutionHR, HR_peak_index, HR_peak_values]=WindowingHR(movmean(HR,5),delta_t,N_frame);
[rate_evolution] = WindowingECG_HR(qrs_i_raw,delta_tECG,length(ECG_signal));
t_HR=linspace(7.5,t(end)-7.5,length(evolutionHR));
t_ECG=linspace(0,tECG(end),length(BPM_ECG));
t_ECG_ev = linspace(7.5,tECG(end)-7.5,length(rate_evolution));

figure(30);
plot(t_HR,evolutionHR, 'r-x')
hold on
plot(t_ECG_ev,rate_evolution, 'k')
hold off

ylabel('Hjärtfrekvens [BPM]', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTex')
xlabel('Tid [s]', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTeX')
%title('\fontsize{18} Hög puls för försöksperson A')

lgd = legend({'Radar','Referens'},'FontSize',16,'TextColor','black');
lgd.NumColumns = 1;

BPM_1 = 60.*f_vec;
BPM_2 = 60.*frequency;

figure(31)
%plot(t_ECG,BPM_ECG)
%hold on
%plot(linspace(8,232,length(BPM_1)),BPM_1, 'b-x', 'LineWidth',1)
hold on
plot(linspace(8,232,length(BPM_2)),BPM_2, 'b-x', 'LineWidth',1)
hold on
plot(t_ECG_ev,rate_evolution, 'k', 'LineWidth',1)
ylabel('Hjärtfrekvens [BPM]', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTex')
xlabel('Tid [s]', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTeX')
lgd = legend({'Uppskattad frekvens','Referens'},'FontSize',16,'TextColor','black');
lgd.NumColumns = 1;
%axis([0 250 50 120])

[accuracy_mean, accuracy_median]=accuracy(evolutionHR,evolution_ECG)

%%
% the code below was used to plot the frequency spectrums before and after
% the harmonics were removed for different part of the phase signal. The
% code is very similar to the code above, hence ýou can only find a few
% comments in the code below.

clf
L = 2500;
number_of_zeros = 2500-400;
%L=300;
f = fsample.*(0:(L/2))/L;
q = fir1(1000, [0.1/fn,3.5/fn], 'bandpass');
number_of_samples = 400;
f_vec = zeros(6000/number_of_samples,1);
m=1;
pm_filter = 0.11;
pm_filter_b = 0.11;

hb = filtfilt(q,1, phase_unwraped);
heart = movmean(hb(1:400),1); % choose interval here!!
hw = hann(length(heart));
r = [hw.*heart;zeros(number_of_zeros,1)];

Y = fft(r);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure(20)
plot(f,P1./(max(P1)), 'LineWidth',1)
ylabel('Amplitud (normaliserad)', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTex')
xlabel('Frekvens [Hz]', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTeX')
xlim([0 4])
ylim([0 1.12])

% 0.12 till 0.7568
f_breathing = f(17:75);
amp_breathing = P1(17:75);
%[max_b, index_b] = maxk(amp_breathing,1);
%[breathing_values, breathing_peaks] = findpeaks(amp_breathing, 'MinPeakHeight', 0.65*max_b);
[breathing_values, breathing_peaks] = findpeaks(amp_breathing);

% if length(breathing_peaks) == 0
%     [breathing_values, breathing_peaks] = findpeaks(amp_breathing);
% end
if length(breathing_peaks) == 1
[b_max_value, b_max_index] = maxk(breathing_values,1);
f_b = f_breathing(breathing_peaks(b_max_index(1)));
else
[b_max_value, b_max_index] = maxk(breathing_values,2);
f_b_max_1 = f_breathing(breathing_peaks(b_max_index(1)));
f_b_max_2 = f_breathing(breathing_peaks(b_max_index(2)));
possible_harmonic_b = 0;

if (2*f_b_max_1 < f_b_max_2 + 0.07) && (2*f_b_max_1 > f_b_max_2 - 0.07)
    possible_harmonic_b = 1;
end

if f_b_max_1 < 0.25 && 0.65*b_max_value(1) < b_max_value(2) && possible_harmonic_b == 0
   f_b = f_b_max_2;
   breathing_max_value = b_max_value(1);
else
   f_b = f_b_max_1;
   breathing_max_value = b_max_value(2);
end
end


if f_b < 0.35
    h_first = 76;
elseif f_b < 0.39 %ändrat
    h_first = 90;
elseif f_b < 0.45
    h_first = 96;
elseif f_b < 0.55
    h_first = 106;
else
    h_first = 113;
end

[heart_peaks_value, heart_peaks_index] = findpeaks(P1(h_first:501));
[heart_max_peak, heart_mx_peak_index] = maxk(heart_peaks_value,6);
f_max_heart = f(h_first+heart_peaks_index(heart_mx_peak_index)-1);
possible_harmonic_heart = 0;

for k = 2:length(f_max_heart)
    if (2*f_max_heart(1) < f_max_heart(k) + 0.07) && (2*f_max_heart(1) > f_max_heart(k) - 0.07)
        possible_harmonic_heart = 1;
        break
    end
end

if (2*f_b < f_max_heart(1) + pm_filter_b+0.01) && (2*f_b > (f_max_heart(1) - pm_filter_b-0.01)) && (heart_max_peak(1) > 0.8*breathing_max_value) && possible_harmonic_heart == 1
    hb_1 = heart;
else
    low_limit = 2*f_b-0.11;
    high_limit = 2*f_b+0.11; % kan behöva ändras
    hb_1 = bandstop(heart,[low_limit, high_limit],25);
end


if  (3*f_b < f_max_heart(1) + pm_filter+0.01) && (3*f_b > (f_max_heart(1) - pm_filter-0.01)) && (heart_max_peak(1) > 0.8*breathing_max_value) && possible_harmonic_heart == 1
    hb_2 = hb_1;
elseif f_b > 0.45
    f_b_harmonic_3 = 3*f_b;
    low_limit_3 = f_b_harmonic_3-pm_filter;
    high_limit_3 = f_b_harmonic_3+pm_filter;
    hb_2 = bandstop(hb_1,[low_limit_3, high_limit_3],25);
    yy=2;
else
    hb_2 = hb_1;
    yy=3;
end


hw_2 = hann(length(hb_2));
r_2 = [hw_2.*hb_2;zeros(number_of_zeros,1)];
Y2 = fft(r_2);
P2_2 = abs(Y2/L);
P1_2 = P2_2(1:L/2+1);
P1_2(2:end-1) = 2*P1_2(2:end-1);
%figure(22)
%plot(f,P1_2)

if f_b < 0.35
    f_limit = 0.75;
    % hjärta 0,75 till 1.5 Hz
    f_h1 = f(76:151);
    a_h1 = P1_2(76:151);
    first = 76;
elseif f_b < 0.39
    f_limit = 0.9;
    % hjärta 0.8911 till 1.4893 Hz
    f_h1 = f(90:151);
    a_h1 = P1_2(90:151);
    first = 90;
elseif f_b < 0.45
    f_limit = 0.95;
    % hjärta 0.9888 till 1.4893 Hz
    f_h1 = f(96:151);
    a_h1 = P1_2(96:151);
    first = 96;
elseif f_b < 0.55
    f_limit = 1.05;
    % hjärta 1.0376 till 1.4893 Hz
    f_h1 = f(106:151);
    a_h1 = P1_2(106:151);
    first = 106;
else
    f_limit = 1.12;
    % hjärta 1.1475 till 1.4893 Hz
    f_h1 = f(113:151);
    a_h1 = P1_2(113:151);
    first = 113;
end

[max_heart, index_max_heart] = maxk([a_h1;P1_2(152:501)],1);
[peaks_value_heart_1, peak_index_heart_1] = findpeaks(a_h1);


[max_heart_1, index_heart_1] = max(peaks_value_heart_1);
if max_heart_1 > 0.5*max_heart
    f_max_1 = f_h1(peak_index_heart_1(index_heart_1));
    hb_3 = bandstop(hb_2, [(2*f_max_1-pm_filter),(2*f_max_1+pm_filter)],25); % kan behöva öka dessa
    h=1;
else
    hb_3 = hb_2;
    h=0;
end

hw_3 = hann(length(hb_3));
r_3 = [hw_3.*hb_3;zeros(number_of_zeros,1)];
Y3 = fft(r_3);
P2_3 = abs(Y3/L);
P1_3 = P2_3(1:L/2+1);
P1_3(2:end-1) = 2*P1_3(2:end-1);
%figure(23)
%plot(f,P1_3)

% hjärta 1.49 till 2.99 Hz
f_h2 = f(150:300);
a_h2 = P1_3(150:300);

[max_heart_2, index_max_heart_2] = maxk(a_h2,1);

if max_heart_2 > 0.5*max_heart
    f_max_2 = f_h2(index_max_heart_2);
    %s = fir1(98, [(2*f_max_2-0.06)/fn,(2*f_max_2+0.06)/fn], 'stop');
    %hb_4 = filtfilt(s,1, hb_3);
    hb_4 = bandstop(hb_3, [(2*f_max_2-pm_filter),(2*f_max_2+pm_filter)],25);
else
    hb_4 = hb_3;
end

hw_4 = hann(length(hb_4));
r_4 = [hw_4.*hb_4;zeros(number_of_zeros,1)];
Y4 = fft(r_4);
P2_4 = abs(Y4/L);
P1_4 = P2_4(1:L/2+1);
P1_4(2:end-1) = 2*P1_4(2:end-1);
figure(24)
plot(f,P1_4./max(P1_4), 'LineWidth',1)
ylabel('Amplitud (normaliserad)', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTex')
xlabel('Frekvens [Hz]', 'FontName', 'Times New Roman','FontSize',20,'Color','k', 'Interpreter', 'LaTeX')
xlim([0 4])
ylim([0 1.12])
hold on
xline(f_limit, 'LineWidth',2)
xline(3.5, 'LineWidth',2)

%
[peak_value,peak_index] = findpeaks(P1_4(first:501));
[max_value, index_value] = maxk(peak_value,1);
freq = f(first:501);
f_vec(m) = freq(peak_index(index_value));
m=m+1;




