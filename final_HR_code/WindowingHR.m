function [rate_evolution, index_of_peaks, values_peaks] = WindowingHR(filtered_HR_signal,delta_t,Nb_frame)
% Note: the majority of this code is from the master thesis 2021

% a window will be overlaping the previous one of half it's length

size_w=15; %seconds, in time domain, to be modified if wanted different

size_w=size_w/delta_t; % number of samples in one window, note: not rounded!
nb_w= floor(Nb_frame/(size_w/2)-2); %total number of windows, each window overlaps the other one by half, floor(X) rounds the elements of X to the nearest integers towards minus infinity.
size_w=floor(size_w); % number of samples in one window, note: floor(X) rounds the elements of X to the nearest integers towards minus infinity.
min_dist=0.15/delta_t; % min time between two hert beats is 0,15 seconds, to be modified if wanted different
half_w=round(size_w/2); % length of half a window 

rate_evolution=zeros(nb_w,1); %will keep in mind the values of heart rate for each window
index_of_peaks_zeros = zeros(length(filtered_HR_signal),1); %will keep in mind the indices of the peaks (in the unwraped phase signal)

k=1; % start point of window
m=1; % index for the array index_of_peaks_zeros
end_w=size_w; % end point of window
for i=0:nb_w 
    if end_w<Nb_frame 
        window=filtered_HR_signal(k:end_w);
    else
        window=filtered_HR_signal(k:Nb_frame);
    end 
    adapt_prom=0.1*rms(window); % adaptable prominence is 10% of rms-value of each window
    [value,index]=findpeaks(window,'MinPeakProminence',adapt_prom,'MinPeakDistance',min_dist); %Find peaks of the unwraped phase signal of one window. 
    % Note 'MinPeakDistance' is set to 0.15s

    index_of_peaks_zeros(m:(m+length(index)-1)) = index+k-1;

    [BPM, IBI]=Calculate_rate_ECG(index,delta_t); % calculate BPM
    rate_evolution(i+1)=mean(BPM); % rate_evolution = mean heart rate in window
    k=k+half_w; %start for the next window
    end_w=k+size_w; % end of the next window
    m = m+length(index);
end
index_of_peaks = unique(nonzeros(index_of_peaks_zeros)); % only keep one of each index
values_peaks = filtered_HR_signal(index_of_peaks); % the values of the phase signal that correspond to the indices
end 