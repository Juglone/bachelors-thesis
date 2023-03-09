function [ECG_signal, sample_nr_signal] = synchronize_ECG_Radar(ECG,sample_nr,delta_t_ECG,delta_t_radar,N_frames)
% Use this function to synchronize ECG data with data fram radar.
% Note: this can only be used if ECG measurement ended at the same time as the
% radar measurement ended and if the ECG measurement was started before the
% radar measurement

radar_time = delta_t_radar*N_frames;
ECG_time = delta_t_ECG*(length(ECG));

time_diff = ECG_time-radar_time;

if time_diff >= 0
samples_to_remove = round(time_diff/delta_t_ECG);
ECG_signal = ECG(samples_to_remove+1:end);
sample_nr_signal = sample_nr(samples_to_remove+1:end)-samples_to_remove;
else 
 ECG_signal = ECG(1:end);
sample_nr_signal = sample_nr(1:end);   
end
end
