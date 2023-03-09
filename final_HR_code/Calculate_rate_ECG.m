function [BPM, IBI] = Calculate_rate_ECG(index_RR,delta_t)
% calculate heart rate

IBI = diff(delta_t.*index_RR);
BPM = 60./IBI;

end 