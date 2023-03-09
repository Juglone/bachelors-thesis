function [heart_rate] = WindowingECG_HR(index_RR,delta_t,Nb_frame) 
% Note: index_RR = index of R-peaks, delta_t = time difference between two samples, 
% Nb_frame = number of samples of the ECG signal

% Note: a window will be overlaping the previous one of half it's length

size_w_time=15; %length of each window in seconds, to be modified if wanted different 
rate_evolution=zeros(length(index_RR),1); % will contain the values of HR for each window

end_w_time=size_w_time; % the total time of the last point in the window. 
% For example, for the first window: end_w_time = 15, for the second window
% end_w_time = 15+7.5 = 22.5, for the third window end_w_time = 30 ...
i=1;
q=0;
index = zeros(length(index_RR),2); % Note: this matrix is never used and can be remoeved if wanted

while  q==0
    if end_w_time<delta_t*index_RR(end) % check if we have at least one R-pak after the end point of the window
        % If we do have at least one R-peak after the end point of the
        % window:
        index_of_first = find(index_RR >= (end_w_time - size_w_time)/delta_t, 1, 'first'); % index of first R-peak in the window
        index_of_last  = find(index_RR <= end_w_time/delta_t, 1, 'last'); % index of last R-peak in the window
        window=index_RR(index_of_first:index_of_last); % the indices of all R-peaks in the window
        index(i,1) = index_RR(index_of_first); 
        index(i,2) = index_RR(index_of_last); 
    else
        % if there is no R-peaks after the end point of the window:
        index_of_first = find(index_RR >= (end_w_time - size_w_time)/delta_t, 1, 'first'); % index of first R-peak in the window
        window=index_RR(index_of_first:end); % the indices of all R-peaks in the window
        index(i,1) = index_RR(index_of_first);
        index(i,2) = index_RR(end);
        q=1;
    end 
   
    [BPM, IBI]=Calculate_rate_ECG(window,delta_t); % calculate BPM and IBI
    rate_evolution(i,1)=mean(BPM); % rate_evolution = mean heart rate in window
   
    end_w_time=size_w_time/2 + end_w_time;
    i=i+1;
  
end

heart_rate = nonzeros(rate_evolution);
end 