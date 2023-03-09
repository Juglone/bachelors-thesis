function [phase_unwraped, finalRangeBins, rangebins, m, selected_rangeBins, L] = selectRangeBin_4(Magnitude, Phase, N_frame)
%phase_unwraped = unwraped phase signal
% finalRangeBins = the range bins that we think are possible range bins.
% After finding these range bins we look at the variance.
% m = number of peaks for the four range bins with highest variance
% selected_rangeBins = the selected range bins (with this function you always select 1-4 range bins). 
% L = number of selected range bins
% (could be suitible to change the names of these variables...)

% find the peaks of the magnitude of one (arbitrary) chirp:
[amplitudeValue,peakIndex]=findpeaks(Magnitude(round(N_frame/2),:)); 

maxPeak=max(amplitudeValue); % the value of the highest peak

% neglect the peaks that have a lower amplitude than 20% of the max value:
threshold=0.2*maxPeak;
RangeBins_threshold = zeros(length(amplitudeValue),1);
for j=1:length(amplitudeValue)
    if amplitudeValue(j) > threshold
        RangeBins_threshold(j,1) = peakIndex(j);
    else
        RangeBins_threshold(j,1) = 0;
    end
end

rangeBinsFromPeaks = nonzeros(RangeBins_threshold);

% possible range bins are the range bins that corresponds to peaks and the
% four closest range bins to each peak. In the for-loop below we find possible range bins.
% Note that possibleRangeBins(n) >=1.

possibleRangeBins = zeros(N_frame,1);
n=1;
for i = 1:length(rangeBinsFromPeaks)
    if rangeBinsFromPeaks(i) - 2 >= 1
        possibleRangeBins(n) = rangeBinsFromPeaks(i)-2;
        possibleRangeBins(n+1) = rangeBinsFromPeaks(i)-1;
        possibleRangeBins(n+2) = rangeBinsFromPeaks(i);
        possibleRangeBins(n+3) = rangeBinsFromPeaks(i)+1;
        possibleRangeBins(n+4) = rangeBinsFromPeaks(i)+2;
        n = n+5;
    elseif rangeBinsFromPeaks(i) - 1 >= 1
        possibleRangeBins(n) = rangeBinsFromPeaks(i)-1;
        possibleRangeBins(n+1) = rangeBinsFromPeaks(i);
        possibleRangeBins(n+2) = rangeBinsFromPeaks(i)+1;
        possibleRangeBins(n+3) = rangeBinsFromPeaks(i)+2;
        n = n+4;
    else 
        possibleRangeBins(n) = rangeBinsFromPeaks(i);
        possibleRangeBins(n+1) = rangeBinsFromPeaks(i)+1;
        possibleRangeBins(n+2) = rangeBinsFromPeaks(i)+2;  
        n = n+3;
    end
end

finalRangeBins = unique(nonzeros(possibleRangeBins));
phase_variance = zeros(length(finalRangeBins),1);

fn = 12.5; % nyquist, change if sampling frequency changes
b = fir1(300, [0.6/fn,3.5/fn]); % bandpass filter, to be modified if wanted different

for k=1:length(finalRangeBins)
    phase_variance(k,1) = var(filtfilt(b,1,unwrap(Phase(:,finalRangeBins(k))))); % calculate the variance
end

% choose the 4 range bins with highest variance:
[variance,index] = maxk(phase_variance,4); 
rangebins = finalRangeBins(index); % rangebins = the four range bins with highest variance

% look at the peaks of each chirp. Calculate how many chirps that have
% peaks corresponding to each of the four range bins with highest variance:
m = zeros(4,1); % m(i) = number of chirps that have a peak for rangebins(i) 
for a = 1:N_frame
    [amplitude, index_range] = findpeaks(Magnitude(a,:));
    if any(ismember(rangebins(1), index_range))
        m(1,1) = m(1,1)+1;
    end
    if any(ismember(rangebins(2), index_range))
       m(2,1) = m(2,1)+1;
    end
     if any(ismember(rangebins(3), index_range))
       m(3,1) = m(3,1)+1;
     end
     if any(ismember(rangebins(4), index_range))
       m(4,1) = m(4,1)+1;
     end
%      if ismember(rangebins(5), index_range)
%        m(5,1) = m(5,1)+1;
%      end
end 

% choose the range bins that have at least 30 % of the maximum number of
% peaks
[max_value, index_of_max] = max(m);
rangebins_with_peaks = zeros(length(m),1);
for s = 1:length(m)
if m(s,1) >= 0.3*max_value
    rangebins_with_peaks(s,1) = rangebins(s,1);
end 
end 
selected_rangeBins = nonzeros(rangebins_with_peaks); 

% uwrap the phase signal:
if length(selected_rangeBins) == 1
phase_unwraped = unwrap(Phase(:,selected_rangeBins(1,1)));
L=1;

elseif length(selected_rangeBins) == 2
phase_unwraped = (unwrap(Phase(:,selected_rangeBins(1,1))) + unwrap(Phase(:,selected_rangeBins(2,1))))/2;
L=2;

elseif length(selected_rangeBins) == 3
phase_unwraped = (unwrap(Phase(:,selected_rangeBins(1,1)))+ unwrap(Phase(:,selected_rangeBins(2,1)))+ unwrap(Phase(:,selected_rangeBins(3,1))))/3;
L=3;

elseif length(selected_rangeBins) == 4
phase_unwraped = (unwrap(Phase(:,selected_rangeBins(1,1))) + unwrap(Phase(:,selected_rangeBins(2,1))) + unwrap(Phase(:,selected_rangeBins(3,1))) + unwrap(Phase(:,selected_rangeBins(4,1))))/4;
L=4;
end

end