% == Find range bin of each chirp ==
% This function will detect the range bin for each chirp. The bin is found
% as the peak of highest magnitude within each chirp. Because there is
% always a peak from the tx antenna the function skips the first few
% data points of each chirp.

function [bin_magnitudes, bin_indices] = find_range_bins(num_frames, fft_magnitude)
    % Create empty lists to be filled by the for loop
    bin_magnitudes = zeros(num_frames,1);
    bin_indices = zeros(num_frames,1);
    
    % Loop through the chirps and find the range bin of each chirp
    for i = 1:num_frames
        Y = fft_magnitude(i,:);
        [M_i,I_i] = max(Y(5:end));  % Skip the peak of tx antenna
        I_i = I_i + 4;              % Adjust for the index skip
        bin_magnitudes(i) = M_i;     % Save the magnitude of each chirp bin  
        bin_indices(i) = I_i;       % Save the index of each chirp bin
    end
end