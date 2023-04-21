%% Remove peaks with wavelet and RMS hard thresholding

function [despiked_signal, spikes] = wavelet_leveldepend(spiked_signal, wtype, wlevel, thresh, movw, plot_flag)
%t = t_radar;
    % Signal to de-spike
    y = spiked_signal;

% Set wavelet type and level
%wtype = 'wtype';
% = 5;

    % Perform a wavelet decomposition of the signal
    [C, L] = wavedec(y, wlevel, wtype);

    % Indexes that make the following code easier
    d_L = flip(L(2:end-1));         % length of each dec
    d_i_start = d_L+1;              % start i of dec
    d_i_end = d_i_start + d_L - 1;   % end i of dec

    % Initialize arrays to store detail coefficients and thresholds for
    % plotting
    det_coeffs_array = cell(wlevel, 1);
    threshold_array = cell(wlevel, 1);

    % Apply thresholding to the detail coefficients
    for i = 1:wlevel
        % Extract detail coefficients for level i
        det_coeffs = detcoef(C, L, i);
        
        % Calculate the RMS value of the detail coefficients
        det_movrms = sqrt(movmean(det_coeffs.^ 2, movw, 'Endpoints','fill'));
        
        for j = 1:numel(det_movrms)
            if isnan(det_movrms(j))
                det_movrms(j) = rms(det_coeffs);
            end
        end
                
        high_threshold = det_movrms.*3;
        
        for j = 1:numel(C(d_i_start(i):d_i_end(i)))
            if C(d_i_start(i)+j-1) > high_threshold(j)
                C(d_i_start(i)+j-1) = high_threshold(j);
            end
        end
        
        % Set a threshold based on the RMS value
        low_threshold = det_movrms.*thresh;
        % Apply hard thresholding to the detail coefficients
        for j = 1:numel(C(d_i_start(i):d_i_end(i)))
            C(d_i_start(i)+j-1) = wthresh(C(d_i_start(i)+j-1), 'h', low_threshold(j));
        end

    % Store the detail coefficients and threshold for level i
    det_coeffs_array{i} = det_coeffs;
    threshold_array{i} = high_threshold;
    end
    
    for i = 1:wlevel
        if i ~= wlevel
            C(d_i_start:d_i_end) = 0;
        end
    end
    
    a_L = L(1);
    a_i_start = 1;
    a_i_end = a_i_start + a_L - 1;
    
    C(a_i_start:a_i_end) = 0;
    
    % Reconstruct the signal with the modified coefficients, this will be
    % the spikes we want to remoove from original signal 
    spikes = waverec(C, L, wtype);
    
    % Remove spikes from original signal
    despiked_signal = y - spikes;
    
    % Plot
    if plot_flag
        % Plot the original and despiked signals
        figure
        subplot(2, 1, 1)
        hold on
        plot(spikes)
        hold off
        legend('Spikes')
        
        subplot(2, 1, 2)
        plot(despiked_signal)
        legend('Despiked signal')

        % Plot the detail coefficients and thresholds for each level
        figure
        for i = 1:wlevel
            det_coeffs = det_coeffs_array{i};
            high_threshold = threshold_array{i};
            low_threshold = high_threshold./3.*thresh;
    
            subplot(wlevel, 1, i)
            plot(det_coeffs)
            hold on
            plot(low_threshold,'r')
            plot(-low_threshold,'r')
            plot(high_threshold,'k')
            plot(-high_threshold,'k')
            hold off
            title(['Level ' num2str(i) ' Detail Coefficients and Threshold'])
            legend('Detail Coefficients', 'pos threshold', 'neg threshold')
        end
    end
end