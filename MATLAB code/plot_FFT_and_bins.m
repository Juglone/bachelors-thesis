% == Plot range FFTs and detected range bins of each chirp ==
% This function will plot the range FFT of each chirp and circle the peak
% that has been detected as the range bin. There is a slider which makes it
% possible to move through all the chirps.

function plot_FFT_and_bins(num_frames, T_f, num_samples, fft_magnitude, bin_indices)
    % Define global variables so they can be used within callback function
    
    F_c = 1/T_f; % Sampling frequency of chirps   
    f = F_c*(0:num_samples-1)/num_samples; % X-axis of frequency values
    
    % Define the number of graphs to be plotted
    n = num_frames;

    % Create the figure window
    fig = figure();

    % Create the slider
    slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', n, 'Value', 1, ...
        'SliderStep', [1/(n-1) , 1/(n-1)], 'Position', [100 10 300 20]);

    % Create the axes
    axes_handle = axes('Units', 'normalized', 'Position', [0.1 0.3 0.8 0.6]);

    % Call the callback function to update plot based on slider value
    slider.Callback = @(~, ~) update_plot(f, fft_magnitude, bin_indices, slider);

    % Initialize the plot with the first graph
    index = 1;
    hold on
    plot(f,fft_magnitude(index,:))
    plot(f(bin_indices(index)),fft_magnitude(index,bin_indices(index)),'o')
    hold off
    title('Frequency spectrum of a chirp')
    legend('FFT', 'detected range bin')
    xlabel('Frequency [Hz]');
    ylabel('Amplitude [arb. unit]');
end