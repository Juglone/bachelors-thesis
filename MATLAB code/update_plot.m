% Callback function for slider
function update_plot(f, fft_magnitude, bin_indices, slider)
    % Get the value of the slider
    slider_value = round(slider.Value);
    plot(f,fft_magnitude(slider_value,:))
    hold on
    plot(f(bin_indices(slider_value)),fft_magnitude(slider_value,bin_indices(slider_value)),'o')
    hold off
    title('Frekvensspektrum för en svepsignal')
    xlabel('frekvens [Hz]');
    ylabel('Amplitud');
end