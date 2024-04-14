% Usage example
audiofile = 'C:\Users\rutge\Downloads/TMS from Compilation A.wav';
[y, FS] = audioread(audiofile);
expected_freq = 15625; % Hz
observed_min = 15890; % Hz
observed_max = 15910.5; % Hz

% Correct the drift
y_corrected = correctFrequencyDrift(y, FS, expected_freq, observed_min, observed_max);

figure;
plot(expected_freq);
hold on
plot(y);
hold on
plot(y_corrected);

function y_shifted = correctFrequencyDrift(y, FS, expected_freq, observed_min, observed_max)
    % Calculate the average observed frequency
    observed_freq = (observed_min + observed_max) / 2;
    
    % Calculate the frequency shift needed
    freq_shift = expected_freq - observed_freq;
    
    % Perform FFT on the original signal
    Y = fft(y);
    
    % Generate frequency vector for the FFT
    n = length(y); % Length of the signal
    f = (0:n-1)*(FS/n); % Frequency range - 0 to Nyquist
    
    % Shift the FFT
    df = freq_shift * (n/FS); % Convert freq_shift in Hz to bins
    Y_shifted = Y .* exp(-1j*2*pi*df*(0:n-1)/n); % Apply shift
    
    % Convert back to time domain
    y_shifted = real(ifft(Y_shifted));
end
