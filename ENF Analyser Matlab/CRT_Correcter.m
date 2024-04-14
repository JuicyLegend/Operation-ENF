clc; clear; close all

% Load your audio file (ensure the path is correct and accessible)
audiofile = 'C:\Users\rutge\Downloads/TMS from Compilation A.wav';
[y, FS] = audioread(audiofile);

% Set the target frequency to align
target_frequency = 15625; % The desired frequency to align

% Shift the frequency of the audio signal
y_corrected = dynamicallyShiftSignal(y, FS, target_frequency);

% Visualization (if needed)
figure;
subplot(3,1,1);
plot(y(:,1)); % Plot first channel of original if stereo
title('Original Signal - Channel 1');
subplot(3,1,2);
plot(y_corrected(:,1)); % Plot first channel of corrected if stereo
title('Corrected Signal - Channel 1');
subplot(3,1,3);
plot(y(:,1) - y_corrected(:,1)); % Difference (error) plot
title('Difference Between Original and Corrected - Channel 1');

function y_corrected = dynamicallyShiftSignal(y, FS, target_frequency)
    % Check if the input is stereo and get the number of channels
    [numSamples, numChannels] = size(y);

    % Initialize the corrected signal
    y_corrected = zeros(numSamples, numChannels);

    % Process each channel individually
    for channel = 1:numChannels
        % Perform FFT on the current channel
        Y = fft(y(:, channel));
        n = length(Y);
        f = (0:n-1) * (FS/n); % Frequency vector

        % Initialize the corrected FFT for this channel
        Y_corrected = zeros(size(Y));

        % Determine the deviation and apply correction
        for i = 1:n
            current_frequency = f(i);
            % Calculate frequency deviation
            deviation = target_frequency - current_frequency;
            % Calculate the index shift
            shift = round(deviation / (FS/n));
            % Calculate new index
            new_index = i + shift;
            % Check for index bounds and apply shift
            if new_index > 0 && new_index <= n
                Y_corrected(new_index) = Y(i);
            end
        end

        % Convert back to time domain
        y_corrected(:, channel) = real(ifft(Y_corrected));
    end
end

% % Example usage:
% audiofile = 'C:\Users\rutge\Downloads/TMS from Compilation A.wav';
% [y, FS] = audioread(audiofile);
% target_frequency = 15625; % The desired frequency alignment
% 
% y_corrected = dynamicallyShiftSignal(y, FS, target_frequency);
% 
% figure;
% plot(target_frequency);
% % % % hold on
% plot(y);
% hold on
% plot(y_corrected);
% 
% function y_corrected = dynamicallyShiftSignal(y, FS, target_frequency)
%     % Perform FFT on the original signal
%     Y = fft(y);
%     n = length(Y);
%     f = (0:n-1) * (FS/n); % Frequency vector
% 
%     % Initialize the corrected FFT
%     Y_corrected = zeros(size(Y));
% 
%     % Determine the deviation and apply correction
%     for i = 1:n
%         current_frequency = f(i);
%         % Calculate frequency deviation
%         deviation = target_frequency - current_frequency;
%         % Calculate the index shift
%         shift = round(deviation / (FS/n));
%         % Calculate new index
%         new_index = i + shift;
%         % Check for index bounds and apply shift
%         if new_index > 0 && new_index <= n
%             Y_corrected(new_index) = Y(i);
%         end
%     end
% 
%     % Convert back to time domain
%     y_corrected = real(ifft(Y_corrected));
% end
