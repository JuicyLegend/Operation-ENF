clc; clear; close all;

% Define the file paths
filepath1 = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\ENF-WHU-Dataset-master\ENF-WHU-Dataset\H1/';
filepath2 = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\ENF-WHU-Dataset-master\ENF-WHU-Dataset\H1_ref/';
filename1 = '001.wav';  % First audio file
filename2 = '001_ref.wav';  % Second audio file (optional)

% Define the common file path and filenames
fullpath1 = fullfile(filepath1, filename1);
fullpath2 = fullfile(filepath2, filename2);

% Hardcoded resampling frequency
FS_resample = 800;

% Initialize an empty cell array for the ENF data
mean_frequencies = cell(1, 2);

% Process the first audio file and store ENF data
mean_frequencies{1} = ENF(fullpath1, FS_resample);

% Check if the second audio file exists and process it
if exist(fullpath2, 'file')
    mean_frequencies{2} = ENF(fullpath2, FS_resample);
end

% Overlay the ENF plots if both files are processed
figure;
plot(mean_frequencies{1}, 'DisplayName', ['ENF from ' filename1]);
hold on;
if ~isempty(mean_frequencies{2})
    plot(mean_frequencies{2}, 'DisplayName', ['ENF from ' filename2]);
end
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Overlay of ENF - Zero Crossings');
legend show;
grid on;

% Function to extract ENF from the audio signal
function mean_frequency = ENF(Sx_in, FS_resample)
    % Extract the samples and the sampling frequency of the audio signal
    [y, FS] = audioread(Sx_in);
    
    % Ensure that it's a row vector
    y = y(:)';
    
    % Resample the signal to the resampling frequency
    y = resample(y, FS_resample, FS);
    
    % Create a butterworth bandpass filter
    [b, a] = butter(4, [(49 / (FS_resample / 2)), (51 / (FS_resample / 2))]);
    
    % Filter the signal
    y = filter(b, a, y);
    
    % Initialize mean_frequency as a zeros row vector
    number_of_blocks = length(y) / FS_resample;
    mean_frequency = zeros(1, ceil(number_of_blocks));
    
    % ENF extraction by the method of - Linear Interpolation -
    for e = 1: ceil(number_of_blocks)
        % Define the start and end of each block
        start_block = (e - 1) * FS_resample + 1;
        end_block = min(start_block + FS_resample - 1, length(y));
        
        % Extract the block from the signal
        block = y(start_block:end_block);
        
        % Initialize aux for this block
        aux = [];
        
        % Traverse the block looking for zero crossings
        for i = 1:length(block) - 1
            if (block(i) * block(i + 1) < 0)
                % Formula for linear interpolation
                aux = [aux, (((-block(i)) * (i + 1 - i)) / (block(i + 1) - block(i))) + i];
            end
        end
        
        % Find the distance between zeros and calculate the mean frequency
        if ~isempty(aux) && length(aux) > 1
            samples_between_zeros = diff(aux);
            t = samples_between_zeros / FS_resample;
            mean_frequency(e) = 1 / mean(t * 2);
        end
    end
end
