% clc; clear; close all;
% 
% Hardcoded filenames (just the names, without the path)
% filename1 = '001.wav';
% filename2 = '001 ref.wav';
% 
% % Define the common file path
% filepath1 = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\ENF-WHU-Dataset-master\ENF-WHU-Dataset\H1/';
% filepath2 = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\ENF-WHU-Dataset-master\ENF-WHU-Dataset\H1_ref/';
% 
% Define the common file path
% filepath = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\ENF Analyser Matlab/';
% 
% Concatenate the filepath with the filenames
% fullpath1 = fullfile(filepath, filename1);
% fullpath2 = fullfile(filepath, filename2);
% 
% List of full file paths to process
% filenames = {fullpath1, fullpath2}; % Add or remove file paths as needed
% 
% Hardcoded resampling frequency
% FS_resample = 800;
% 
% Initialize an empty cell array for the ENF data
% mean_frequencies = cell(1, length(filenames));
% 
% Process each audio file and store ENF data
% for i = 1:length(filenames)
%     Check if file exists
%     if exist(filenames{i}, 'file')
%         Call the ENF function and store the result
%         mean_frequencies{i} = ENF(filenames{i}, FS_resample);
%     else
%         disp(['File does not exist: ' filenames{i}]);
%     end
% end
% 
% Overlay the ENF plots if both files are processed
% if all(cellfun(@(x) ~isempty(x), mean_frequencies))
%     figure;
%     plot(mean_frequencies{1}, 'DisplayName', ['Signal 1: ' filename1]);
%     hold on;
%     plot(mean_frequencies{2}, 'DisplayName', ['Signal 2: ' filename2]);
%     xlabel('Time (s)');
%     ylabel('Frequency (Hz)');
%     title('Overlay of ENF - Zero Crossings');
%     legend show;
%     grid on;
% end
% 
% Calculate and display the correlation if both ENF data are available
% if length(filenames) == 2 && all(cellfun(@(f) ~isempty(f), mean_frequencies))
%     Ensure the lengths are equal for correlation calculation
%     minLength = min(length(mean_frequencies{1}), length(mean_frequencies{2}));
%     correlation = corrcoef(mean_frequencies{1}(1:minLength), mean_frequencies{2}(1:minLength));
% 
%     Display the correlation coefficient
%     disp(['Correlation between the two signals: ', num2str(correlation(1,2))]);
% end
% 
% Auto Arrange Figures
% 
% autoArrangeFigures(2,3,1,1)

% % Main function to run ENF analysis with hardcoded parameters
clc; clear; close all;

% Hardcoded filenames (just the names, without the path)
filename1 = 'BASF 4-1 TMS.wav';
filename2 = 'TMS-new 32-bit PCM.wav ';

% Define the common file path
filepath = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\TapeRecordings-WAV/';

% Concatenate the filepath with the filenames
fullpath1 = fullfile(filepath, filename1);
fullpath2 = fullfile(filepath, filename2);

% List of full file paths to process
filenames = {fullpath1, fullpath2}; % Add or remove file paths as needed

% Hardcoded resampling frequency
FS_resample = 800;

nfft = 2048;

% Process each audio file
% for i = 1:length(filenames)
%     % Check if file exists
%     if exist(filenames{i}, 'file')
%         % Call the ENF function with hardcoded inputs
%         mean_frequency = ENF(filenames{i}, FS_resample);
% 
%         % Plot the ENF signals
%         figure();
%         plot(mean_frequency, 'DisplayName', ['Signal ' num2str(i)]);
%         xlabel('Time (s)');
%         ylabel('Frequency (Hz)');
%         title(['ENF - Zero Crossings for file ' num2str(i)]);
%         legend show;
%         grid on;
% 
%         % Plot the spectrogram for each file
%         figure;
%         [y, FS] = audioread(filenames{i});
%         y = y(:)';
%         y = resample(y, FS_resample, FS);
%         spectrogram(y, nfft, nfft / 4, 2 * nfft, FS_resample, 'yaxis');
%         title(['Spectrogram for file ' num2str(i)]);
%     else
%         disp(['File does not exist: ' filenames{i}]);
%     end
% end

if exist(filenames{1,2}, 'file')
        % Call the ENF function with hardcoded inputs
        mean_frequency1 = ENF(filenames{1}, FS_resample);
        mean_frequency2 = ENF(filenames{2}, FS_resample);

        % Plot the ENF signals
        figure();
        plot(mean_frequency1, 'DisplayName', ['Signal ' num2str(1)]);
        hold on
        plot(mean_frequency2, 'DisplayName', ['Signal ' num2str(2)]);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title(['ENF - Zero Crossings for file ' num2str(1),num2str(2)]);
        legend show;
        grid on;

        % Plot the spectrogram for each file
        figure;
        [y, FS] = audioread(filenames{1});
        y = y(:)';
        y = resample(y, FS_resample, FS);
        spectrogram(y, nfft, nfft / 4, 2 * nfft, FS_resample, 'yaxis');
        title(['Spectrogram for file ' num2str(1)]);

        figure;
        [y, FS] = audioread(filenames{2});
        y = y(:)';
        y = resample(y, FS_resample, FS);
        spectrogram(y, nfft, nfft / 4, 2 * nfft, FS_resample, 'yaxis');
        title(['Spectrogram for file ' num2str(2)]);

elseif ~exist(filenames{2}, 'file')

        % Call the ENF function with hardcoded inputs
        mean_frequency = ENF(filenames{1}, FS_resample);

        % Plot the ENF signals
        figure();
        plot(mean_frequency, 'DisplayName', ['Signal ' num2str(1)]);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title(['ENF - Zero Crossings for file ' num2str(1)]);
        legend show;
        grid on;

        % Plot the spectrogram for each file
        figure;
        [y, FS] = audioread(filenames{1});
        y = y(:)';
        y = resample(y, FS_resample, FS);
        spectrogram(y, nfft, nfft / 4, 2 * nfft, FS_resample, 'yaxis');
        title(['Spectrogram for file ' num2str(1)]);
else
        disp(['File does not exist: ' filenames{1,2}]);
end

% If two files are processed, calculate and display the correlation
if length(filenames) == 2 && all(cellfun(@(f) exist(f, 'file'), filenames))
    mean_frequency1 = ENF(filenames{1}, FS_resample);
    mean_frequency2 = ENF(filenames{2}, FS_resample);

    % Ensure the lengths are equal for correlation calculation
    minLength = min(length(mean_frequency1), length(mean_frequency2));
    correlation = corrcoef(mean_frequency1(1:minLength), mean_frequency2(1:minLength));

    % Display the correlation coefficient
    disp(['Correlation between the two signals: ', num2str(correlation(1,2))]);
end


% clc; clear; close all;
% 
% filename1 = 'BASF 4-1.wav';
% % filename2 = 'TMS-new 32-bit PCM.wav';
% 
% filepath = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\TapeRecordings-WAV/';
% 
% 
% 
% % Hardcode the path to the input audio file
% Sx_in1 = fullfile(filepath, filename1);  % Replace with your audio file path
% Sx_in2 = fullfile(filepath, filename2);  % Replace with your audio file path
% 
% % Hardcode the desired resampling frequency
% FS_resample = 800;  % Replace with your desired resampling frequency
% 
% % Call the ENF function with hardcoded inputs
% mean_frequency1 = ENF(Sx_in1, FS_resample);
% mean_frequency2 = ENF(Sx_in2, FS_resample);
% 
% % Plot the ENF signals overlaying each other
% figure;
% plot(mean_frequency1, 'DisplayName', 'Signal 1');
% hold on;
% plot(mean_frequency2, 'DisplayName', 'Signal 2');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('ENF Comparison');
% legend show;
% grid on;
% 
% % Calculate the correlation between the two ENF signals
% % Ensure both signals are the same length for correlation calculation
% len = min(length(mean_frequency1), length(mean_frequency2));
% correlation = corrcoef(mean_frequency1(1:len), mean_frequency2(1:len));
% 
% % Display the correlation coefficient
% disp(['Correlation between the two signals: ', num2str(correlation(1,2))]);


% % Main function to run ENF analysis with hardcoded parameters
% clc;clear;close all;
% 
% filename1 = 'BASF 4-1 TMS.wav';
% filename2 = 'TMS-new 32-bit PCM.wav';
% 
% filepath = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\TapeRecordings-WAV/';
% 
% % Hardcode the path to the input audio file
% Sx_in1 = join(filename1, filepath);  % Replace with your audio file path
% Sx_in2 = join(filename2, filepath);  % Replace with your audio file path
% 
% % Hardcode the desired resampling frequency
% FS_resample = 1000;  % Replace with your desired resampling frequency
% 
% % Call the ENF function with hardcoded inputs
% ENF(Sx_in1, FS_resample);
% hold on
% ENF(Sx_in2, FS_resample);
% 
% % Auto Arrange Figures
% 
% autoArrangeFigures(2,3,1,1)