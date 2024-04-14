clc; clear; close all

audiofile1 = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\TapeRecordings-WAV/BASF 4-1 TMS.wav'; % Define the path to the audio file
audiofile2 = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\TapeRecordings-WAV/TMS-new 32-bit PCM.wav';
audiofile3 = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\TapeRecordings-WAV/TMS from Compilation A.wav';


Fs = 44100;



% [total_seconds1, mean_frequency1, y1, nfft1] = ENF(audiofile1,Fs);
% [total_seconds2, mean_frequency2, y2, nfft2] = ENF(audiofile2,Fs);
[total_seconds3, mean_frequency3, y3, nfft3] = ENF(audiofile3,Fs);
%ENF(audiofile1,Fs);

figure
subplot(2,1,1)
% %subplot(4,1,1)
% % Plot_ENF(total_seconds1, mean_frequency1)
% % hold on
% % Plot_ENF(total_seconds2, mean_frequency2)
% % hold on
Plot_ENF(total_seconds3, mean_frequency3)
% %subplot(4,1,2)
% Spectrogram_ENF(y1,nfft1,Fs)
% %subplot(4,1,3)
% Spectrogram_ENF(y2,nfft2,Fs)
%subplot(4,1,4)
subplot(2,1,2)
Spectrogram_ENF(y3,nfft3,Fs)

autoArrangeFigures(2,2,1,1)

function [total_seconds, mean_frequency, y, nfft] = ENF(Sx_in, FS_resample)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_order = 4;

% lower_bound = 49; % Hz
% upper_bound = 51; % Hz

lower_bound = 15500; %Hz
upper_bound = 16500; %Hz

% 1. Extract the samples and the sampling frequency of the audio signal
[y, FS] = audioread(Sx_in);
nfft = 4096;

% 2. Ensure that it's a row vector
y = y(:)';

% 3. Calculate the total duration of the recording in seconds
total_seconds = length(y) / FS;

% 4. Resample the signal to the resampling frequency
y = resample(y, FS_resample, FS);

% 5. Create a butterworth bandpass filter
[b, a] = butter(b_order, [(lower_bound / (FS_resample / 2)), (upper_bound / (FS_resample / 2))]); % 4*8 dB = 24 dB/octave

% 6. Filter the signal
y = filter(b, a, y);

% 7. Calculate the number of blocks
number_of_blocks = length(y) / FS_resample;

% 8. Create an auxiliary variable to traverse the "aux" vector
k = 1; mean_frequency = zeros(1, ceil(number_of_blocks));

%%%%%%%%%%% ENF by the method of - Linear Interpolation - %%%%%%%%%%%%%%
for e = 1: ceil(number_of_blocks)
    
    % 9. Define the start of each block
    start_block = (e - 1) * FS_resample + 1;
    
    % 10. Define the end of each block
    end_block = (start_block + FS_resample) - 1;
    
    % 11. For the case of the last block, take up to where the samples reach
    if e > number_of_blocks
        end_block = length(y); % Define the end of this block
    end
    
    % 12. Extract from vector "y" the portion of each block in question
    block = y(1, start_block:end_block);
    
    % 13. Traverse the block looking for zero crossings
    for i = 1:length(block) - 1
        if (block(i) * block(i + 1) < 0)
            
            % Formula for linear interpolation
            aux(k) = (((-block(i)) * (i + 1 - i)) / (block(i + 1) - block(i))) + i;
            % Advance k
            k = k + 1;
        end
    end

    % 14. Find the distance between zeros
    for i = 2:length(aux)
        samples_between_zeros(i - 1) = aux(i) - aux(i - 1);
    end
    
    % 15. Convert the interval from samples to seconds
    t = samples_between_zeros / FS_resample;
    
    % 16. Calculate the mean frequency of the block
    mean_frequency(e) = 1 / mean(t * 2);
    
    % Reset k to its initial value
    k = 1;
end

return

% % 17. Display ENF according to zero crossings
% axis = linspace(0, total_seconds, length(mean_frequency));
% figure
% plot(axis, mean_frequency)
% axis([0 total_seconds 48 52])
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Zero Crossings');
% grid
% 
% % 18. Display ENF with the spectrogram
% figure
% spectrogram(y, 512, 512 / 4, 2 * 512, FS_resample, 'yaxis');

% % 17. Display ENF according to zero crossings
% time_axis = linspace(0, total_seconds, length(mean_frequency));
% figure;
% plot(time_axis, mean_frequency);
% if total_seconds > 0  % Ensure total_seconds is positive
%     axis([0 total_seconds 48 52]);
%     %axis([0 total_seconds 45 55]);
% else
%     warning('Total seconds is non-positive, cannot set axis limits.');
% end
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Zero Crossings');
% grid on;

% % 18. Display ENF with the spectrogram
% figure;
% spectrogram(y, nfft, nfft / 4, 2 * nfft, FS_resample, 'yaxis');

end

function Plot_ENF(total_seconds, mean_frequency)
zoom = 100; % Hz
% lower_bound = 49; % Hz
% upper_bound = 51; % Hz

lower_bound = 15500; %Hz
upper_bound = 16500; %Hz

% 17. Display ENF according to zero crossings
time_axis = linspace(0, total_seconds, length(mean_frequency));
%figure;
plot(time_axis, mean_frequency);
if total_seconds > 0  % Ensure total_seconds is positive
    axis([0 total_seconds (lower_bound-zoom) (upper_bound+zoom)]);
    %axis([0 total_seconds 45 55]);
else
    warning('Total seconds is non-positive, cannot set axis limits.');
end
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Zero Crossings');
grid on;

end

function Spectrogram_ENF(y,nfft,FS_resample)

% 18. Display ENF with the spectrogram
%figure;
spectrogram(y, nfft, nfft / 4, 2 * nfft, FS_resample, 'yaxis');

end