function mainENFAnalysis()
    % Define the audio file path
    audioFilePath = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\TapeRecordings-WAV/BASF 4-1 TMS.wav'; % Update this path to your audio file's location

    % Define parameters
    FS_resample = 1000; % Resampling frequency
    targetFreq = 50; % Target ENF frequency, e.g., 50 Hz or 60 Hz
    frameSize = 1024; % Frame size for STFT
    overlapSize = 512; % Overlap size for STFT
    freqRange = 1; % Frequency range to search around the target frequency

    % Load and resample the audio signal
    [audio, fs] = audioread(audioFilePath);
    audio = resample(audio(:,1), FS_resample, fs); % Resampling and taking the first channel if stereo

    % Extract the ENF signal using the STFT-based method
    enf_signal = refineENF_STFT(audio, FS_resample, targetFreq, frameSize, overlapSize, freqRange);

    % Plot the extracted ENF signal
    figure;
    plot(enf_signal);
    xlabel('Frame Index');
    ylabel('Frequency (Hz)');
    title('Extracted ENF Signal');
    grid on;
end

function enf = refineENF_STFT(signal, fs, targetFreq, frameSize, overlapSize, freqRange)
    % signal: Input signal
    % fs: Sampling frequency
    % targetFreq: Nominal frequency of the ENF signal (e.g., 50 Hz or 60 Hz)
    % frameSize: The number of samples in each frame
    % overlapSize: The number of samples to overlap between frames
    % freqRange: Range to search for peak frequency around the target frequency

    % Initialize the output ENF signal
    enf = [];

    % Apply STFT
    [S, F, T] = spectrogram(signal, frameSize, overlapSize, [], fs);

    % Convert freqRange to index
    freqIndices = find(F >= targetFreq - freqRange & F <= targetFreq + freqRange);

    % Iterate over all time bins
    for t = 1:length(T)
        % Extract the segment of interest
        segment = abs(S(freqIndices, t));

        % Find the peak within the specified frequency range
        [~, peakIndex] = max(segment);
        peakFreq = F(freqIndices(peakIndex));

        % Append the peak frequency to the ENF signal
        enf = [enf, peakFreq];
    end
end
