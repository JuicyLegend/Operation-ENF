% Main function to run ENF analysis with hardcoded parameters
clc; clear; close all;

% Hardcoded filenames (just the names, without the path)
filename1 = 'BASF 4-1 TMS.wav';
filename2 = 'TMS-new 32-bit PCM.wav';

% Define the common file path
filepath = 'G:\Documents\The Mysterious Song\Spectrogram Analyzer\ENF Analysis\TapeRecordings-WAV/';

% Concatenate the filepath with the filenames
fullpath1 = fullfile(filepath, filename1);
fullpath2 = fullfile(filepath, filename2);

% List of full file paths to process
filenames = {fullpath1, fullpath2}; % Add or remove file paths as needed

% Hardcoded resampling frequency
FS_resample = 1000;

nfft = 512;
%%

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
%%
    % Ensure the lengths are equal for correlation calculation
    minLength = min(length(mean_frequency1), length(mean_frequency2));
    correlation = corrcoef(mean_frequency1(1:minLength), mean_frequency2(1:minLength))
    [R,P,RL,RU] = corrcoef(mean_frequency1(1:minLength), mean_frequency2(1:minLength),'Alpha',0.05)

    % Display the correlation coefficient
    disp(['Correlation between the two signals: ', num2str(correlation(1,2))]);
    
    Fs = FS_resample;         % Sample Rate
    sig1 = mean_frequency1;
    sig2 = mean_frequency2;

    % sig1 = alignsignals(sig1,sig2);
    %[sig1,sig2,D] = alignsignals(sig1,sig2, Method="xcorr");

    figure
    ax(1) = subplot(2,1,1);
    plot(sig1)
    grid on 
    title("BASF 4-1 TMS")
    axis tight
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    ax(2) = subplot(2,1,2);
    plot(sig2)
    grid on 
    title("TMS-new 32-bit PCM")
    axis tight
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    linkaxes(ax,"xy")

     %% Auto-Covariance Sig1 & Sig2
    maxlags1 = round(numel(sig1)*0.5);
    [xc1,lag1] = xcov(sig1,maxlags1);  

    maxlags2 = round(numel(sig2)*0.5);
    [xc2,lag2] = xcov(sig2,maxlags2);  
    
    [~,df1] = findpeaks(xc1,MinPeakDistance=5*2);
    [~,mf1] = findpeaks(xc1);

    [~,df2] = findpeaks(xc2,MinPeakDistance=5*2);
    [~,mf2] = findpeaks(xc2);
    
    figure
    ax(1) = subplot(2,1,1);
    plot(lag1/(2*24),xc1,"k",...
         lag1(df1)/(2*24),xc1(df1),"kv",MarkerFaceColor="r")
    grid on
    xlim([-3 3])
    xlabel("Time")
    title("Auto-Covariance sig1")
    
    ax(2) = subplot(2,1,2);
    plot(lag2/(2*24),xc2,"k",...
         lag2(df2)/(2*24),xc2(df2),"kv",MarkerFaceColor="r")
    grid on
    xlim([-3 3])
    xlabel("Time")
    title("Auto-Covariance sig2")
    linkaxes(ax,"xy")

    %% Cross-covariance
    
    %[c,lags] = xcov(sig1,sig2, 'normalized');
    [c,lags] = xcov(sig1,sig2);
    stem(lags,c)
    
    %% Periodograms
    
    %sig2 = alignsignals(sig2,sig1);
    %length(sig2)

    sig1 = mean_frequency1;
    sig2 = mean_frequency2;

    [P1,f1] = periodogram(sig1,[],[],Fs,"power");
    [P2,f2] = periodogram(sig2,[],[],Fs,"power");
    
    figure
    t = (0:numel(sig1)-1)/Fs;
    subplot(2,2,1)
    plot(sig1,"k")
    xlabel("Time (s)")
    ylabel('Frequency (Hz)');
    grid on
    title("Time Series")
    subplot(2,2,3)
    plot(sig2)
    ylabel('Frequency (Hz)');
    grid on
    xlabel("Time (s)")
    subplot(2,2,2)
    plot(f1,log(P1),"k")
    xlabel("Frequency (Hz)")
    ylabel("Decibels (dB)")
    grid on
    axis tight
    title("Power Spectrum")
    subplot(2,2,4)
    plot(f2,log(P2))
    ylabel("Decibels (dB)")
    grid on
    axis tight
    xlabel("Frequency (Hz)")
    
    %% Coherence
    sig1 = mean_frequency1;
    sig2 = mean_frequency2;

    [Cxy,f] = mscohere(sig1,sig2,[],[],[],Fs);
    %[Cxy,f] = mscohere(sig1,sig2,[],[],2048, Fs);
    Pxy = cpsd(sig1,sig2,[],[],[],Fs);
    phase = -angle(Pxy)/pi*180;
    [pks,locs] = findpeaks(Cxy,MinPeakHeight=0.75);
    
    figure
    subplot(2,1,1)
    plot(f,Cxy)
    title("Coherence Estimate")
    grid on
    hgca = gca;
    hgca.XTick = f(locs);
    hgca.YTick = 0.75;
    axis([0 200 0 1])
    subplot(2,1,2)
    plot(f,phase)
    title("Cross-Spectrum Phase (deg)")
    grid on
    hgca = gca;
    hgca.XTick = f(locs); 
    hgca.YTick = round(phase(locs));
    xlabel("Frequency (Hz)")
    axis([0 200 -180 180])

    %% Magnitude and phase
    sig1 = mean_frequency1;
    sig2 = mean_frequency2;

    NFFT1 = length(sig1);
    Y1 = fft(sig2,NFFT1);
    F1 = ((0:1/NFFT1:1-1/NFFT1)*Fs).';

    magnitudeY1 = abs(Y1);        % Magnitude of the FFT
    phaseY1 = unwrap(angle(Y1));  % Phase of the FFT

    NFFT2 = length(sig2);
    Y2 = fft(sig2,NFFT2);
    F2 = ((0:1/NFFT2:1-1/NFFT2)*Fs).';

    magnitudeY2 = abs(Y2);        % Magnitude of the FFT
    phaseY2 = unwrap(angle(Y2));  % Phase of the FFT
    
    figure
    dB_mag1=mag2db(magnitudeY1);
    subplot(2,2,1);
    plot(F1(1:NFFT1/2),dB_mag1(1:NFFT1/2),"k");
    title('Magnitude response of signal 1');
    ylabel('Magnitude(dB)');
    subplot(2,2,2);
    plot(F1(1:NFFT1/2),phaseY1(1:NFFT1/2),"k");
    title('Phase response of signal 1');
    xlabel('Frequency in kHz');
    ylabel('radians');

    dB_mag2=mag2db(magnitudeY2);
    subplot(2,2,3);
    plot(F2(1:NFFT2/2),dB_mag2(1:NFFT2/2));
    title('Magnitude response of signal 2');
    ylabel('Magnitude(dB)');
    subplot(2,2,4);
    plot(F2(1:NFFT2/2),phaseY2(1:NFFT2/2));
    title('Phase response of signal 2');
    xlabel('Frequency in kHz');
    ylabel('radians');

    %% Measuring power
    y = sig1;

    Fs = 44100;
    NFFT = length(y);
    Fs2 = Fs/2;
    
    % Power spectrum is computed when you pass a 'power' flag input
    [P,F] = periodogram(y,[],NFFT,Fs,'power');
    
    helperFrequencyAnalysisPlot2(F,10*log10(P),'Frequency in Hz',...
      'Power spectrum (dBW)',[],[],[-0.5 Fs2])

    PdBW = 10*log10(P);
    power_at_DC_dBW = PdBW(F==0)   % dBW
    
    [peakPowers_dBW, peakFreqIdx] = findpeaks(PdBW);
    peakFreqs_Hz = F(peakFreqIdx)
    peakPowers_dBW

    SegmentLength = NFFT;

    % Power spectrum is computed when you pass a 'power' flag input
    [P,F] = pwelch(y,ones(SegmentLength,1),0,NFFT,Fs,'power');
    
    helperFrequencyAnalysisPlot2(F,10*log10(P),'Frequency in Hz',...
      'Power spectrum (dBW)',[],[],[-0.5 Fs2])

    pwr = sum(y.^2)/length(y) % in watts

    pwr1 = sum(P) % in watts

    pwr_band = bandpower(y,Fs,[50 70]);
    pwr_band_dBW = 10*log10(pwr_band) % dBW

    % Power spectral density is computed when you specify the 'psd' option
    [PSD,F]  = pwelch(y,ones(SegmentLength,1),0,NFFT,Fs,'psd');
    pwr_band1 = bandpower(PSD,F,[50 70],'psd');
    pwr_band_dBW1 = 10*log10(pwr_band1) % dBW

    y = sig2;

    Fs = 44100;
    NFFT = length(y);
    
    % Power spectrum is computed when you pass a 'power' flag input
    [P,F] = periodogram(y,[],NFFT,Fs,'power');
    
    helperFrequencyAnalysisPlot2(F,10*log10(P),'Frequency in Hz',...
      'Power spectrum (dBW)',[],[],[-0.5 Fs2])

    PdBW = 10*log10(P);
    power_at_DC_dBW = PdBW(F==0)   % dBW
    
    [peakPowers_dBW, peakFreqIdx] = findpeaks(PdBW);
    peakFreqs_Hz = F(peakFreqIdx)
    peakPowers_dBW

    SegmentLength = NFFT;

    % Power spectrum is computed when you pass a 'power' flag input
    [P,F] = pwelch(y,ones(SegmentLength,1),0,NFFT,Fs,'power');
    
    helperFrequencyAnalysisPlot2(F,10*log10(P),'Frequency in Hz',...
      'Power spectrum (dBW)',[],[],[-0.5 Fs2])

    pwr = sum(y.^2)/length(y) % in watts

    pwr1 = sum(P) % in watts

    pwr_band = bandpower(y,Fs,[50 70]);
    pwr_band_dBW = 10*log10(pwr_band) % dBW

    % Power spectral density is computed when you specify the 'psd' option
    [PSD,F]  = pwelch(y,ones(SegmentLength,1),0,NFFT,Fs,'psd');
    pwr_band1 = bandpower(PSD,F,[50 70],'psd');
    pwr_band_dBW1 = 10*log10(pwr_band1) % dBW
end

%% Auto Arrange Figures

autoArrangeFigures(2,3,1,1)
