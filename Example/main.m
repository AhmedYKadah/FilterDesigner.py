clear; close all; clc;

% Create results folder if it doesn't exist
if ~exist('results', 'dir')
    mkdir('results');
    fprintf('Created results folder\n');
end

% 1. Load EEG Signal from MAT file
% Specify your MAT file path
eeg_file = 'data/chb12_29_data.mat';  % Your EEG data file

% Load MAT file
eeg_data = load(eeg_file);

% Display available variables to help identify the signal
disp('Variables in EEG data file:');
disp(fieldnames(eeg_data));

% Extract signal - common variable names in EEG .mat files
% Adjust these based on your actual variable names
if isfield(eeg_data, 'data')
    eeg_signal = eeg_data.data;
elseif isfield(eeg_data, 'signal')
    eeg_signal = eeg_data.signal;
elseif isfield(eeg_data, 'eeg')
    eeg_signal = eeg_data.eeg;
elseif isfield(eeg_data, 'val')
    eeg_signal = eeg_data.val;
else
    % Get the first field if none of the above exist
    field_names = fieldnames(eeg_data);
    eeg_signal = eeg_data.(field_names{1});
end

% Display signal dimensions
fprintf('Original signal dimensions: %d x %d\n', size(eeg_signal, 1), size(eeg_signal, 2));

% Handle multi-channel data correctly
% EEG data is typically: channels x samples or samples x channels
if size(eeg_signal, 1) > size(eeg_signal, 2)
    % More rows than columns - likely samples x channels, transpose it
    fprintf('Data appears to be samples x channels. Transposing...\n');
    eeg_signal = eeg_signal';
end

if size(eeg_signal, 1) > 1
    fprintf('Multi-channel data detected (%d channels). Using channel 1.\n', size(eeg_signal, 1));
    eeg_signal = eeg_signal(1, :);  % First channel
else
    fprintf('Single channel data detected.\n');
end

% Ensure row vector initially for easier processing
eeg_signal = eeg_signal(:)';
fprintf('Signal length before cleaning: %d samples\n', length(eeg_signal));

% Check for and handle non-finite values (NaN, Inf)
num_nans = sum(isnan(eeg_signal));
num_infs = sum(isinf(eeg_signal));
num_non_finite = num_nans + num_infs;

fprintf('\n--- Data Quality Check ---\n');
fprintf('Total samples: %d\n', length(eeg_signal));
fprintf('NaN values: %d (%.4f%%)\n', num_nans, 100 * num_nans / length(eeg_signal));
fprintf('Inf values: %d (%.4f%%)\n', num_infs, 100 * num_infs / length(eeg_signal));
fprintf('Total non-finite: %d (%.4f%%)\n', num_non_finite, ...
        100 * num_non_finite / length(eeg_signal));

if num_non_finite > 0
    fprintf('\nRemoving non-finite values from signal...\n');
    
    % Remove all non-finite values (NaN and Inf)
    valid_idx = isfinite(eeg_signal);
    eeg_signal = eeg_signal(valid_idx);
    
    fprintf('Cleaned signal length: %d samples\n', length(eeg_signal));
    fprintf('Removed: %d samples\n', num_non_finite);
else
    fprintf('Signal is clean - no non-finite values detected\n');
end

% Ensure column vector for filtering
eeg_signal = eeg_signal(:);

% Extract or set sampling frequency
if isfield(eeg_data, 'fs')
    fs = eeg_data.fs;
elseif isfield(eeg_data, 'Fs')
    fs = eeg_data.Fs;
elseif isfield(eeg_data, 'samplingRate')
    fs = eeg_data.samplingRate;
else
    % Default for CHB-MIT dataset is 256 Hz
    fs = 256;
    fprintf('Sampling frequency not found. Using default: %.2f Hz\n', fs);
end

fprintf('EEG data loaded successfully\n');
fprintf('Sampling frequency: %.2f Hz\n', fs);
fprintf('Signal length: %d samples (%.2f seconds)\n', ...
        length(eeg_signal), length(eeg_signal)/fs);

% 2. Load IIR Filter Coefficients
% Specify your MAT file path
filter_file = 'filters/Notching.mat';  % Your filter coefficients file

try
    filter_data = load(filter_file);
    
    % Extract filter coefficients
    % Assuming the filter coefficients are stored as 'b' and 'a'
    % Modify these variable names based on your .mat file structure
    if isfield(filter_data, 'b') && isfield(filter_data, 'a')
        b = filter_data.b;
        a = filter_data.a;
    elseif isfield(filter_data, 'Num') && isfield(filter_data, 'Den')
        b = filter_data.Num;
        a = filter_data.Den;
    else
        % Display available variables in the .mat file
        disp('Variables in filter file:');
        disp(fieldnames(filter_data));
        error('Please modify the code to match your filter coefficient variable names');
    end
    
    fprintf('Filter coefficients loaded successfully\n');
    fprintf('Numerator order: %d\n', length(b)-1);
    fprintf('Denominator order: %d\n', length(a)-1);
catch
    error('Failed to load filter file. Check filename and path.');
end

% 3. Apply IIR Filter to EEG Signal
% Apply zero-phase filtering to avoid phase distortion
eeg_filtered = filtfilt(b, a, eeg_signal);

fprintf('Filtering completed\n');

% 4. Plot Original vs Filtered Signal
t = (0:length(eeg_signal)-1) / fs;  % Time vector

fig1 = figure('Position', [100, 100, 1200, 600]);

% Time domain comparison
subplot(2,1,1);
plot(t, eeg_signal, 'b', 'LineWidth', 1);
hold on;
plot(t, eeg_filtered, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title('EEG Signal: Original vs Filtered');
legend('Original', 'Filtered');
grid on;
xlim([0, min(10, t(end))]);  % Show first 10 seconds

% Zoomed view
subplot(2,1,2);
zoom_duration = 2;  % seconds
zoom_samples = min(zoom_duration * fs, length(eeg_signal));
plot(t(1:zoom_samples), eeg_signal(1:zoom_samples), 'b', 'LineWidth', 1);
hold on;
plot(t(1:zoom_samples), eeg_filtered(1:zoom_samples), 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title(sprintf('Zoomed View (First %d seconds)', zoom_duration));
legend('Original', 'Filtered');
grid on;

% Save figure
saveas(fig1, 'results/time_domain_comparison.png');

% 5. Frequency Domain Analysis - Welch's Method

% Welch's method parameters - adjust based on signal length
signal_length = length(eeg_signal);
window_length = min(2 * fs, floor(signal_length / 4));  % 2 seconds or 1/4 signal length
overlap = floor(window_length / 2);  % 50% overlap
nfft = max(256, 2^nextpow2(window_length));

fprintf('\nWelch method parameters:\n');
fprintf('Window length: %d samples (%.2f seconds)\n', window_length, window_length/fs);
fprintf('Overlap: %d samples\n', overlap);
fprintf('NFFT: %d\n', nfft);

% Compute Power Spectral Density
[pxx_orig, f_welch] = pwelch(eeg_signal, hamming(window_length), ...
                              overlap, nfft, fs);
[pxx_filt, ~] = pwelch(eeg_filtered, hamming(window_length), ...
                       overlap, nfft, fs);

% Convert to dB scale
pxx_orig_db = 10*log10(pxx_orig);
pxx_filt_db = 10*log10(pxx_filt);

fig2 = figure('Position', [100, 100, 1200, 800]);

% PSD comparison
subplot(2,1,1);
plot(f_welch, pxx_orig_db, 'b', 'LineWidth', 1.5);
hold on;
plot(f_welch, pxx_filt_db, 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density - Welch Method');
legend('Original', 'Filtered');
grid on;
xlim([0, min(130, fs/2)]);  % Show up to 130 Hz or Nyquist frequency

% Filter frequency response
subplot(2,1,2);
[h, f_response] = freqz(b, a, nfft, fs);
plot(f_response, 20*log10(abs(h)), 'k', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('IIR Filter Frequency Response');
grid on;
xlim([0, min(130, fs/2)]);

% Save figure
saveas(fig2, 'results/frequency_analysis_welch.png');

% 6. Frequency Domain Analysis - Spectrogram (STFT)

% Check if signal is long enough for spectrogram
min_samples_for_spectrogram = fs * 1;  % At least 1 second

if length(eeg_signal) < min_samples_for_spectrogram
    fprintf('\nWarning: Signal too short for spectrogram analysis (%.2f seconds)\n', ...
            length(eeg_signal)/fs);
    fprintf('Skipping spectrogram...\n');
else
    % Spectrogram parameters - adaptive based on signal length
    window_spec = hamming(min(round(1 * fs), floor(length(eeg_signal)/4)));  % 1-second window or 1/4 signal
    overlap_spec = floor(0.9 * length(window_spec));  % 90% overlap
    nfft_spec = max(512, 2^nextpow2(length(window_spec)));
    
    fprintf('\nSpectrogram parameters:\n');
    fprintf('Window length: %d samples (%.2f seconds)\n', length(window_spec), length(window_spec)/fs);
    fprintf('Overlap: %d samples\n', overlap_spec);
    fprintf('NFFT: %d\n', nfft_spec);
    
    fig3 = figure('Position', [100, 100, 1200, 800]);
    
    % Original signal spectrogram
    subplot(2,1,1);
    spectrogram(eeg_signal, window_spec, overlap_spec, nfft_spec, fs, 'yaxis');
    title('Spectrogram - Original EEG Signal');
    colorbar;
    caxis([-60, 20]);  % Adjust color scale as needed
    ylim([0, min(130, fs/2)]);  % Show up to 130 Hz
    
    % Filtered signal spectrogram
    subplot(2,1,2);
    spectrogram(eeg_filtered, window_spec, overlap_spec, nfft_spec, fs, 'yaxis');
    title('Spectrogram - Filtered EEG Signal');
    colorbar;
    caxis([-60, 20]);  % Adjust color scale as needed
    ylim([0, min(130, fs/2)]);  % Show up to 130 Hz
    
    % Save figure
    saveas(fig3, 'results/spectrogram_analysis.png');
end

% 7. Statistical Comparison
fprintf('\n--- Signal Statistics ---\n');
fprintf('Original Signal:\n');
fprintf('  Mean: %.4f µV\n', mean(eeg_signal));
fprintf('  Std Dev: %.4f µV\n', std(eeg_signal));
fprintf('  RMS: %.4f µV\n', rms(eeg_signal));

fprintf('\nFiltered Signal:\n');
fprintf('  Mean: %.4f µV\n', mean(eeg_filtered));
fprintf('  Std Dev: %.4f µV\n', std(eeg_filtered));
fprintf('  RMS: %.4f µV\n', rms(eeg_filtered));

% Compute differences
mean_diff = mean(eeg_filtered) - mean(eeg_signal);
std_diff = std(eeg_filtered) - std(eeg_signal);
rms_diff = rms(eeg_filtered) - rms(eeg_signal);

fprintf('\n--- Differences (Filtered - Original) ---\n');
fprintf('Mean difference: %.4f µV\n', mean_diff);
fprintf('Std Dev difference: %.4f µV\n', std_diff);
fprintf('RMS difference: %.4f µV\n', rms_diff);

% Compute Mean Squared Error
mse = mean((eeg_filtered - eeg_signal).^2);
fprintf('\nMean Squared Error: %.4f µV²\n', mse);
fprintf('Root Mean Squared Error: %.4f µV\n', sqrt(mse));

% Compute noise reduction (assuming filtering removes noise)
noise_power_orig = bandpower(eeg_signal, fs, [50, min(100, fs/2)]);
noise_power_filt = bandpower(eeg_filtered, fs, [50, min(100, fs/2)]);
noise_reduction_db = 10*log10(noise_power_orig / noise_power_filt);

fprintf('\nNoise Reduction (50-100 Hz): %.2f dB\n', noise_reduction_db);

fprintf('\nAnalysis completed successfully!\n');
fprintf('All figures saved to results/ folder\n');
