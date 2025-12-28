
clear; close all; clc;

% Configuration
% Directory containing .mat files
filter_dir = '../Filters'; 

% Specifications
specs.stopband_atten = 10;    % dB (minimum attenuation in stopband)
specs.passband_ripple = 1.5;  % dB (maximum ripple in passband)
specs.fc = 60;                % Hz (center frequency to reject)
specs.stopband = [50, 70];    % Hz (stopband region)
specs.fs = 256;              % Hz (sampling frequency)

% Frequency ranges for analysis
freq_ranges.passband_low = [0, 45];      % Hz (lower passband)
freq_ranges.passband_high = [75, specs.fs/2];  % Hz (upper passband)
freq_ranges.stopband = specs.stopband;    % Hz

% Number of frequency points for analysis (higher = more accurate)
nfft = 8192;  % Increased for better frequency resolution

% Find all .mat files
mat_files = dir(fullfile(filter_dir, '*.mat'));
num_filters = length(mat_files);

if num_filters == 0
    error('No .mat files found in directory: %s', filter_dir);
end

fprintf('Found %d filter files to validate\n\n', num_filters);

% Initialize results storage
results = struct('filename', {}, 'passed', {}, 'stopband_met', {}, ...
                 'passband_met', {}, 'max_stopband_atten', {}, ...
                 'max_passband_ripple', {}, 'comments', {});

% Iterate through each filter
for idx = 1:num_filters
    filename = mat_files(idx).name;
    filepath = fullfile(filter_dir, filename);
    
    fprintf('===== Validating Filter %d/%d: %s =====\n', idx, num_filters, filename);
    
    try
        % Load filter coefficients
        filter_data = load(filepath);
        
        % Extract coefficients (try common variable names)
        % Initialize a and b
        b = [];
        a = [];
        
        if isfield(filter_data, 'Num') && isfield(filter_data, 'Den')
            % IIR format with Num/Den
            b = filter_data.Num;
            a = filter_data.Den;
        elseif isfield(filter_data, 'Num') && ~isfield(filter_data, 'Den')
            % FIR format with only Num (no Den field)
            b = filter_data.Num;
            a = 1;
        elseif isfield(filter_data, 'b') && isfield(filter_data, 'a')
            % IIR format with b/a
            b = filter_data.b;
            a = filter_data.a;
        elseif isfield(filter_data, 'b') && ~isfield(filter_data, 'a')
            % FIR format with only b (no a field)
            b = filter_data.b;
            a = 1;
        elseif isfield(filter_data, 'SOS')
            % Second-order sections format (IIR)
            sos = filter_data.SOS;
            if isfield(filter_data, 'G')
                g = filter_data.G;
            else
                g = 1;
            end
            [b, a] = sos2tf(sos, g);
        elseif isfield(filter_data, 'Numerator')
            % FIR format with Numerator field
            b = filter_data.Numerator;
            a = 1;
        else
            % Try to find coefficient arrays automatically
            field_names = fieldnames(filter_data);
            if length(field_names) >= 2
                b = filter_data.(field_names{1});
                a = filter_data.(field_names{2});
            elseif length(field_names) == 1
                % Single field - assume FIR
                b = filter_data.(field_names{1});
                a = 1;
            else
                error('Could not find filter coefficients in expected format');
            end
        end
        
        % Ensure b and a are row vectors
        b = b(:).';
        if length(a) == 1
            a = a(:).';
        else
            a = a(:).';
        end
        
        % Determine filter type
        if length(a) == 1 && a == 1
            filter_type = 'FIR';
        else
            filter_type = 'IIR';
        end
        
        % Calculate filter complexity metrics
        num_coeffs_b = length(b);
        num_coeffs_a = length(a);
        
        if strcmp(filter_type, 'FIR')
            % FIR: y[n] = b0*x[n] + b1*x[n-1] + ... + bN*x[n-N]
            % Multiplications: N (one per coefficient, excluding zero coefficients)
            % Additions: N-1 (to sum all products)
            non_zero_coeffs = sum(abs(b) > eps);
            num_multipliers = non_zero_coeffs;
            num_adders = max(0, non_zero_coeffs - 1);
            total_ops = num_multipliers + num_adders;
            
            complexity_str = sprintf('FIR: %d coeffs, %d mults, %d adds, %d total ops', ...
                                    num_coeffs_b, num_multipliers, num_adders, total_ops);
        else
            % IIR: y[n] = (b0*x[n] + b1*x[n-1] + ... ) - (a1*y[n-1] + a2*y[n-2] + ...)
            % Numerator (feedforward): length(b) multiplications, length(b)-1 additions
            % Denominator (feedback): length(a)-1 multiplications (a0 is normalized to 1), length(a)-2 additions
            non_zero_b = sum(abs(b) > eps);
            non_zero_a = sum(abs(a(2:end)) > eps);  % Exclude a[0] = 1
            
            num_multipliers = non_zero_b + non_zero_a;
            num_adders = max(0, non_zero_b - 1) + max(0, non_zero_a - 1) + 1; % +1 for combining feedforward and feedback
            total_ops = num_multipliers + num_adders;
            
            complexity_str = sprintf('IIR: %d num coeffs, %d den coeffs, %d mults, %d adds, %d total ops', ...
                                    num_coeffs_b, num_coeffs_a, num_multipliers, num_adders, total_ops);
        end
        
        fprintf('  Filter type: %s\n', filter_type);
        fprintf('  Complexity: %s\n', complexity_str);
        
        % Compute frequency response with high resolution
        [H, f] = freqz(b, a, nfft, specs.fs);
        H_dB = 20*log10(abs(H));
        
        % Check filter stability (only for IIR filters)
        is_stable = true;  % FIR filters are always stable
        if strcmp(filter_type, 'IIR')
            % Also compute group delay to check for stability issues
            [gd, f_gd] = grpdelay(b, a, nfft, specs.fs);
            
            % Check filter stability
            poles = roots(a);
            is_stable = all(abs(poles) < 1);
            if ~is_stable
                warning('Filter %s may be unstable (poles outside unit circle)', filename);
            end
        end
        
        % Analyze stopband performance at 60 Hz specifically
        % Find the exact response at 60 Hz using interpolation
        H_at_60Hz_dB = interp1(f, H_dB, specs.fc, 'spline');
        stopband_atten = -H_at_60Hz_dB;  % Attenuation at 60 Hz (positive value)
        
        % Also get general stopband statistics for reference
        stopband_idx = (f >= specs.stopband(1)) & (f <= specs.stopband(2));
        stopband_response = H_dB(stopband_idx);
        max_stopband_level = max(stopband_response);
        min_stopband_level = min(stopband_response);
        
        % Report worst case in stopband
        [~, worst_idx] = max(stopband_response);
        stopband_freqs = f(stopband_idx);
        worst_freq = stopband_freqs(worst_idx);
        
        fprintf('  Stopband Analysis:\n');
        fprintf('    Level at 60 Hz: %.2f dB (Attenuation: %.2f dB)\n', H_at_60Hz_dB, stopband_atten);
        fprintf('    Worst case in stopband: %.2f dB @ %.2f Hz\n', max_stopband_level, worst_freq);
        fprintf('    Best attenuation in stopband: %.2f dB\n', -min_stopband_level);
        
        % Analyze passband performance (both lower and upper passbands)
        passband_low_idx = (f >= freq_ranges.passband_low(1)) & (f <= freq_ranges.passband_low(2));
        passband_high_idx = (f >= freq_ranges.passband_high(1)) & (f <= freq_ranges.passband_high(2));
        passband_idx = passband_low_idx | passband_high_idx;
        
        passband_response = H_dB(passband_idx);
        passband_max = max(passband_response);
        passband_min = min(passband_response);
        passband_ripple_magnitude = passband_max - passband_min;
        
        % Check specifications
        % Ripple constraint: response must stay within [-1.5dB, +1.5dB] bounds
        stopband_met = (stopband_atten >= specs.stopband_atten);
        passband_met = (passband_max <= specs.passband_ripple) && (passband_min >= -specs.passband_ripple);
        all_passed = stopband_met && passband_met;
        
        % Store results
        results(idx).filename = filename;
        results(idx).passed = all_passed;
        results(idx).stopband_met = stopband_met;
        results(idx).passband_met = passband_met;
        results(idx).max_stopband_atten = stopband_atten;
        results(idx).passband_max = passband_max;
        results(idx).passband_min = passband_min;
        results(idx).passband_ripple_magnitude = passband_ripple_magnitude;
        results(idx).is_stable = is_stable;
        results(idx).H_at_60Hz_dB = H_at_60Hz_dB;
        results(idx).worst_stopband_freq = worst_freq;
        results(idx).worst_stopband_level = max_stopband_level;
        results(idx).filter_type = filter_type;
        results(idx).num_coeffs_b = num_coeffs_b;
        results(idx).num_coeffs_a = num_coeffs_a;
        results(idx).num_multipliers = num_multipliers;
        results(idx).num_adders = num_adders;
        results(idx).total_ops = total_ops;
        
        % Print results
        fprintf('  Stopband Attenuation @ 60 Hz: %.2f dB (Required: %.2f dB) - %s\n', ...
                stopband_atten, specs.stopband_atten, ...
                iif(stopband_met, 'PASS', 'FAIL'));
        fprintf('  Passband Ripple Bounds: [%.2f, %.2f] dB (Required: [%.2f, %.2f] dB) - %s\n', ...
                passband_min, passband_max, -specs.passband_ripple, specs.passband_ripple, ...
                iif(passband_met, 'PASS', 'FAIL'));
        fprintf('    (Ripple magnitude: %.2f dB)\n', passband_ripple_magnitude);
        fprintf('  Overall: %s\n', iif(all_passed, '✓ PASS', '✗ FAIL'));
        if strcmp(filter_type, 'IIR') && ~is_stable
            fprintf('  WARNING: IIR filter may be unstable!\n');
        end
        
        % Compute phase response
        H_phase = angle(H);
        H_phase_deg = unwrap(H_phase) * 180 / pi;
        
        % Generate individual plot
        figure('Position', [100, 100, 1400, 1000]);
        
        subplot(4,2,1);
        plot(f, H_dB, 'b', 'LineWidth', 1.5);
        hold on;
        
        % Mark specification regions
        yline(0, 'k--', 'LineWidth', 1);
        yline(-specs.stopband_atten, 'r--', 'LineWidth', 1.5, ...
              'Label', sprintf('Min Stopband: -%.1f dB', specs.stopband_atten));
        
        % Shade stopband region
        ylims = ylim;
        fill([specs.stopband(1) specs.stopband(2) specs.stopband(2) specs.stopband(1)], ...
             [ylims(1) ylims(1) ylims(2) ylims(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        title(sprintf('Filter: %s (%s) - %s | Ops: %d mults + %d adds = %d total', ...
                      filename, filter_type, iif(all_passed, 'PASS', 'FAIL'), ...
                      num_multipliers, num_adders, total_ops));
        legend('Frequency Response', '0 dB', 'Stopband Spec', 'Stopband Region');
        xlim([0 specs.fs/2]);
        
        subplot(4,2,2);
        plot(f, H_phase_deg, 'b', 'LineWidth', 1.5);
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Phase (degrees)');
        title('Phase Response - Full Range');
        xlim([0 specs.fs/2]);
        
        subplot(4,2,3);
        plot(f, H_dB, 'b', 'LineWidth', 1.5);
        hold on;
        
        % Mark the specific frequencies
        plot(worst_freq, max_stopband_level, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        plot(60, H_at_60Hz_dB, 'mx', 'MarkerSize', 12, 'LineWidth', 3);
        
        yline(0, 'k--', 'LineWidth', 1);
        yline(-specs.stopband_atten, 'r--', 'LineWidth', 1.5);
        
        % Shade stopband region
        ylims = ylim;
        fill([specs.stopband(1) specs.stopband(2) specs.stopband(2) specs.stopband(1)], ...
             [ylims(1) ylims(1) ylims(2) ylims(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        title(sprintf('Stopband Detail - 60Hz: %.2f dB | Worst: %.2f dB @ %.2f Hz', ...
                      H_at_60Hz_dB, max_stopband_level, worst_freq));
        legend('Response', 'Worst Case', '60 Hz (TEST POINT)', '0 dB', 'Stopband Spec');
        xlim([40 80]);
        
        subplot(4,2,4);
        plot(f, H_phase_deg, 'b', 'LineWidth', 1.5);
        hold on;
        
        % Shade stopband region
        ylims_phase = ylim;
        fill([specs.stopband(1) specs.stopband(2) specs.stopband(2) specs.stopband(1)], ...
             [ylims_phase(1) ylims_phase(1) ylims_phase(2) ylims_phase(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Phase (degrees)');
        title('Phase Response - Stopband Detail');
        xlim([40 80]);
        
        subplot(4,2,5);
        plot(f, H_dB, 'b', 'LineWidth', 1.5);
        hold on;
        yline(0, 'k--', 'LineWidth', 1);
        yline(-specs.passband_ripple, 'g--', 'LineWidth', 1.5, 'Label', sprintf('-%.1f dB', specs.passband_ripple));
        yline(specs.passband_ripple, 'g--', 'LineWidth', 1.5, 'Label', sprintf('+%.1f dB', specs.passband_ripple));
        
        % Shade passband regions
        ylims_pb = [-3 3];
        fill([freq_ranges.passband_low(1) freq_ranges.passband_low(2) freq_ranges.passband_low(2) freq_ranges.passband_low(1)], ...
             [ylims_pb(1) ylims_pb(1) ylims_pb(2) ylims_pb(2)], 'g', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
        fill([freq_ranges.passband_high(1) freq_ranges.passband_high(2) freq_ranges.passband_high(2) freq_ranges.passband_high(1)], ...
             [ylims_pb(1) ylims_pb(1) ylims_pb(2) ylims_pb(2)], 'g', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
        
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        title(sprintf('Passband Detail - Bounds: [%.2f, %.2f] dB | Ripple Mag: %.2f dB', ...
                      passband_min, passband_max, passband_ripple_magnitude));
        ylim(ylims_pb);
        xlim([0 specs.fs/2]);
        
        subplot(4,2,6);
        plot(f, H_phase_deg, 'b', 'LineWidth', 1.5);
        hold on;
        
        % Shade passband regions
        ylims_phase_pb = ylim;
        fill([freq_ranges.passband_low(1) freq_ranges.passband_low(2) freq_ranges.passband_low(2) freq_ranges.passband_low(1)], ...
             [ylims_phase_pb(1) ylims_phase_pb(1) ylims_phase_pb(2) ylims_phase_pb(2)], 'g', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
        fill([freq_ranges.passband_high(1) freq_ranges.passband_high(2) freq_ranges.passband_high(2) freq_ranges.passband_high(1)], ...
             [ylims_phase_pb(1) ylims_phase_pb(1) ylims_phase_pb(2) ylims_phase_pb(2)], 'g', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
        
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Phase (degrees)');
        title('Phase Response - Passband Detail');
        xlim([0 specs.fs/2]);
        
        % Compute and plot group delay
        if strcmp(filter_type, 'IIR')
            [gd, f_gd] = grpdelay(b, a, nfft, specs.fs);
        else
            % For FIR, group delay is constant (linear phase)
            [gd, f_gd] = grpdelay(b, a, nfft, specs.fs);
        end
        
        subplot(4,2,7);
        plot(f_gd, gd, 'b', 'LineWidth', 1.5);
        hold on;
        
        % Shade stopband region
        ylims_gd = ylim;
        fill([specs.stopband(1) specs.stopband(2) specs.stopband(2) specs.stopband(1)], ...
             [ylims_gd(1) ylims_gd(1) ylims_gd(2) ylims_gd(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Group Delay (samples)');
        title('Group Delay - Full Range');
        xlim([0 specs.fs/2]);
        
        subplot(4,2,8);
        plot(f_gd, gd, 'b', 'LineWidth', 1.5);
        hold on;
        
        % Shade passband regions
        ylims_gd_pb = ylim;
        fill([freq_ranges.passband_low(1) freq_ranges.passband_low(2) freq_ranges.passband_low(2) freq_ranges.passband_low(1)], ...
             [ylims_gd_pb(1) ylims_gd_pb(1) ylims_gd_pb(2) ylims_gd_pb(2)], 'g', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
        fill([freq_ranges.passband_high(1) freq_ranges.passband_high(2) freq_ranges.passband_high(2) freq_ranges.passband_high(1)], ...
             [ylims_gd_pb(1) ylims_gd_pb(1) ylims_gd_pb(2) ylims_gd_pb(2)], 'g', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
        
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Group Delay (samples)');
        title('Group Delay - Passband Detail');
        xlim([0 specs.fs/2]);
        
        % Save figure
        saveas(gcf, sprintf('../outputs/filter_validation_%s.png', strrep(filename, '.mat', '')));
        
    catch ME
        fprintf('  ERROR: Failed to process filter\n');
        fprintf('  Error message: %s\n', ME.message);
        results(idx).filename = filename;
        results(idx).passed = false;
        results(idx).comments = ME.message;
    end
    
    fprintf('\n');
end

% Summary Report
fprintf('\n========== VALIDATION SUMMARY ==========\n');
fprintf('Total filters tested: %d\n', num_filters);
passed_count = sum([results.passed]);
fprintf('Filters passed: %d\n', passed_count);
fprintf('Filters failed: %d\n', num_filters - passed_count);
fprintf('\nFilter Type Summary:\n');
fir_count = sum(strcmp({results.filter_type}, 'FIR'));
iir_count = sum(strcmp({results.filter_type}, 'IIR'));
fprintf('FIR filters: %d\n', fir_count);
fprintf('IIR filters: %d\n', iir_count);
fprintf('\nDetailed Results:\n');
fprintf('%-30s | %-6s | %-10s | %-18s | %-15s\n', 'Filename', 'Type', 'Status', '60Hz Atten (dB)', 'Ripple (dB)');
fprintf('%s\n', repmat('-', 1, 95));
for idx = 1:length(results)
    fprintf('%-30s | %-6s | %-10s | %18.2f | %15.2f\n', ...
            results(idx).filename, ...
            iif(results(idx).passed, 'PASS', 'FAIL'), ...
            results(idx).max_stopband_atten, ...
            results(idx).max_passband_ripple);
end

% Helper function for conditional output
function out = iif(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end
