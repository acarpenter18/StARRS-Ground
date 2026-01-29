% matlab_0210_6.m 
% Note for Yan/Emmanuel: Updates are pending: this is not the final version of the SAR processor, but feel free to have a go with it/make improvements.
clear; close all; clc;

%% Configuration
% Force all figures to dock as tabs in MATLAB window
set(0, 'DefaultFigureWindowStyle', 'docked');

config.data_path = 'C:\Users\amc2918\Desktop\Raw_Data';  % Update this path to where the raw data is located, after offloading from the Samsung SSD
config.timestamp = '20250918_214635'; % Update this to match the flightpath timestamp (timings might not be accurate if RPi date/time is wrong - this happens when the RPi is not connected to the Internet)

% File paths
data_file = fullfile(config.data_path, sprintf('burst_stepped_%s.dat', config.timestamp));
header_file = fullfile(config.data_path, sprintf('headers_%s.dat', config.timestamp));
meta_file = fullfile(config.data_path, sprintf('metadata_%s.txt', config.timestamp));

%% Parse Metadata
fprintf('Parsing metadata...\n');
metadata = parse_metadata_robust(meta_file);

% Extract parameters
radar.fc = metadata.center_freq;
radar.step_bw = metadata.step_bandwidth;
radar.synthetic_bw = metadata.synthetic_bandwidth;
radar.num_steps = metadata.num_freq_steps;
radar.fs = metadata.sample_rate;
radar.burst_duration = metadata.burst_duration;
radar.burst_interval = metadata.burst_interval;
radar.chirps_per_step = metadata.chirps_per_step;
radar.total_chirps_per_burst = metadata.total_chirps_per_burst;
radar.range_bins = metadata.range_bins;
radar.range_start_bin = metadata.range_start_bin;
radar.chirp_duration = 1e-3;  % 1ms chirps
radar.platform_velocity = 1.0; % m/s (set in DJI Ground Station Pro iPad application. May use higher velocity in the future)

% Display configuration
fprintf('\n=== Configuration ===\n');
fprintf('Mode: Burst + Frequency Stepping\n');
fprintf('Center frequency: %.3f GHz\n', radar.fc/1e9);
fprintf('Synthetic bandwidth: %.1f MHz\n', radar.synthetic_bw/1e6);
fprintf('Frequency steps: %d × %.1f MHz\n', radar.num_steps, radar.step_bw/1e6);

%% Load Burst Headers
fprintf('\nLoading burst headers...\n');
fid = fopen(header_file, 'rb');
if fid == -1
    warning('Cannot open header file, proceeding without headers');
    num_bursts = 0;
    burst_headers = [];
else
    header_data = fread(fid, 'uint8');
    fclose(fid);
    
    header_size = 24;
    num_bursts = floor(length(header_data) / header_size);
    
    if mod(length(header_data), header_size) ~= 0
        warning('Header file size not exact multiple of header size');
    end
    
    fprintf('Loaded %d burst headers\n', num_bursts);
    
    burst_headers = struct();
    for i = 1:num_bursts
        offset = (i-1) * header_size;
        burst_headers(i).burst_id = typecast(header_data(offset+1:offset+4), 'uint32');
    end
end

%% Load and Parse Burst Data
fprintf('\nLoading burst data...\n');
fid = fopen(data_file, 'rb');
if fid == -1
    error('Cannot open data file: %s', data_file);
end

fseek(fid, 0, 'eof');
file_size = ftell(fid);
fseek(fid, 0, 'bof');
fprintf('Data file size: %.2f MB\n', file_size/1e6);

% If we don't have num_bursts from headers, estimate from file size
if num_bursts == 0
    chirp_meta_size = 8;
    range_data_size = radar.range_bins * 8;
    bytes_per_chirp = chirp_meta_size + range_data_size;
    bytes_per_burst = radar.total_chirps_per_burst * bytes_per_chirp;
    num_bursts = floor(file_size / bytes_per_burst);
    fprintf('Estimated %d bursts from file size\n', num_bursts);
end

num_bursts = max(1, floor(num_bursts));
all_bursts = cell(num_bursts, 1);

chirp_meta_size = 8;
range_data_size = radar.range_bins * 8;
bytes_per_chirp = chirp_meta_size + range_data_size;
bytes_per_burst = radar.total_chirps_per_burst * bytes_per_chirp;

actual_bursts = 0;
for burst_idx = 1:num_bursts
    remaining = file_size - ftell(fid);
    if remaining < bytes_per_burst
        break;
    end
    
    burst_raw = fread(fid, bytes_per_burst, 'uint8');
    
    burst_data = zeros(radar.range_bins, radar.total_chirps_per_burst, 'single');
    freq_steps = zeros(radar.total_chirps_per_burst, 1);
    center_freqs = zeros(radar.total_chirps_per_burst, 1);
    
    for chirp_idx = 1:radar.total_chirps_per_burst
        offset = (chirp_idx - 1) * bytes_per_chirp;
        
        if offset + bytes_per_chirp > length(burst_raw)
            break;
        end
        
        meta_bytes = burst_raw(offset+1:offset+8);
        if length(meta_bytes) >= 8
            freq_steps(chirp_idx) = double(typecast(uint8(meta_bytes(1:4)), 'uint32'));
            center_freqs(chirp_idx) = double(typecast(uint8(meta_bytes(5:8)), 'single'));
        end
        
        range_start = offset + 9;
        range_end = offset + bytes_per_chirp;
        if range_end <= length(burst_raw)
            range_raw = burst_raw(range_start:range_end);
            if length(range_raw) >= radar.range_bins * 8
                range_floats = typecast(uint8(range_raw), 'single');
                if length(range_floats) >= radar.range_bins * 2
                    real_part = range_floats(1:2:radar.range_bins*2);
                    imag_part = range_floats(2:2:radar.range_bins*2);
                    burst_data(:, chirp_idx) = complex(real_part, imag_part);
                end
            end
        end
    end
    
    all_bursts{burst_idx}.data = burst_data;
    all_bursts{burst_idx}.freq_steps = freq_steps;
    all_bursts{burst_idx}.center_freqs = center_freqs;
    
    actual_bursts = actual_bursts + 1;
end

fclose(fid);

if actual_bursts < num_bursts
    all_bursts = all_bursts(1:actual_bursts);
    num_bursts = actual_bursts;
end

fprintf('Data loading complete: %d bursts loaded\n', num_bursts);

%% Check for Static Test Mode and Data Cropping. Static mode is for lab-based testing of the system when the antennas aren't moving.
is_static = input('Is this static desk data? (y/n): ', 's');
if strcmpi(is_static, 'y')
    fprintf('Static mode: Will show range profiles only (no SAR processing)\n');
    radar.platform_velocity = 0;
    crop_start = 0;
    crop_end = 0;
else
    % Analyze burst power for cropping assistance ... update this in the
    % future. Probably not the best way to crop the data.
    fprintf('\n=== Data Quality Check for Cropping ===\n');
    fprintf('How cropping assistance works:\n');
    fprintf('1. Analyzes power levels in each burst\n');
    fprintf('2. Identifies stable flight region (middle 50%% of data)\n');
    fprintf('3. Detects power variations indicating ascent/descent\n');
    fprintf('4. Suggests crop points where power deviates from stable level\n\n');
    
    burst_powers = zeros(num_bursts, 1);
    for b = 1:num_bursts
        burst_data = all_bursts{b}.data(:);
        burst_data = burst_data(~isnan(burst_data) & burst_data ~= 0);
        if ~isempty(burst_data)
            burst_powers(b) = mean(abs(burst_data).^2);
        end
    end
    
    % Find stable region (middle 50% of data)
    mid_start = round(num_bursts * 0.25);
    mid_end = round(num_bursts * 0.75);
    stable_power = mean(burst_powers(mid_start:mid_end));
    power_std = std(burst_powers(mid_start:mid_end));
    
    % Identify potential ascent/descent regions
    power_threshold_low = stable_power - 2*power_std;
    power_threshold_high = stable_power + 2*power_std;
    
    % Find first stable burst
    ascent_end = 1;
    for b = 1:mid_start
        if burst_powers(b) > power_threshold_low && burst_powers(b) < power_threshold_high
            ascent_end = b;
            break;
        end
    end
    
    % Find last stable burst
    descent_start = num_bursts;
    for b = num_bursts:-1:mid_end
        if burst_powers(b) > power_threshold_low && burst_powers(b) < power_threshold_high
            descent_start = b;
            break;
        end
    end
    
    fprintf('Power analysis results:\n');
    fprintf('  Stable power level: %.2e ± %.2e\n', stable_power, power_std);
    fprintf('  Potential ascent phase: bursts 1-%d (%.1f seconds)\n', ...
            ascent_end, ascent_end * radar.burst_interval);
    fprintf('  Potential descent phase: bursts %d-%d (%.1f seconds)\n', ...
            descent_start, num_bursts, (num_bursts-descent_start) * radar.burst_interval);
    
    % Show power plot for visual inspection
    figure('Name', 'Burst Power Analysis for Cropping');
    plot(1:num_bursts, 10*log10(burst_powers), 'b.-');
    hold on;
    plot([1, num_bursts], [10*log10(stable_power), 10*log10(stable_power)], 'g--', 'LineWidth', 2);
    plot([1, num_bursts], [10*log10(power_threshold_low), 10*log10(power_threshold_low)], 'r--');
    plot([1, num_bursts], [10*log10(power_threshold_high), 10*log10(power_threshold_high)], 'r--');
    plot([ascent_end, ascent_end], ylim, 'r:', 'LineWidth', 2);
    plot([descent_start, descent_start], ylim, 'r:', 'LineWidth', 2);
    xlabel('Burst Number'); ylabel('Power (dB)');
    title('Burst Power - Identify Ascent/Descent for Cropping');
    legend('Burst Power', 'Stable Level', 'Thresholds', '', 'Suggested Crop Points', 'Location', 'best');
    grid on;
    
    % Ask about cropping
    fprintf('\n=== Data Cropping Options ===\n');
    fprintf('Total bursts available: %d\n', num_bursts);
    fprintf('Programmed burst interval: %.1f seconds\n', radar.burst_interval);
    fprintf('Minimum collection time: %.1f seconds\n', num_bursts * radar.burst_interval);
    
    % Check for actual duration in metadata
    actual_duration = 0;
    if isfield(metadata, 'total_duration')
        actual_duration = metadata.total_duration;
        fprintf('Actual collection duration from metadata: %.1f seconds\n', actual_duration);
        actual_burst_interval = actual_duration / num_bursts;
        fprintf('Actual burst interval: %.2f seconds (vs programmed %.2f s)\n', ...
                actual_burst_interval, radar.burst_interval);
    else
        % Estimate actual duration (typical is ~2.25x programmed for current setup at WWS)
        actual_burst_interval = 2.25;
        actual_duration = num_bursts * actual_burst_interval;
        fprintf('Estimated actual duration: %.1f seconds (no metadata found)\n', actual_duration);
        fprintf('Estimated actual burst interval: %.2f seconds\n', actual_burst_interval);
    end
    
    fprintf('\nSuggested cropping based on power analysis:\n');
    fprintf('  Crop first %.1f seconds (ascent)\n', ascent_end * actual_burst_interval);
    fprintf('  Crop last %.1f seconds (descent)\n', (num_bursts-descent_start) * actual_burst_interval);
    
    crop_start_sec = input('Seconds to crop from start (0 for none): ');
    crop_end_sec = input('Seconds to crop from end (0 for none): ');
    
    % Convert seconds to bursts using actual interval
    crop_start = round(crop_start_sec / actual_burst_interval);
    crop_end = round(crop_end_sec / actual_burst_interval);
    
    % Apply cropping
    if crop_start > 0 || crop_end > 0
        original_bursts = num_bursts;
        
        if crop_end > 0
            all_bursts = all_bursts(crop_start+1:end-crop_end);
        else
            all_bursts = all_bursts(crop_start+1:end);
        end
        
        num_bursts = length(all_bursts);
        fprintf('\nCropped data: %d bursts removed from start, %d from end\n', ...
                crop_start, crop_end);
        fprintf('Remaining: %d bursts (%.1f seconds actual time)\n', ...
                num_bursts, num_bursts * actual_burst_interval);
        
        % Update actual duration after cropping
        actual_duration = num_bursts * actual_burst_interval;
    end
    
    % Flight parameters
    fprintf('\n=== Flight Path Estimation ===\n');
    fprintf('After cropping: %d bursts over %.1f seconds actual time\n', ...
            num_bursts, actual_duration);
    
    fprintf('\nOptions:\n');
    fprintf('1. Enter known flight distance (most accurate)\n'); %most accurate, if you know where the crop occured in the flight path. Option 3 is probably best for now.
    fprintf('2. Enter estimated velocity\n');
    fprintf('3. Use default velocity (1 m/s)\n');
    
    option = input('Select option (1/2/3): ', 's');
    
    switch option
        case '1'
            actual_distance = input('Enter flight distance in meters: ');
            if isempty(actual_distance) || actual_distance <= 0
                actual_distance = 135;
                fprintf('Using default: %.1f m\n', actual_distance);
            end
            radar.platform_velocity = actual_distance / actual_duration;
            fprintf('Calculated velocity: %.2f m/s\n', radar.platform_velocity);
            
        case '2'
            radar.platform_velocity = input('Enter velocity in m/s: ');
            if isempty(radar.platform_velocity) || radar.platform_velocity <= 0
                radar.platform_velocity = 1.0;
                fprintf('Using default: %.2f m/s\n', radar.platform_velocity);
            end
            actual_distance = radar.platform_velocity * actual_duration;
            fprintf('Estimated flight distance: %.1f m\n', actual_distance);
            
        otherwise
            radar.platform_velocity = 1.0;
            actual_distance = radar.platform_velocity * actual_duration;
            fprintf('Using default velocity: %.2f m/s\n', radar.platform_velocity);
            fprintf('Estimated flight distance: %.1f m\n', actual_distance);
    end
    
    % Store actual timing for backprojection
    radar.actual_duration = actual_duration;
    radar.actual_burst_interval = actual_duration / num_bursts;
end

%% Calculate Range Parameters
c = 3e8;
range_res_step = c / (2 * radar.step_bw);
range_res_synth = c / (2 * radar.synthetic_bw);

RANGE_CORRECTION_FACTOR = 10; %Range scaling issue fix, due to step to synthetic bandwidth problem.

% fprintf('\n=== Range Parameters ===\n');
% fprintf('Range resolution:\n');
% fprintf('  Per step (%.0f MHz): %.2f m\n', radar.step_bw/1e6, range_res_step);
% fprintf('  Synthetic (%.0f MHz): %.2f m (%.1fx improvement)\n', ...
%         radar.synthetic_bw/1e6, range_res_synth, range_res_step/range_res_synth);
% fprintf('Applying 10x Range Correction\n');

% Calculate range coverage (fixed I think).
actual_range_start = (radar.range_start_bin * range_res_step) / RANGE_CORRECTION_FACTOR; %Start bin is set in the E312 script. Usually 2. The first 2 bins are just saturated returns around the drone.
range_res_effective = range_res_step / RANGE_CORRECTION_FACTOR;

% Create the range axis
synthetic_range_axis = actual_range_start + (0:radar.range_bins-1) * range_res_effective;

% Limit maximum range for realistic display
MAX_DISPLAY_RANGE = 500; % change if necessary. Can go more than 500m, but depends on altitude of the drone and the resultant viewing angle. Schematic in EDF ppt.

fprintf(' Range coverage:\n');
fprintf('  Start range (bin %d): %.1f m\n', radar.range_start_bin, actual_range_start);
fprintf('  Effective range resolution: %.1f m/bin\n', range_res_effective);
fprintf('  Max range in data: %.1f m\n', synthetic_range_axis(end));
fprintf('  Display limited to: %.1f m\n', MAX_DISPLAY_RANGE);

% fprintf('\n=== Understanding Range Axis ===\n');
% fprintf('The range axis does NOT start at 0m.\n');
% fprintf('  - Leftmost data starts at: %.1f m slant range\n', actual_range_start);
% fprintf('  - This is the distance from the radar to the target\n');
% fprintf('  - Early range bins (< %.1fm) excluded to avoid nearfield effects\n', actual_range_start);

%% Process Based on Mode
if strcmpi(is_static, 'y')
    %% Static Mode Analysis
    fprintf('\n=== Static Mode Analysis ===\n');
    
    if num_bursts > 0
        range_axis = synthetic_range_axis(1:min(length(synthetic_range_axis), radar.range_bins));
        range_mask = range_axis <= MAX_DISPLAY_RANGE;
        
        burst1 = all_bursts{1}.data;
        
        figure('Name', 'Static Test - Range Profiles');
        
        subplot(2,2,1);
        avg_profile = mean(abs(burst1), 2);
        plot(range_axis(range_mask), 20*log10(avg_profile(range_mask) + eps));
        xlabel('Range (m)'); ylabel('Power (dB)');
        title('Average Range Profile');
        grid on;
        xlim([range_axis(1), MAX_DISPLAY_RANGE]);
        
        subplot(2,2,2);
        imagesc(1:size(burst1,2), range_axis(range_mask), ...
                20*log10(abs(burst1(range_mask,:)) + eps));
        colormap gray; colorbar;
        xlabel('Chirp Number'); ylabel('Range (m)');
        title('Range Profiles Over Chirps');
        ylim([range_axis(1), MAX_DISPLAY_RANGE]);
        
        subplot(2,2,3);
        hold on;
        colors = lines(radar.num_steps);
        for step = 0:(radar.num_steps-1)
            step_mask_chirps = (all_bursts{1}.freq_steps == step);
            if any(step_mask_chirps)
                step_data = burst1(:, step_mask_chirps);
                step_profile = mean(abs(step_data), 2);
                plot(range_axis(range_mask), 20*log10(step_profile(range_mask) + eps), ...
                     'Color', colors(step+1,:), ...
                     'DisplayName', sprintf('Step %d', step));
            end
        end
        xlabel('Range (m)'); ylabel('Power (dB)');
        title('Range Profile by Frequency Step');
        legend('Location', 'best');
        grid on;
        xlim([range_axis(1), MAX_DISPLAY_RANGE]);
        
        subplot(2,2,4);
        snr_estimate = zeros(num_bursts, 1);
        for b = 1:num_bursts
            burst_data = abs(all_bursts{b}.data(:));
            burst_data = burst_data(~isnan(burst_data));
            if ~isempty(burst_data)
                signal = max(burst_data);
                noise = median(burst_data);
                if noise > 0
                    snr_estimate(b) = 20*log10(signal/noise);
                end
            end
        end
        plot(1:num_bursts, snr_estimate, 'b.-');
        xlabel('Burst Number'); ylabel('SNR (dB)');
        title('Estimated SNR per Burst');
        grid on;
    end
    
else
    %% Dynamic Mode - Full SAR Processing
    fprintf('\n=== SAR Processing Mode ===\n');
    
    % Organize data by frequency steps
    fprintf('Organizing data by frequency steps...\n');
    data_by_step = cell(radar.num_steps, 1);
    for i = 1:radar.num_steps
        data_by_step{i} = [];
    end
    
    for burst_idx = 1:num_bursts
        burst = all_bursts{burst_idx};
        
        for step = 0:(radar.num_steps-1)
            step_mask = (burst.freq_steps == step);
            step_data = burst.data(:, step_mask);
            data_by_step{step+1} = [data_by_step{step+1}, step_data];
        end
    end
    
    for step = 1:radar.num_steps
        fprintf('  Step %d: %d chirps collected\n', step, size(data_by_step{step}, 2));
    end
    
    % Combine frequency steps coherently
    fprintf('\nCombining frequency steps for synthetic bandwidth...\n');
    combined_data = combine_frequency_steps_coherent(data_by_step, radar);
    
    % Find peaks in range profile
    range_profile = mean(abs(combined_data), 2);
    [peaks, peak_locs] = findpeaks(range_profile, 'MinPeakHeight', max(range_profile)*0.3);
    
    if ~isempty(peak_locs)
        fprintf('\nStrong returns detected at (with 10x correction applied):\n');
        for i = 1:min(3, length(peak_locs))
            fprintf('  Peak %d: %.1f m\n', i, synthetic_range_axis(peak_locs(i)));
        end
    end
    
    % Analyze signal strength - improved noise floor estimation
    fprintf('\nAnalyzing signal strength...\n');
    range_profile_avg = mean(abs(combined_data), 2);
    range_profile_db = 20*log10(range_profile_avg + eps);
    
    % Better noise floor: use the lowest 10% of range bins
    sorted_profile = sort(range_profile_avg);
    noise_samples = sorted_profile(1:round(length(sorted_profile)*0.1));
    noise_floor = median(noise_samples);
    noise_floor_db = 20*log10(noise_floor);
    
    fprintf('  Noise floor (from lowest 10%% of bins): %.1f dB\n', noise_floor_db);
    fprintf('  Max signal: %.1f dB\n', max(range_profile_db));
    fprintf('  Peak SNR: %.1f dB\n', max(range_profile_db) - noise_floor_db);
    
    % SAR imaging parameters
    total_flight_distance = radar.platform_velocity * radar.actual_duration;
    
    % Calculate actual resolutions
    wavelength = c / radar.fc;
    sar.scene_center_range = (synthetic_range_axis(1) + MAX_DISPLAY_RANGE) / 2;
    azimuth_res_center = wavelength * sar.scene_center_range / (2 * total_flight_distance);
    
    fprintf('\n=== Resolution Analysis ===\n');
    fprintf('Range resolution: %.2f m\n', range_res_effective);
    fprintf('Azimuth resolution at scene center (%.0fm): %.2f m\n', ...
            sar.scene_center_range, azimuth_res_center);
    
    % Scene parameters
    sar.scene_width = MAX_DISPLAY_RANGE - synthetic_range_axis(1);
    sar.scene_height = total_flight_distance;
    
    % Resolution-based image size
    num_range_pixels = round(sar.scene_width / range_res_effective);
    num_azimuth_pixels = round(sar.scene_height / azimuth_res_center);
    num_range_pixels = max(64, min(2048, num_range_pixels));
    num_azimuth_pixels = max(64, min(2048, num_azimuth_pixels));
    
    fprintf('\n=== SAR Scene Parameters ===\n');
    fprintf('Scene width (range): %.1f m\n', sar.scene_width);
    fprintf('Scene height (azimuth): %.1f m\n', sar.scene_height);
    
    % Azimuth resolution at different ranges
    min_detection_range = synthetic_range_axis(1);
    azimuth_res_min = wavelength * min_detection_range / (2 * total_flight_distance);
    azimuth_res_250m = wavelength * 250 / (2 * total_flight_distance);
    azimuth_res_500m = wavelength * 500 / (2 * total_flight_distance);
    
    fprintf('\n=== Azimuth Resolution vs Range ===\n');
    fprintf('Formula: δ_azimuth = λR / (2L) where λ=wavelength, R=range, L=aperture\n');
    fprintf('  At %.0fm: %.2f m\n', min_detection_range, azimuth_res_min);
    fprintf('  At 250m: %.2f m\n', azimuth_res_250m);
    fprintf('  At 500m: %.2f m\n', azimuth_res_500m);
    fprintf('Note: Azimuth resolution degrades (gets worse) as range increases\n');
    
    % Fast or slow processing
    use_interp = input('\nUse fast processing with interpolation? (y/n, default=y): ', 's');
    if isempty(use_interp) || strcmpi(use_interp, 'y')
        use_interpolation = true;
        fprintf('Using fast processing with interpolation\n');
    else
        use_interpolation = false;
        fprintf('Using full resolution processing (slower but more accurate)\n');
    end
    
    % Generate SAR images for comparison (512x512 grid, and true
    % resolution)
    
    % 1. Original 512x512 image
    fprintf('\n=== Generating 512x512 SAR Image (original) ===\n');
    sar_512 = sar;
    sar_512.image_size = [512, 512];
    [sar_image_512, sar_512] = backprojection_synthetic(combined_data, radar, sar_512, ...
                                                synthetic_range_axis, num_bursts, use_interpolation, MAX_DISPLAY_RANGE);
    
    % 2. Resolution-based image
    fprintf('\n=== Generating Resolution-Based SAR Image ===\n');
    fprintf('Image size: %d (azimuth) × %d (range) pixels\n', num_azimuth_pixels, num_range_pixels);
    fprintf('Pixel spacing:\n');
    fprintf('  Range: %.2f m/pixel (matches resolution: %.2f m)\n', ...
            sar.scene_width/num_range_pixels, range_res_effective);
    fprintf('  Azimuth: %.2f m/pixel (matches resolution: %.2f m)\n', ...
            sar.scene_height/num_azimuth_pixels, azimuth_res_center);
    
    sar_resbased = sar;
    sar_resbased.image_size = [num_azimuth_pixels, num_range_pixels];
    [sar_image_resbased, sar_resbased] = backprojection_synthetic(combined_data, radar, sar_resbased, ...
                                                synthetic_range_axis, num_bursts, use_interpolation, MAX_DISPLAY_RANGE);
    
    % Apply speckle filtering to both
    if max(abs(sar_image_512(:))) > 0
        fprintf('Applying speckle reduction to 512x512 image...\n');
        sar_image_512_filtered = apply_speckle_filter(sar_image_512);
    else
        sar_image_512_filtered = sar_image_512;
    end
    
    if max(abs(sar_image_resbased(:))) > 0
        fprintf('Applying speckle reduction to resolution-based image...\n');
        sar_image_resbased_filtered = apply_speckle_filter(sar_image_resbased);
    else
        sar_image_resbased_filtered = sar_image_resbased;
    end
    
    % Calculate SNR metrics
    burst_snr = zeros(num_bursts, 1);
    for b = 1:num_bursts
        burst_data = abs(all_bursts{b}.data(:));
        burst_data = burst_data(~isnan(burst_data));
        if ~isempty(burst_data)
            signal = max(burst_data);
            noise = median(burst_data);
            if noise > 0
                burst_snr(b) = 20*log10(signal/noise);
            end
        end
    end
    mean_snr = mean(burst_snr(burst_snr > 0));
    
    fprintf('\n=== SNR Analysis ===\n');
    fprintf('Mean burst SNR: %.1f dB (compares max to median within each burst)\n', mean_snr);
    fprintf('Range profile peak SNR: %.1f dB (compares max to noise floor)\n', ...
            max(range_profile_db) - noise_floor_db);
    fprintf('Note: Different SNR metrics can give different values\n');
    
    if mean_snr < 20
        fprintf('\nWARNING: Low SNR detected. Try more aggressive cropping.\n');
    end
    
    %% DISPLAY RESULTS
    % Figure 1 is the data cropping figure.
    % Figure 2: 512x512 SAR Image (Original) Probably get rid of this, as
    % interpolation creates weird curved lines.
    figure('Name', 'SAR Image - 512x512 (Original Interpolated)');
    sar_db_512 = 20*log10(abs(sar_image_512_filtered) + 1e-10);
    valid_pixels = sar_db_512(sar_db_512 > -100 & ~isinf(sar_db_512));
    if ~isempty(valid_pixels)
        max_val = max(valid_pixels);
        min_val = max_val - 40;
    else
        max_val = max(sar_db_512(:));
        min_val = max_val - 40;
    end
    
    h1 = imagesc(sar_512.range_vec, sar_512.azimuth_vec, sar_db_512);
    colormap gray; colorbar;
    xlabel('Range (m)'); ylabel('Azimuth (m)');
    title('SAR Image - 512×512 Grid (Interpolated)');
    caxis([min_val, max_val]);
    axis equal; axis tight;
    xlim([sar_512.range_vec(1), min(sar_512.range_vec(end), MAX_DISPLAY_RANGE)]);
    imcontrast(h1);
    
    % Figure 3: Resolution-Based SAR Image
    figure('Name', 'SAR Image - Resolution-Based Grid');
    sar_db_res = 20*log10(abs(sar_image_resbased_filtered) + 1e-10);
    
    h2 = imagesc(sar_resbased.range_vec, sar_resbased.azimuth_vec, sar_db_res);
    colormap gray; colorbar;
    xlabel('Range (m)'); ylabel('Azimuth (m)');
    title(sprintf('SAR Image - %d×%d Grid (Resolution-Based)', num_azimuth_pixels, num_range_pixels));
    caxis([min_val, max_val]);
    axis equal; axis tight;
    xlim([sar_resbased.range_vec(1), min(sar_resbased.range_vec(end), MAX_DISPLAY_RANGE)]);
    imcontrast(h2);
    
    fprintf('\nInteractive contrast adjustment windows opened.\n');
    
    % Figure 4: Side-by-side comparison
    figure('Name', 'SAR Comparison - 512x512 vs Resolution-Based');
    
    subplot(1,2,1);
    imagesc(sar_512.range_vec, sar_512.azimuth_vec, sar_db_512);
    colormap gray; colorbar;
    xlabel('Range (m)'); ylabel('Azimuth (m)');
    title('512×512 (Interpolated)');
    caxis([min_val, max_val]);
    axis equal; axis tight;
    xlim([sar_512.range_vec(1), min(sar_512.range_vec(end), MAX_DISPLAY_RANGE)]);
    
    subplot(1,2,2);
    imagesc(sar_resbased.range_vec, sar_resbased.azimuth_vec, sar_db_res);
    colormap gray; colorbar;
    xlabel('Range (m)'); ylabel('Azimuth (m)');
    title(sprintf('%d×%d (Resolution-Based)', num_azimuth_pixels, num_range_pixels));
    caxis([min_val, max_val]);
    axis equal; axis tight;
    xlim([sar_resbased.range_vec(1), min(sar_resbased.range_vec(end), MAX_DISPLAY_RANGE)]);
    
    % Figure 5: Range Analysis with improved noise floor
    figure('Name', 'Range Analysis');
    
    range_mask_display = synthetic_range_axis <= MAX_DISPLAY_RANGE;
    
    subplot(1,2,1);
    plot(synthetic_range_axis(range_mask_display), range_profile_db(range_mask_display), 'b-', 'LineWidth', 1.5);
    hold on;
    plot([synthetic_range_axis(1), MAX_DISPLAY_RANGE], ...
         [noise_floor_db, noise_floor_db], 'r--', 'DisplayName', 'Noise Floor', 'LineWidth', 2);
    plot([synthetic_range_axis(1), MAX_DISPLAY_RANGE], ...
         [noise_floor_db+10, noise_floor_db+10], 'g--', 'DisplayName', '+10 dB above noise', 'LineWidth', 1.5);
    xlabel('Range (m)'); ylabel('Power (dB)');
    title('Range Profile');
    legend('Location', 'best');
    grid on;
    xlim([synthetic_range_axis(1), MAX_DISPLAY_RANGE]);
    
    subplot(1,2,2);
    near_range = min(250, MAX_DISPLAY_RANGE);
    near_mask = synthetic_range_axis <= near_range;
    plot(synthetic_range_axis(near_mask), range_profile_db(near_mask), 'b-', 'LineWidth', 1.5);
    hold on;
    plot([synthetic_range_axis(1), near_range], ...
         [noise_floor_db, noise_floor_db], 'r--', 'DisplayName', 'Noise Floor', 'LineWidth', 2);
    xlabel('Range (m)'); ylabel('Power (dB)');
    title('Near Range Profile (<250m)');
    legend('Location', 'best');
    grid on;
    xlim([synthetic_range_axis(1), near_range]);
    
    % Figure 6: Resolution Comparison
    figure('Name', 'Resolution Verification');
    
    single_step_profile = mean(abs(data_by_step{1}), 2);
    synthetic_profile = mean(abs(combined_data), 2);
    
    valid_range = min(length(single_step_profile), length(synthetic_range_axis));
    range_mask_comp = synthetic_range_axis(1:valid_range) <= MAX_DISPLAY_RANGE;
    
    plot(synthetic_range_axis(range_mask_comp), ...
         20*log10(single_step_profile(range_mask_comp) + eps), 'b-', ...
         'DisplayName', 'Single Step', 'LineWidth', 1.5);
    hold on;
    plot(synthetic_range_axis(range_mask_comp), ...
         20*log10(synthetic_profile(range_mask_comp) + eps), 'r-', ...
         'DisplayName', 'Synthetic BW', 'LineWidth', 1.5);
    xlabel('Range (m)'); ylabel('Power (dB)');
    title('Resolution Comparison');
    legend('Location', 'best');
    grid on;
    xlim([synthetic_range_axis(1), MAX_DISPLAY_RANGE]);
    
    % Figure 7: SNR Analysis
    figure('Name', 'SNR Analysis');
    
    subplot(2,1,1);
    plot(1:num_bursts, burst_snr, 'b.-', 'LineWidth', 1.5);
    hold on;
    plot([1, num_bursts], [mean_snr, mean_snr], 'r--', ...
         'DisplayName', sprintf('Mean: %.1f dB', mean_snr), 'LineWidth', 1.5);
    xlabel('Burst Number'); ylabel('SNR (dB)');
    title(sprintf('Burst SNR (Mean: %.1f dB)', mean_snr));
    grid on;
    legend('Location', 'best');
    
    subplot(2,1,2);
    histogram(burst_snr, 20, 'FaceColor', [0.5 0.5 0.5]);
    xlabel('SNR (dB)'); ylabel('Count');
    title('SNR Distribution');
    grid on;
    
    % Figure 8: Waterfall Display
    figure('Name', 'Range-Compressed Data');
    
    valid_idx = synthetic_range_axis <= MAX_DISPLAY_RANGE;
    imagesc(1:size(combined_data,2), synthetic_range_axis(valid_idx), ...
            20*log10(abs(combined_data(valid_idx,:)) + eps));
    colormap gray; colorbar;
    xlabel('Chirp Number'); ylabel('Range (m)');
    title('Range-Compressed Data');
    caxis([max(range_profile_db)-50, max(range_profile_db)]);
    ylim([synthetic_range_axis(1), MAX_DISPLAY_RANGE]);
    
    % Figure 9: Platform Trajectory
    figure('Name', 'Platform Trajectory');
    
    n_chirps_total = size(combined_data, 2);
    chirps_per_burst_combined = round(n_chirps_total / num_bursts);
    
    platform_pos_plot = zeros(n_chirps_total, 1);
    chirp_times = zeros(n_chirps_total, 1);
    
    for chirp_idx = 1:n_chirps_total
        burst_num = floor((chirp_idx - 1) / chirps_per_burst_combined) + 1;
        chirp_in_burst = mod(chirp_idx - 1, chirps_per_burst_combined) + 1;
        
        if burst_num <= num_bursts
            burst_start_time = (burst_num - 1) * radar.actual_burst_interval;
            chirp_time = burst_start_time + (chirp_in_burst - 1) * radar.chirp_duration;
            
            chirp_times(chirp_idx) = chirp_time;
            platform_pos_plot(chirp_idx) = radar.platform_velocity * chirp_time;
        end
    end
    
    subplot(2,1,1);
    plot(1:n_chirps_total, platform_pos_plot, 'b-', 'LineWidth', 1.5);
    xlabel('Chirp Number'); ylabel('Platform Position (m)');
    title(sprintf('Platform Trajectory - Burst Mode (v=%.2f m/s)', radar.platform_velocity));
    grid on;
    % text(0.5, 0.95, 'Stepped pattern is NORMAL for burst mode:', ...
    %      'Units', 'normalized', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
    %text(0.5, 0.88, '• Small steps = motion during burst', ...
         %'Units', 'normalized', 'VerticalAlignment', 'top');
    %text(0.5, 0.83, '• Large jumps = motion during gap', ...
         %'Units', 'normalized', 'VerticalAlignment', 'top');
    
    subplot(2,1,2);
    plot(chirp_times, platform_pos_plot, 'b-', 'LineWidth', 1.5);
    hold on;
    plot([0, chirp_times(end)], [0, platform_pos_plot(end)], 'r--', ...
         'DisplayName', 'Expected (constant velocity)', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Platform Position (m)');
    title('Position vs Time (Should be Linear for Constant Velocity)');
    legend('Actual', 'Constant velocity reference', 'Location', 'northwest');
    grid on;
    
    % Figure 9: Near Range Detail
    figure('Name', 'SAR Near Range Detail (Resolution-Based)');
    near_range_limit = min(100, sar_resbased.range_vec(end));
    range_mask = sar_resbased.range_vec >= sar_resbased.range_vec(1) & sar_resbased.range_vec <= near_range_limit;
    
    if any(range_mask)
        imagesc(sar_resbased.range_vec(range_mask), sar_resbased.azimuth_vec, sar_db_res(:,range_mask));
        colormap gray; colorbar;
        axis image;
        caxis([min_val, max_val]);
        xlabel('Range (m)'); ylabel('Azimuth (m)');
        title(sprintf('Near Range Detail (%.0f-%.0fm)', sar_resbased.range_vec(1), near_range_limit));
        xlim([sar_resbased.range_vec(1), near_range_limit]);
    end
    
    % Figure 10: Mid Range Detail
    figure('Name', 'SAR Mid Range Detail (Resolution-Based)');
    mid_range_mask = sar_resbased.range_vec >= 100 & sar_resbased.range_vec <= min(300, sar_resbased.range_vec(end));
    
    if any(mid_range_mask)
        imagesc(sar_resbased.range_vec(mid_range_mask), sar_resbased.azimuth_vec, sar_db_res(:,mid_range_mask));
        colormap gray; colorbar;
        axis image;
        caxis([min_val, max_val]);
        xlabel('Range (m)'); ylabel('Azimuth (m)');
        title('Mid Range Detail (100-300m)');
        xlim([100, min(300, sar_resbased.range_vec(end))]);
    end
    
    % Figure 11: Ground Range SAR Image
    drone_altitude = input('\nEnter drone altitude for ground range conversion (m, 0 to skip): ');
    
    if ~isempty(drone_altitude) && drone_altitude > 0
        fprintf('\n=== Ground Range Conversion ===\n');
        fprintf('Drone altitude: %.1f m\n', drone_altitude);
        
        % Convert slant to ground range
        slant_range_vec = sar_resbased.range_vec;
        ground_range_vec = sqrt(max(0, slant_range_vec.^2 - drone_altitude^2));
        
        % Find valid ranges (slant >= altitude)
        valid_mask = slant_range_vec >= drone_altitude;
        
        fprintf('Valid ground range data: %.1f%% (%d/%d range bins)\n', ...
                100*sum(valid_mask)/length(valid_mask), sum(valid_mask), length(valid_mask));
        
        if sum(valid_mask) < 0.5*length(valid_mask)
            fprintf('WARNING: Over 50%% of range bins have slant < altitude\n');
            fprintf('These ranges cannot represent ground targets\n');
        end
        
        % Figure 12: SAR Image with Ground Range X-Axis
        figure('Name', 'SAR Image Ground Range');
        imagesc(ground_range_vec, sar_resbased.azimuth_vec, sar_db_res);
        colormap gray; colorbar;
        xlabel('Ground Range (m)'); ylabel('Azimuth (m)');
        title(sprintf('SAR Image - Ground Range (Altitude=%.0fm)', drone_altitude));
        caxis([min_val, max_val]);
        axis equal; axis tight;
        
        % Mark boundary between invalid/valid ground range if needed
        hold on;
        if ~all(valid_mask)
            % Find the ground range corresponding to altitude (where slant = altitude)
            boundary_ground = 0;  % At slant=altitude, ground range = 0
            plot([boundary_ground, boundary_ground], ylim, 'r--', 'LineWidth', 2);
            text(10, max(sar_resbased.azimuth_vec)*0.95, ...
                 'Valid ground range →', 'Color', 'red', 'FontWeight', 'bold', 'FontSize', 10);
        end
    end

    % Need to update figures to show true beam spreading, and pixel size
    % change with range.
end
%% Save Results
output_file = fullfile(config.data_path, sprintf('processed_%s.mat', config.timestamp));
if strcmpi(is_static, 'y')
    save(output_file, 'all_bursts', 'radar', 'synthetic_range_axis');
else
    if exist('sar_image_512', 'var') && exist('sar_image_resbased', 'var')
        save(output_file, 'sar_image_512', 'sar_image_512_filtered', 'sar_512', ...
             'sar_image_resbased', 'sar_image_resbased_filtered', 'sar_resbased', ...
             'combined_data', 'all_bursts', 'radar', 'synthetic_range_axis');
    else
        save(output_file, 'combined_data', 'all_bursts', 'radar', 'synthetic_range_axis');
        warning('SAR images were not created');
    end
end
fprintf('\nResults saved to: %s\n', output_file);

%% Functions

function metadata = parse_metadata_robust(filename)
    metadata = struct();
    
    % Set defaults
    metadata.center_freq = 5.4e9;
    metadata.step_bandwidth = 2e6;
    metadata.synthetic_bandwidth = 10e6;
    metadata.num_freq_steps = 5;
    metadata.sample_rate = 2e6;
    metadata.burst_duration = 0.1;
    metadata.burst_interval = 1.0;
    metadata.chirps_per_step = 20;
    metadata.total_chirps_per_burst = 100;
    metadata.range_bins = 256;
    metadata.range_start_bin = 10;
    
    if ~exist(filename, 'file')
        warning('Metadata file not found, using defaults');
        return;
    end
    
    fid = fopen(filename, 'r');
    while ~feof(fid)
        line = fgetl(fid);
        
        if ~ischar(line) || isempty(line)
            continue;
        end
        
        if contains(line, '=')
            parts = split(line, '=');
            
            if length(parts) >= 2
                key = strtrim(parts{1});
                value = strtrim(parts{2});
                
                if isempty(key) || isempty(value)
                    continue;
                end
                
                key = matlab.lang.makeValidName(key);
                
                num_value = str2double(value);
                if ~isnan(num_value)
                    metadata.(key) = num_value;
                else
                    metadata.(key) = value;
                end
            end
        end
    end
    fclose(fid);
end

function combined = combine_frequency_steps_coherent(data_by_step, radar)
    num_steps = length(data_by_step);
    max_chirps = 0;
    
    for i = 1:num_steps
        if ~isempty(data_by_step{i})
            max_chirps = max(max_chirps, size(data_by_step{i}, 2));
        end
    end
    
    fprintf('  Coherent frequency combination:\n');
    fprintf('  Original bins: %d\n', radar.range_bins);
    
    combined = zeros(radar.range_bins, max_chirps);
    
    for step = 1:num_steps
        if isempty(data_by_step{step})
            continue;
        end
        
        step_data = data_by_step{step};
        n_chirps = size(step_data, 2);
        
        for chirp_idx = 1:min(n_chirps, max_chirps)
            if chirp_idx <= size(step_data, 2)
                combined(:, chirp_idx) = combined(:, chirp_idx) + step_data(:, chirp_idx);
            end
        end
    end
    
    combined = combined / num_steps;
    
    fprintf('  NOTE: Power reduction of ~%.1f dB is expected due to averaging\n', ...
            20*log10(num_steps));
    
    if any(isnan(combined(:)))
        fprintf('    Warning: NaN values detected, setting to zero\n');
        combined(isnan(combined)) = 0;
    end
end

function [img, sar] = backprojection_synthetic(data, radar, sar, range_axis, num_bursts, use_interpolation, max_range)
    [n_range, n_azimuth] = size(data);
    
    fprintf('  Input data: %d range bins x %d azimuth samples\n', n_range, n_azimuth);
    
    range_vec = linspace(range_axis(1), min(max_range, range_axis(end)), sar.image_size(2));
    azimuth_vec = linspace(0, sar.scene_height, sar.image_size(1));
    
    fprintf('  Processing range: %.1f to %.1f m\n', range_vec(1), range_vec(end));
    fprintf('  Image dimensions: %d x %d pixels\n', sar.image_size(1), sar.image_size(2));
    
    sar.range_vec = range_vec;
    sar.azimuth_vec = azimuth_vec;
    
    % Platform positions - linear motion
    if isfield(radar, 'actual_burst_interval')
        actual_interval = radar.actual_burst_interval;
    else
        actual_interval = radar.burst_interval;
    end
    
    platform_pos = zeros(n_azimuth, 1);
    chirps_per_burst_combined = round(n_azimuth / num_bursts);
    
    fprintf('  Chirps per burst in combined data: %d (original was %d total, %d per step)\n', ...
            chirps_per_burst_combined, radar.total_chirps_per_burst, radar.chirps_per_step);
    
    for chirp_idx = 1:n_azimuth
        burst_num = floor((chirp_idx - 1) / chirps_per_burst_combined) + 1;
        chirp_in_burst = mod(chirp_idx - 1, chirps_per_burst_combined) + 1;
        
        if burst_num <= num_bursts
            burst_start_time = (burst_num - 1) * actual_interval;
            chirp_time = burst_start_time + (chirp_in_burst - 1) * radar.chirp_duration;
            platform_pos(chirp_idx) = radar.platform_velocity * chirp_time;
        end
    end
    
    fprintf('  Platform trajectory: %.1f to %.1f m (linear)\n', platform_pos(1), platform_pos(end));
    
    % Initialise image
    img = zeros(sar.image_size);
    
    if use_interpolation
        fprintf('  Fast backprojection with interpolation...\n');
        fprintf('  Progress: ');
        
        step_az = 4;
        step_rg = 4;
        step_pulse = 10;
        
        for az_idx = 1:step_az:sar.image_size(1)
            if mod(az_idx, 40) == 0
                fprintf('.');
            end
            
            target_y = azimuth_vec(az_idx);
            
            for rg_idx = 1:step_rg:sar.image_size(2)
                target_x = range_vec(rg_idx);
                
                pixel_sum = 0;
                valid_samples = 0;
                
                for pulse_idx = 1:step_pulse:n_azimuth
                    if pulse_idx > n_azimuth
                        break;
                    end
                    
                    platform_y = platform_pos(pulse_idx);
                    R = sqrt(target_x^2 + (target_y - platform_y)^2);
                    
                    [min_diff, bin_idx] = min(abs(range_axis - R));
                    
                    if bin_idx >= 1 && bin_idx <= n_range && min_diff < (range_axis(2)-range_axis(1))
                        data_val = data(bin_idx, pulse_idx);
                        
                        if ~isnan(data_val) && data_val ~= 0
                            phase_corr = exp(-1j * 4 * pi * radar.fc * R / 3e8);
                            pixel_sum = pixel_sum + data_val * phase_corr;
                            valid_samples = valid_samples + 1;
                        end
                    end
                end
                
                if valid_samples > 0
                    img(az_idx, rg_idx) = pixel_sum / sqrt(valid_samples);
                end
            end
        end
        
        % Interpolate
        [X, Y] = meshgrid(1:step_rg:sar.image_size(2), 1:step_az:sar.image_size(1));
        [Xq, Yq] = meshgrid(1:sar.image_size(2), 1:sar.image_size(1));
        
        img_sparse = img(1:step_az:end, 1:step_rg:end);
        if numel(img_sparse) > 4 && sum(abs(img_sparse(:))) > 0
            img = interp2(X, Y, img_sparse, Xq, Yq, 'linear', 0);
        end
        
        fprintf(' Done!\n');
        
    else
        fprintf('  Full resolution backprojection...\n');
        fprintf('  Progress: ');
        
        for az_idx = 1:sar.image_size(1)
            if mod(az_idx, 50) == 0
                fprintf('%d%%..', round(100*az_idx/sar.image_size(1)));
            end
            
            target_y = azimuth_vec(az_idx);
            
            for rg_idx = 1:sar.image_size(2)
                target_x = range_vec(rg_idx);
                
                pixel_sum = 0;
                valid_samples = 0;
                
                for pulse_idx = 1:2:n_azimuth
                    platform_y = platform_pos(pulse_idx);
                    R = sqrt(target_x^2 + (target_y - platform_y)^2);
                    
                    [min_diff, bin_idx] = min(abs(range_axis - R));
                    
                    if bin_idx >= 1 && bin_idx <= n_range && min_diff < (range_axis(2)-range_axis(1))/2
                        data_val = data(bin_idx, pulse_idx);
                        
                        if ~isnan(data_val) && data_val ~= 0
                            phase_corr = exp(-1j * 4 * pi * radar.fc * R / 3e8);
                            pixel_sum = pixel_sum + data_val * phase_corr;
                            valid_samples = valid_samples + 1;
                        end
                    end
                end
                
                if valid_samples > 0
                    img(az_idx, rg_idx) = pixel_sum / sqrt(valid_samples);
                end
            end
        end
        
        fprintf(' Done!\n');
    end
end

function filtered_img = apply_speckle_filter(img)
    if max(abs(img(:))) == 0
        filtered_img = img;
        return;
    end
    
    intensity = abs(img).^2;
    intensity_filt = medfilt2(intensity, [2, 2]);
    
    window_size = 3;
    h = ones(window_size) / window_size^2;
    mean_img = imfilter(intensity_filt, h, 'replicate');
    
    weight = 0.5;
    intensity_lee = (1-weight) * intensity_filt + weight * mean_img;
    
    phase_orig = angle(img);
    filtered_img = sqrt(max(intensity_lee, 0)) .* exp(1j * phase_orig);
    
    if max(abs(filtered_img(:))) == 0
        filtered_img = img;
    end
    
    fprintf('  Gentle speckle reduction applied\n');
end
