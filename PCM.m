clc; clear; close all;

% Parameters
f = 100;                  % Analog signal frequency
A = 3;                    % Amplitude ±3V (6V peak-to-peak)
fs = 2.4 * f;             % Sampling frequency (20% above Nyquist)
T = 1;                    % Total time in seconds
t = 0:1e-4:T;             % Time vector for continuous signal
ts = 0:1/fs:T;            % Time vector for sampled signal

% Original Analog Signal
x = A * sin(2*pi*f*t);    % Continuous-time signal
xs = A * sin(2*pi*f*ts);  % Sampled signal

% --- Quantization ---
q_levels = linspace(-A, A, 8);   % 8 uniform levels from -3V to 3V
num_samples = length(xs);
xq = zeros(1, num_samples);     % Quantized signal
q_idx = zeros(1, num_samples);  % Level index for PCM

for i = 1:num_samples
    min_diff = abs(xs(i) - q_levels(1));
    index = 1;
    for j = 2:length(q_levels)
        this_diff = abs(xs(i) - q_levels(j));
        if this_diff < min_diff
            min_diff = this_diff;
            index = j;
        end
    end
    xq(i) = q_levels(index);  % Save quantized value
    q_idx(i) = index - 1;     % Save index for binary (0 to 7)
end

% --- Reconstruction using Zero-Order Hold and LPF ---
xr = zeros(1, length(t));
for i = 1:length(ts)-1
    xr(t >= ts(i) & t < ts(i+1)) = xq(i);
end
xr(t >= ts(end)) = xq(end);

[b, a] = butter(4, f/(fs/2));     % LPF design (cutoff 100 Hz)
x_rec = filter(b, a, xr);        % Reconstructed signal

% --- Plot signals ---
figure;
subplot(3,1,1); plot(t, x); title('Original Signal'); grid on;
subplot(3,1,2); stem(ts, xq, 'filled'); title('Quantized Samples'); grid on;
subplot(3,1,3); plot(t, x_rec); title('Reconstructed Signal'); grid on;

% --- PCM Encoding ---
% Convert indices to 3-bit binary
bitstream = [];
for i = 1:num_samples
    bits = dec2bin(q_idx(i), 3) - '0';  % Convert to binary
    bitstream = [bitstream, bits];     % Append to stream
end

% --- Line Coding (Manchester Encoding) ---
% Each bit: 1 -> [1 0], 0 -> [0 1]
manchester = zeros(1, 2*length(bitstream));
for i = 1:length(bitstream)
    if bitstream(i) == 1
        manchester(2*i-1:2*i) = [1 0];
    else
        manchester(2*i-1:2*i) = [0 1];
    end
end

% Plot zoomed-in portion of Manchester encoded signal
figure;

% Let's show the first 60 Manchester bits (i.e., 30 original bits)
N = 60;
if length(manchester) < N
    N = length(manchester);
end

stairs(0:N-1, manchester(1:N), 'b', 'LineWidth', 2);
ylim([-0.2 1.2]);
xlabel('Time (bit steps)');
ylabel('Voltage');
title('Manchester Line Code (Zoomed View)');
grid on;

% --- Bit Rate Calculation ---
%bits_per_sample = 3;             % 8 levels = 3 bits
%theoretical_bitrate = fs * bits_per_sample;

% Simulated bitrate:
simulated_bitrate = length(bitstream) / T;
%fprintf('Theoretical Bit Rate: %.2f bps\n', theoretical_bitrate);

fprintf('Simulated Bit Rate:   %.2f bps\n', simulated_bitrate);

% --- TDM with 5 signals ---
% Multiplex 5 bitstreams (identical for simulation)
tdm_stream = [];
for i = 1:length(bitstream)
    for ch = 1:5
        tdm_stream = [tdm_stream bitstream(i)];
    end
end


% Plot Zoomed-in TDM Multiplexed Signal (first 200 bits)
figure;
N = 200;  % Show only the first 200 bits
if length(tdm_stream) < N
    N = length(tdm_stream);
end

stairs(0:N-1, tdm_stream(1:N), 'b', 'LineWidth', 2);

ylim([-0.2 1.2]);
xlabel('Time (bit steps)');
ylabel('Voltage');
title('TDM Multiplexed Signal (Zoomed View, 5 Channels)');
grid on;

% TDM Bitrate
tdm_bitrate = simulated_bitrate * 5;
fprintf('TDM Bit Rate (5 signals): %.2f bps\n', tdm_bitrate);
