originalSignal = stim(:,1);
windowDuration = 0.003; % Window duration in seconds
windowSize = round(fs * windowDuration); % Number of samples in the window: increase if insufficient, decrease if too aggressive smoothing

%% v1 - Gaussian

sigma = windowSize / 6; % Standard deviation of the Gaussian
gaussianWindow = gausswin(windowSize, sigma); % Create Gaussian window
gaussianWindow = gaussianWindow / sum(gaussianWindow); % Normalize
smoothedSignal = conv(stim(:,1), gaussianWindow, 'same');
% smoothedSignal = movmean(stim(:,1), windowSize); % Apply moving average

%% v2 - median filter
alpha = 0.5; % Smoothing factor (0 < alpha < 1)
smoothedSignal = filter(alpha, [1 alpha-1], stim(:,1));

%% v3 - Hanning
hanningWindow = hann(windowSize); % Generate Hanning window
hanningWindow = hanningWindow / sum(hanningWindow); % Normalize the window

paddedSignal = [zeros(windowSize,1); originalSignal; zeros(windowSize,1)];

% Apply the window for smoothing using convolution
smoothedPaddedSignal = filtfilt(hanningWindow, 1, paddedSignal);
smoothedSignal = smoothedPaddedSignal(windowSize+1:end-windowSize);

%% Plot the result of smoothing
time = (0:length(stim)-1) / fs; % Time vector
% plot(time, stim(:,1)); hold on;
plot(smoothedSignal);
% legend('Original Signal', 'Smoothed Signal');
xlabel('Time (s)');
ylabel('Amplitude');
title('Signal Smoothing with Sliding Window');
hold off;

