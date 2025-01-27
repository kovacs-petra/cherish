function sig = sig_triwave(f0,fs,d,dramp)
%SIG_TRIWAVE - Triangular waveform 
%   Usage: sig = sig_triwave(f0,fs,d,dramp)
%
%   Input parameters:
%     f0     : Fundamental frequency in Hz
%     fs     : Sampling rate in Hz
%     d      : Duration in seconds
%     dramp  : On-/offset ramp (Tukey window) duration in seconds. [Default: 0.01 s]
%
%   Output parameters:
%     sig    : signal wave form

%   #Author: Robert Baumgartner (2025), Acoustics Research Institute, Vienna, Austria

if not(exist('dramp','var'))
  dramp = 0.01;
end

dursamp = round(d*fs); % duration in samples
t = linspace(0,d,dursamp); % Time vector from 0 to d with sampling interval 1/fs

tri_wave = sawtooth(2 * pi * f0 * t, 0.5); % 0.5 for a symmetric triangular wave

win = tukeywin(dursamp,dramp/d);

sig = win'.*tri_wave;

%% Plot the waveform (optional)
% figure;
% plot(t, sig);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Triangular Waveform');
% grid on;

%% Play the sound (optional)
% sound(sig, fs);