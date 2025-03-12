function [ild,time] = extractILD(signal,fs,dt)
%   Usage: [ild,time] = extractILD(signal,fs,dt)
%
%   Input parameters:
%       signal: binaural signal (dimensions: time x channel)
%       fs:     sampling rate
%       dt:     time steps in sec
%
%   Output parameters:
%       ild:    broadband ILD per time step in dB
%       time:   time vector in sec

%   #Author: Robert Baumgartner

if size(signal,1) == 2 && size(signal,2) > 2
  signal = signal';
end

dur = size(signal,1)/fs;
windur = round(dt*fs);
win = tukeywin(windur,.01);
time = 0:dt:dur-dt;

ild = zeros(length(time),1);
for tt=1:length(time)
  idt = (tt-1)*windur + (1:windur);
  sig = [win(:),win(:)].*signal(idt,:);
  ild(tt) = mag2db(rms(sig(:,1))) - mag2db(rms(sig(:,2)));
end