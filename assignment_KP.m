%% Programming assignment for Cherish DC6 PhD application
% Petra Kovács

%% Task 1
% Create a Gaussian white noise burst with a duration of 1 sec and smooth
% on/offset ramps of 10 msec. 

% Load HRTF set to obtain sampling rate
database = 'scut';       
HRTFfilename = 'SCUT_KEMAR_radius_all.sofa';
fullfn = fullfile(SOFAdbPath, 'database', database, HRTFfilename);
Obj = SOFAload(fullfn);

sampleFreq = Obj.Data.SamplingRate; 
dur = 1; % 1 sec duration
mu = 0; % mean
sigma = 1; % std
rampLength = 0.01; % 10 msec
numberOfSamples = sampleFreq * dur;

% creating a cosine ramp
numberOfOnsetSamples = sampleFreq * rampLength;
onsetRamp = sin(linspace(0, 1, numberOfOnsetSamples) * pi / 2);
onsetOffsetRamp = [onsetRamp, ones(1, numberOfSamples - 2*numberOfOnsetSamples), fliplr(onsetRamp)];

X = sigma * randn(sampleFreq*dur,1) + mu; % the white noise
signal = X .* onsetOffsetRamp'; % with onset/offset ramps

%% Task 2
% Create 3 spatial trajectories for this noise burst by using the HRTF set
% SCUT_KEMAR_radius_all.sofa (available at
% https://sofacoustics.org/data/database/scut/)

% 1. Horizontal clockwise rotation around the listener at a distance of 1 m
% (varying azimuth, 0° elevation, 1 m radius)
azi = [0 360];
ele = [0 0];
r = [1 1];

outOne = SOFAspat(signal,Obj,azi,ele,r);

% 2. Approach from left (90° azimuth, 0° elevation, radius decreasing from 1 m to 0.2 m) 
azi = [90 90];
ele = [0 0];
r = [1 0.2];

outTwo = SOFAspat(signal,Obj,azi,ele,r);

% 3. Approach from front (0° azimuth, 0° elevation, radius decreasing from 1 m to 0.2 m)
azi = [0 0];
ele = [0 0];
r = [1 0.2];

outThree = SOFAspat(signal,Obj,azi,ele,r);

%% Task 3
% Plot the spectrograms of the generated sounds. Limit the plotted dynamic
% range to 60 dB.
ltfatstart;

subplot(3,1,1);
sgram(outOne(:,1),sampleFreq,'dynrange',60);
title('Horizontal clockwise rotation');

subplot(3,1,2);
sgram(outTwo(:,1),sampleFreq,'dynrange',60);
title('Approach from left');

subplot(3,1,3);
sgram(outThree(:,1),sampleFreq,'dynrange',60);
title('Approach from front');
