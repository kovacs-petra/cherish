%% Test 1: change reverb path settings only and compare
u = pnet('udpsocket',10017); % Listen port in BRT
oscsend(u, '/control/connect', 'si', 'localhost',10019);

oscsend(u,'/removeAllSources', 's', '');

SWname = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\20-May-2025-1333_400\20-May-2025-1333_400-1.wav';
speechname = 'C:\Program Files\University of Malaga\BeRTA-Renderer App\data\resources\MusArch_Sample_48kHz_Anechoic_FemaleSpeech.wav';

oscsend(u, '/source/loadSource', 'sss', 'SW', SWname, 'DirectivityModel');
oscsend(u, '/source/loadSource', 'sss', 'speech', speechname, 'DirectivityModel');

oscsend(u, '/enableModel', 'sB', 'DirectPath',1);
oscsend(u, '/enableModel', 'sB', 'ReverbPath',0);
oscsend(u, '/enableModel', 'sB', 'ReverbPath',1);

oscsend(u, '/source/play', 's', 'SW');
oscsend(u, '/source/play', 's', 'speech');

speechR = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\speechR.mat';
speechNR = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\speechNR.mat';
SWR = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\SWR.mat';
SWNR = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\SWNR.mat';

oscsend(u, '/source/playAndRecord', 'sssf', 'speech', speechR, 'mat', 10); 
oscsend(u, '/source/playAndRecord', 'sssf', 'speech', speechNR, 'mat', 10); 
oscsend(u, '/source/playAndRecord', 'sssf', 'SW', SWR, 'mat', 10); 
oscsend(u, '/source/playAndRecord', 'sssf', 'SW', SWNR, 'mat', 10); 

fs = 48000;
soundSWNR = Data.Receiver';
soundSWR = Data.Receiver';
soundsc(soundSWNR,fs)
soundsc(soundSWR,fs)

audiowrite([pwd, '\speechR.wav'], Data.Receiver',Data.SamplingRate);
audiowrite([pwd, '\speechNR.wav'], Data.Receiver',Data.SamplingRate);
audiowrite([pwd, '\SWR.wav'], Data.Receiver',Data.SamplingRate);
audiowrite([pwd, '\SWNR.wav'], Data.Receiver',Data.SamplingRate);

%% Test 2: change all settings I use in stimulus generation, including reverb path, and compare
u = pnet('udpsocket',10017); % Listen port in BRT
oscsend(u, '/control/connect', 'si', 'localhost',10019);
oscsend(u, '/listener/enableSpatialization', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableInterpolation', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableNearFieldEffect', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableITD', 'sB', 'DefaultListener',1);
oscsend(u, '/enableModel', 'sB', 'FreeField',0);
%%%
oscsend(u, '/enableModel', 'sB', 'ReverbPath',0);
oscsend(u, '/environment/enableReverbPath', 'sB', 'SDN',0);
%%%
oscsend(u, '/environment/enableDirectPath', 'sB', 'SDN',1);
oscsend(u, '/listener/enableParallaxCorrection', 'sB', 'DefaultListener',1);
oscsend(u, '/enableModel', 'sB', 'DirectPath', 1);
% oscsend(u, '/environment/enableDistanceAttenuation', 'sB', 'FreeField', 1);
oscsend(u, '/environment/enableDistanceAttenuation', 'sB', 'SDN', 1);

oscsend(u,'/removeAllSources', 's', '');
oscsend(u, '/source/loadSource', 'sss', 'SW', SWname, 'DirectivityModel');
oscsend(u, '/source/play', 's', 'SW');

SWR2 = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\SWR2.mat';
SWNR2 = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\SWNR2.mat';

oscsend(u, '/source/playAndRecord', 'sssf', 'SW', SWR2, 'mat', 10); 
oscsend(u, '/source/playAndRecord', 'sssf', 'SW', SWNR2, 'mat', 10); 

soundSWNR2 = Data.Receiver';
soundSWR2 = Data.Receiver';
soundsc(soundSWNR2,fs)
soundsc(soundSWR2,fs)

% Compare reverb stimulus to original stimulus I used in the pilot
origSW = stimStruct.stim{1,1,1}{1};

soundsc(origSW,fs)

%% Test 3: Stimuli with amplitude modulation
t = linspace(0,2,2*fs); % Time vector from 0 to 2 s with sampling interval 1/fs
triwave = sawtooth(2 * pi * 400 * t, 0.5);
AM40 = sin(2*pi*40*t); % 40 Hz sine wave
AM2 = sin(2*pi*2*t); % 2 Hz sine wave
triwave40 = triwave.*AM40; 
triwave2 = triwave.*AM2; 
win = tukeywin(size(triwave,2),(0.1*fs)/size(triwave,2)); % 0.1 s on-offset ramp
triwave40 = triwave40.*win'; audiowrite("triwave_AM40Hz.wav", triwave40,fs);
triwave2 = triwave2.*win'; audiowrite("triwave_AM2Hz.wav", triwave2,fs);

tw40path = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\triwave_AM40Hz.wav';
tw2path = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\triwave_AM2Hz.wav';

oscsend(u, '/source/loadSource', 'sss', 'tw40', tw40path, 'DirectivityModel');
oscsend(u, '/source/loadSource', 'sss', 'tw2', tw2path, 'DirectivityModel');

oscsend(u, '/enableModel', 'sB', 'ReverbPath',0);
oscsend(u, '/enableModel', 'sB', 'ReverbPath',1);

oscsend(u, '/source/play', 's', 'tw40');
oscsend(u, '/source/play', 's', 'tw2');

tw40mat = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\tw40.mat';
tw2mat = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\tw2.mat';
tw40NR = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\tw40NR.mat';
tw2NR = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\tw2NR.mat';
oscsend(u, '/source/playAndRecord', 'sssf', 'tw40', tw40NR, 'mat', 10); 
oscsend(u, '/source/playAndRecord', 'sssf', 'tw2', tw2mat, 'mat', 10); 

%% Wall absorption thingies
oscsend(u, '/environment/setShoeBoxRoom', 'sfff', 'SDN', 30,40,30);
oscsend(u, '/environment/setWallAbsorption', 'sif', 'SDN', 5, 0.1);
oscsend(u, '/source/play', 's', 'source');

%% Near AND far sounds with reverb
u = pnet('udpsocket',10017); % Listen port in BRT
oscsend(u, '/control/connect', 'si', 'localhost',10019);
oscsend(u,'/removeAllSources', 's', '');

SWname = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration\externalization_pilot\21-May-2025-151_400\21-May-2025-151_400-1.wav';
oscsend(u, '/source/loadSource', 'sss', 'SW', SWname, 'DirectivityModel');

