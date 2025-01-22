function cueSpat = BRTspat(cue, totalDur, azi, ele, r, savename, updateRate, u)

% Spatialize stimuli by sending OSC commands to the BRT App.
% Required:
% - oscsend.m
% - BRT Renderer Application 
%
% Inputs:
% cue   	    - path for the .wav file to be spatialized
% totalDur      - duration of the sound stimulus (cue) in s
% azi           - azimuth trajectory vector 
% ele           - elevation trajectory vector
% r             - distance trajectory vector 
% savename      - name of the .mat files to save the cue into
% updateRate    - at what rate should the location be updated (ms)
% u             - udpport to use for the communication with BRT

%% Load and record cue
oscsend(u, '/removeAllSources', 's', '');
oscsend(u, '/source/loadSource', 'sss', 'Cue', cue, 'DirectivityModel'); 

% Set trajectory
[x,y,z] = sph2cart(deg2rad(azi),deg2rad(ele),r);
steps = round(totalDur/updateRate);

% Set start location
start = 1;
oscsend(u, '/source/location', 'sfff', 'Cue', x(start), y(start), z(start));
pause(1); % 1 s delay

% Start recording, then rapidly change location in a for loop
oscsend(u, '/source/playAndRecord', 'sssf', 'Cue', savename, 'mat', totalDur);

startTime = tic;
for i = start+1 : steps
    oscsend(u, '/source/location', 'sfff', 'Cue', x(i), y(i), z(i));
    elapsedTime = toc(startTime);
    targetTime = i * updateRate;
    if elapsedTime < targetTime
        pause(targetTime - elapsedTime);
    end
end
pause(1);
oscsend(u, '/stop', 's', '');

%% Create output
load(savename,"Data");
cueSpat = Data.Receiver;


