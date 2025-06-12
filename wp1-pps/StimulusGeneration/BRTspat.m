function [cueSpat, cueParams] = BRTspat(cue, totalDur, azi, ele, r, savename, updateRate, u)
% Usage: [cueSpat, cueParams] = BRTspat(cue, totalDur, azi, ele, r, savename, updateRate, u)
%
% Spatialize stimuli by sending OSC commands to the BRT Renderer App. To be called
% by the main stimulus generation script (generateStimuli.m), usually.
%
% Inputs:
% cue   	    - path for the .wav file to be spatialized
% totalDur      - duration of the sound stimulus (cue) in s
% azi           - azimuth trajectory vector 
% ele           - elevation trajectory vector
% r             - distance trajectory vector 
% savename      - name of the .mat file to save the cue into
% updateRate    - at what rate should the location be updated (ms)
% u             - udpport to use for the communication with BRT
%
% Outputs:
% cueSpat       - the spatialized sound in .mat format
% cueParams     - parameters about the timing and spatial position of the
%               stimulus
%
% Required:
% - oscsend.m
% - BRT Renderer Application 
% - PsychToolbox for timing precision during recording
%
% #Author: Petra Kovacs

%% Initialize parameters output
cueParams = struct; 

%% Load and record cue
oscsend(u,'/removeAllSources', 's', '');
oscsend(u, '/source/loadSource', 'sss', 'Cue', cue, 'DirectivityModel');
WaitSecs(3);

% Set trajectory
[x,y,z] = sph2cart(deg2rad(azi),deg2rad(ele),r);
steps = round(totalDur/updateRate);

% Set start location
start = 1;
oscsend(u, '/source/location', 'sfff', 'Cue', x(start), y(start), z(start));
WaitSecs(3); 

% Start recording, then rapidly change location in a for loop
oscsend(u, '/source/playAndRecord', 'sssf', 'Cue', savename, 'mat', totalDur);

startTime = GetSecs();
for i = start+1 : steps
    oscsend(u, '/source/location', 'sfff', 'Cue', x(i), y(i), z(i));
    targetTime = startTime + i * updateRate; 
    WaitSecs('UntilTime', targetTime);
end
WaitSecs(3);
oscsend(u, '/stop', 's', '');
WaitSecs(3);

%% Create output
cueParams.sofa = load(savename);
cueSpat = cueParams.sofa.Data.Receiver;
cueParams.rGoal = r;
cueParams.aziGoal = azi;
cueParams.timeGoal = linspace(0,totalDur,length(r));
[aziActual,~,rActual] = cart2sph(cueParams.sofa.EmitterPosition(1,1,:), cueParams.sofa.EmitterPosition(1,2,:),cueParams.sofa.EmitterPosition(1,3,:));
cueParams.rActual = squeeze(rActual);
cueParams.aziActual = squeeze(aziActual);
cueParams.timeActual = cueParams.sofa.M;


