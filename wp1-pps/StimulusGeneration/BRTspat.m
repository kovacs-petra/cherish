function cueSpat = BRTspat(cue, totalDur, aziC, eleC, rC, savenameC, updateRate, u)

% Spatialize stimuli by sending OSC commands to the BRT App.
% Required:
% - oscsend.m
% - BRT Renderer Application 
%
% Inputs:
% cue, target   - paths for cue and target .wav files to be spatialized
% totalDur      - duration of the sound stimulus in s
% aziC, T       - azimuth trajectory vector for the Cue and the Target
% eleC, T       - elevation trajectory vector
% rC, rT        - distance trajectory vector for the Cue and the Target
% savenameC, T  - name of the .mat files to save the Cue and Target into
% updateRate    - at what rate should the location be updated (ms)
% u             - udpport to use for the communication with BRT

%% Load and record cue
oscsend(u, '/removeAllSources', 's', '');
oscsend(u, '/source/loadSource', 'sss', 'Cue', cue, 'DirectivityModel');

% Set trajectory
[xC,yC,zC] = sph2cart(deg2rad(aziC),deg2rad(eleC),rC);
steps = round(totalDur/updateRate);

for j = 1:3 % The first iterations don't catch the beginning distance precisely enough, but the rest pretty much does
    % Set start location
    start = 1;
    oscsend(u, '/source/location', 'sfff', 'Cue', xC(start), yC(start), zC(start));
    pause(1); % 1 s delay

    % Start recording, then rapidly change location in a for loop
    oscsend(u, '/source/playAndRecord', 'sssf', 'Cue', savenameC, 'mat', totalDur);

    startTime = tic;
    for i = start+1 : steps
        oscsend(u, '/source/location', 'sfff', 'Cue', xC(i), yC(i), zC(i));
        elapsedTime = toc(startTime);
        targetTime = i * updateRate;
        if elapsedTime < targetTime
            pause(targetTime - elapsedTime);
        end
    end
end
oscsend(u, '/stop', 's', '');

%% Create output
load(savenameC,"Data");
cueSpat = Data.Receiver;


