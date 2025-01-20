function BRTspat(source, totalDur, azi, ele, r, savename, updateRate, u)

% Spatialize stimuli by sending OSC commands to the BRT App.
% Required:
% - oscsend.m
% - BRT Renderer Application 
%
% Inputs:
% source        - path for a .wav file to be spatialized
% totalDur      - duration of the sound stimulus in s
% azi           - azimuth trajectory vector
% ele           - elevation trajectory vector
% r             - distance trajectory vector
% savename      - name of the .mat file to save the results into
% updateRate    - at what should the location be updated (ms)
% u             - udpport to use for the communication with BRT

%% Load source
oscsend(u, '/removeAllSources', 's', '');
oscsend(u, '/source/loadSource', 'sss', 'SoundSource1', source, 'DirectivityModel'); 

%% Send out dummy location updates for later timing precision
% Set trajectory
testSteps = 12000; 
testAzi = deg2rad(linspace(90,90,testSteps)); 
testEle = deg2rad(linspace(0,0,testSteps)); 
testDist = linspace(2,0.2,testSteps); 
[xT,yT,zT] = sph2cart(testAzi,testEle,testDist);

for update = 1:testSteps
    oscsend(u, '/source/location', 'sfff', 'SoundSource1', xT(update), yT(update), zT(update));
end

%% Set real trajectory
[x,y,z] = sph2cart(azi,ele,r);
steps = round(totalDur/updateRate);

for j = 1:4 % The first iterations don't catch the beginning distance precisely enough, but the rest pretty much does
    % Set start location
    start = 1;
    oscsend(u, '/source/location', 'sfff', 'SoundSource1', x(start), y(start), z(start));
    pause(1); % 1 s delay

    % Start recording, then rapidly change location in a for loop
    oscsend(u, '/source/playAndRecord', 'sssf', 'SoundSource1', savename, 'mat', totalDur);

    startTime = tic;
    for i = start+1 : steps
        oscsend(u, '/source/location', 'sfff', 'SoundSource1', x(i), y(i), z(i));
        elapsedTime = toc(startTime);
        targetTime = i * updateRate;
        if elapsedTime < targetTime
            pause(targetTime - elapsedTime);
        end
    end
    pause(1); % maybe then it crashes less often
end
oscsend(u, '/stop', 's', '');

