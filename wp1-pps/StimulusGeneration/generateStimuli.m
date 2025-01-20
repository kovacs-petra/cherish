function generateStimuli(nStimuli,trajectory)
% Generate sound stimuli for wp1 pilot in three possible trajectories: looming,
% receding, rotating near, rotating far
%
% Inputs:
% nStimuli      - integer, number of stimuli to generate in the given trajectory
%                 (usually 120; has to be divisible by 8 for counterbalancing)
% trajectory    - 1: looming, 2: receding, 3: rotating near, 4: rotating
%                 far
%
% Output:
% Directory with the generated .wav files, named after the date and time
% created
%
% Author: Petra Kovacs

%% Initial user message and input check
disp(['The BRT Renderer App should be open for this. If it isn''t, ', newline,...
    'it''s not too late to press Ctrl + C and rethink your life.',newline,newline]);

if mod(nStimuli,8) ~= 0
    error(['nStimuli has to be divisible by 8 to allow to counterbalance some parameters.', ...
        ' Try 8 for debugging or 120 for a typical batch of stimuli.'])
elseif trajectory > 4
    error('trajectory has to be an integer between 1 and 4');
end

% Create directory for saving audio data + parameters files
c = clock;  %#ok<CLOCK> % dir name based on current time
wavDir = strcat(date, '-', num2str(c(4)), num2str(c(5))); %#ok<DATE>
dircount = 0;
while exist(wavDir, 'dir')
    dircount = dircount + 1;
    if dircount > 1
        wavDir = strsplit(wavDir, '_');
        wavDir = wavDir{1};
    end
    wavDir = strcat(wavDir, '_', num2str(dircount));
end
mkdir(wavDir);
fs = 48000;

% Set radii for all stimuli
% Representative distances of each space in meter:
PPS = 0.2;
EPS = 2.0;
v = 1; % velocity in m/s

% Set values which have to be counterbalanced
offsetAzimuthOptions = [...
    repmat(90,1,nStimuli/2), ...
    repmat(-90,1,nStimuli/2)]; % left or right
targetTrials = [...
    zeros(1,nStimuli/4), ...
    ones(1,nStimuli/4), ...
    zeros(1,nStimuli/4), ...
    ones(1,nStimuli/4)]; % target or no target
congruence = [...
    nan(1,nStimuli/4), ... % nan where targetTrials is 0
    zeros(1,nStimuli/8), ...
    ones(1,nStimuli/8), ...
    nan(1,nStimuli/4), ...
    zeros(1,nStimuli/8), ...
    ones(1,nStimuli/8)];

%% Create a home for saved parameters
outFields = {
    'filename', ...
    'frequency', ...         % Cue frequency in Hz
    'totalDur', ...          % Cue duration in s
    'durStatOnset', ...      % Duration of stationary onset
    'durStatOffset', ...     % Duration of stationary offset
    'onsetDistance', ...     % Cue onset distance in m
    'offsetDistance', ...    % Cue offset distance in m
    'trajectory', ...        % 1 - loom, 2 - rec, 3 - rotate near, 4 - rotate far
    'offsetAzimuth', ...     % Side of the cue (at offset): 90 - left, -90 - right
    'target', ...            % 1 - target trial, 0 - nontarget trial
    'targetAzimuth', ...     % Side of the target: 90 - left, -90 - right               
    'congruence', ...        % 1 - congruent target, 0 - incongruent target
    };

outCsv = cell(nStimuli+1,length(outFields));
outCsv(1,:) = outFields;

%% Set connection between Matlab and BRT through OSC; set some settings
u = pnet('udpsocket',10017); % Listen port in BRT
oscsend(u, '/control/connect', 'si', 'localhost',10019);
oscsend(u, '/listener/enableSpatialization', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableInterpolation', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableNearFieldEffect', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableITD', 'sB', 'DefaultListener',1);
oscsend(u, '/environment/enableModel', 'sB', 'FreeField',1);
oscsend(u, '/environment/enableDirectPath', 'sB', 'FreeField',1);
oscsend(u, '/environment/enableReverbPath', 'sB', 'FreeField',1);
oscsend(u, '/listener/enableParallaxCorrection', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableModel', 'sB', 'DirectPath', 1);
oscsend(u, '/listener/enableModel', 'sB', 'ReverbPath', 1);
oscsend(u, '/environment/enablePropagationDelay', 'sB', 'FreeField', 1);
oscsend(u, '/environment/enableDistanceAttenuation', 'sB', 'FreeField', 1);

%% Stimulus generation loop
for stimNo = 1:nStimuli
    frequency = (randi(4,1)+5)*100; % 600-900 Hz with round 100 values

    % Set radii
    switch trajectory
        case 1 % Looming
            rOnset = [EPS EPS];  % stationary portion at the onset
            rOffset = [PPS PPS]; % stationary portion at the offset
            rMoving = [EPS PPS]; % moving portion

        case 2 % Receding
            rOnset = [PPS PPS];  % stationary portion at the onset
            rOffset = [EPS EPS]; % stationary portion at the offset
            rMoving = [PPS EPS]; % moving portion

        case 3 % Rotating near
            rOnset = [PPS PPS];  % stationary portion at the onset
            rOffset = rOnset;    % stationary portion at the offset
            rMoving = rOnset;    % moving portion

        case 4 % Rotating far
            rOnset = [EPS EPS];  % stationary portion at the onset
            rOffset = rOnset;    % stationary portion at the offset
            rMoving = rOnset;    % moving portion
    end

    %% Set durations
    % Stationary portions: 600-900 ms with round 100 values, in s
    durStatOnset = (randi(4,1)+5)/10; 
    durStatOffset = (randi(4,1)+5)/10;

    % Moving portion: duration calculated from distance and velocity
    durMoving = abs(EPS-PPS)/v;

    % Calculate total duration
    totalDur = durStatOnset+durMoving+durStatOffset; % in s

    updateRate = 2/1000; % every 2 ms
    ele = linspace(0,0,round(totalDur/updateRate));

    % Set trajectory vector based on the durations
    rMain = [linspace(rOnset(1), rOnset(2), round(durStatOnset/updateRate)), ...
            linspace(rMoving(1), rMoving(2), round(durMoving/updateRate)), ...
            linspace(rOffset(1), rOffset(2), round(durStatOffset/updateRate))];

    %% Set azimuth
    offsetAzimuth = offsetAzimuthOptions(stimNo);

    if trajectory <= 2 % looming or receding
        % aziMain = [offsetAzimuth offsetAzimuth];
        aziMain = [linspace(offsetAzimuth, offsetAzimuth, round(durStatOnset/updateRate)), ...
            linspace(offsetAzimuth, offsetAzimuth, round(durMoving/updateRate)), ...
            linspace(offsetAzimuth, offsetAzimuth, round(durStatOffset/updateRate))];
    elseif trajectory >= 3 % rotating near or far
        % aziMain = [0 offsetAzimuth];
        aziMain = [linspace(0, 0, round(durStatOnset/updateRate)), ...
            linspace(0, offsetAzimuth, round(durMoving/updateRate)), ...
            linspace(offsetAzimuth, offsetAzimuth, round(durStatOffset/updateRate))];        
    end

    %% Generate square waves for each portion
    % Generate square waves for stationary portions
    tStatOnset = linspace(0,durStatOnset,durStatOnset*fs);
    statOnset = square(2*pi*frequency*tStatOnset);
    tStatOffset = linspace(0,durStatOffset,durStatOffset*fs);
    statOffset  = square(2*pi*frequency*tStatOffset);

    % Ramp the onset and the offset
    rampSamples = round(fs*0.01);
    rampOnset = sin(linspace(0,pi/2,rampSamples));
    rampOffset = cos(linspace(0,pi/2,rampSamples));
    statOnset(1:rampSamples) = statOnset(1:rampSamples)'.*rampOnset';
    statOffset((end-rampSamples+1):end) = statOffset((end-rampSamples+1):end)'.*rampOffset';

    % Generate square waves for moving portion
    tMoving = linspace(0,durMoving,durMoving*fs);
    moving = square(2*pi*frequency*tMoving);

    %% Generate targets
    targetITI = 0.1; % 100 ms
    target = sin(2*pi*2400*(0:(1/fs):targetITI-1/fs)); % 100 ms, 2400 Hz
    gap = zeros(targetITI*fs,1);

    % Ramp the target as well
    target(1:rampSamples) = target(1:rampSamples)'.*rampOnset';
    target((end-rampSamples+1):end) = target((end-rampSamples+1):end)'.*rampOffset';

    % Spatialize target at the same distance where the cue ends
    rTarget = rOffset;
    rMain = [linspace(rOnset(1),rOnset(2),round(durStatOnset/updateRate)), ...
        linspace(rMoving(1),rMoving(2),round(durMoving/updateRate)), ...
        linspace(rOffset(1),rOffset(2),round(durStatOffset/updateRate)), ...
        linspace(rTarget(1),rTarget(2),round((targetITI/updateRate))*2)];

    % Target azimuth is congruent (same as offset azimuth) half the time
    % and incongruent (opposite of offset azimuth) half the time
    if targetTrials(stimNo) == 1 % if it is a target trial
        if congruence(stimNo) == 1
            targetAzimuth = offsetAzimuth; 
        else
            targetAzimuth = -1*offsetAzimuth;
        end
    else
        targetAzimuth = nan;
    end

    aziMain = [aziMain, ...
        linspace(targetAzimuth,targetAzimuth,round((targetITI/updateRate))*2)];

    ele = [ele, linspace(0,0,round((targetITI/updateRate))*2)];

    %% Concatenate the sound portions and spatialize the stimuli
    if targetTrials(stimNo) == 1
        allPortions = [statOnset moving statOffset gap' target];
    else
        allPortions = [statOnset moving statOffset gap' gap'];
    end

    % Save as wav (BRT needs it)
    digits = ceil(log10(nStimuli + 1));
    stimulusdigits = ceil(log10(stimNo + 1));
    temp = char('');
    for digind = 1:(digits-stimulusdigits)
        temp = strcat(temp, '0');
    end
    filename = fullfile(pwd, strcat(wavDir, '/', wavDir, '-', temp, num2str(stimNo)));
    wavname = strcat(filename,'.wav');
    audiowrite(wavname, int32(allPortions), fs);
    savename = strcat(filename,'.mat');

    % Spatialize with BRT
    BRTspat(wavname,totalDur,aziMain,ele,rMain,savename,updateRate,u);

    % Add parameters to cell array for later saving
    outCsv(stimNo+1, :) = {filename, frequency, totalDur, durStatOnset, durStatOffset,  rMoving(1), ...
        rMoving(2), trajectory, offsetAzimuth, targetTrials(stimNo), targetAzimuth, congruence(stimNo)};

    pause(1); % maybe then it crashes less often
end % stimulus generation loop

%% Save output
% Convert cell to a table and use first row as variable names
T = cell2table(outCsv(2:end, :), 'VariableNames', outCsv(1, :));

% Write the table to a CSV file, final user message
writetable(T,strcat('./', wavDir, '/', strcat(wavDir, '-', 'StimuliData.csv')));
disp([newline, 'Task done, files and parameters are saved to directory ', wavDir, newline]);
% 
% % Delete unnecessary files
% delete(wavname); % only BRT needed the wav
% if ~strcmp(savename(end-5:end), '_4.mat')
%     delete(savename);
% end

end % function
