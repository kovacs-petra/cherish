function generateStimuli(nStimuli,trajectory,offsetAzimuth,targetSpat)
% Generate sound stimuli for wp1 pilot in four possible trajectories: looming,
% receding, rotating near, rotating far
%
% Inputs:
% nStimuli      - integer, number of stimuli to generate in the given trajectory
%                 (has to be divisible by 4 for counterbalancing)
% trajectory    - 1: looming, 2: receding, 3: rotating near, 4: rotating
%                 far
% offsetAzimuth - of the cue, integer, e.g 90 or -90
% targetSpat    - 1: spatialized target, 0: diotic target
%
% Output:
% Directory with the generated .mat files, named after the date and time
% created
%
% Author: Petra Kovacs, 2025

%% Initial user message and input check
if mod(nStimuli,4) ~= 0
    error(['nStimuli has to be divisible by 4 to allow to counterbalance some parameters.', ...
        ' Try 4 for debugging or 120 for a typical batch of stimuli.'])
elseif trajectory > 4
    error('trajectory has to be an integer between 1 and 4');
end

%% Initialize AMT 
if ~exist("amt_start.m","file")
    amtPath = '\\KFS\Fileserver\ProjektDaten\CherISH\code\amtoolbox-full-1.4.0';
    addpath(amtPath);
    amt_start;
end

%% Load HRTF datasets for spatialization
if not(exist("SOFAdbPath.m","file"))
    sofaPath = '\\kfs\fileserver\ProjektDaten\CherISH\code\SOFAtoolbox\SOFAtoolbox';
    addpath(sofaPath);
    SOFAstart;
end
database = 'scut';
switch trajectory
    case 1 % looming
        switch offsetAzimuth
            case 90
                HRTFfilename = 'HRTF_left.sofa';
            case -90
                HRTFfilename = 'HRTF_right.sofa';
        end
    case 2 % receding
        switch offsetAzimuth
            case 90
                HRTFfilename = 'HRTF_left.sofa';
            case -90
                HRTFfilename = 'HRTF_right.sofa';
        end
    case 3 % rotating PPS
        HRTFfilename = 'HRTF_PPS.sofa';
    case 4 % rotating EPS
        HRTFfilename = 'HRTF_EPS.sofa';
end

fullfn = fullfile(SOFAdbPath, 'database', database, HRTFfilename);
ObjCue = SOFAload(fullfn);
fs = ObjCue.Data.SamplingRate;

%% Create directory for saving audio data + parameters files
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

% Set radii for all stimuli
% Representative distances of each space in meter:
PPS = 0.2;
EPS = 2.0;
v = 1; % velocity in m/s

% Set values which have to be counterbalanced
switch targetSpat
    case 0 % Counterbalance target and no target
        targetTrial = [...
            zeros(1,nStimuli/4), ...
            ones(1,nStimuli/4), ...
            zeros(1,nStimuli/4), ...
            ones(1,nStimuli/4)]; 
        congruence = nan(1,nStimuli);
    case 1 % Counterbalance congruent and incongruent
        congruence = [...
            zeros(1,nStimuli/4), ...
            ones(1,nStimuli/4), ...
            zeros(1,nStimuli/4), ...
            ones(1,nStimuli/4)];
        targetTrial = ones(1,nStimuli);
end

%% Create a home for saved parameters
outFields = {
    'filename', ...
    'frequency', ...         % Cue frequency in Hz
    'totalDur', ...          % Cue duration in s
    'durStatOnset', ...      % Duration of stationary onset in the cue
    'durStatOffset', ...     % Duration of stationary offset in the cue
    'onsetDistance', ...     % Cue onset distance in m
    'offsetDistance', ...    % Cue offset distance in m
    'trajectory', ...        % 1 - loom, 2 - rec, 3 - rotate near, 4 - rotate far
    'offsetAzimuth', ...     % Side of the cue (at offset): 90 - left, -90 - right
    'targetTrial', ...       % 1 - target trial, 0 - nontarget trial
    'congruence',...         % 1 - congruent, 0 - incongruent
    };

outCsv = cell(nStimuli+1,length(outFields));
outCsv(1,:) = outFields;

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
    % Stationary portions: 600-800 ms with round 100 values, in s
    durStatOnset = ((randi(4,1)+4)/10); 
    durStatOffset = ((randi(4,1)+4)/10);

    % Moving portion: duration calculated from distance and velocity
    durMoving = abs(EPS-PPS)/v;

    % Calculate total duration 
    totalDur = durStatOnset+durMoving+durStatOffset; % in s

    % Elevation always 0
    ele = linspace(0,0,round(totalDur*fs));

    %% Set distance trajectory vector based on the durations
    r = [linspace(rOnset(1), rOnset(2), durStatOnset*fs), ...
            linspace(rMoving(1), rMoving(2), durMoving*fs), ...
            linspace(rOffset(1), rOffset(2), durStatOffset*fs)];

    %% Set azimuth trajectory
    if trajectory <= 2 % looming or receding
        azi = [linspace(offsetAzimuth, offsetAzimuth, durStatOnset*fs), ...
            linspace(offsetAzimuth, offsetAzimuth, durMoving*fs), ...
            linspace(offsetAzimuth, offsetAzimuth, durStatOffset*fs)];
    elseif trajectory >= 3 % rotating near or far
        azi = [linspace(-offsetAzimuth, -offsetAzimuth, durStatOnset*fs), ...
            linspace(-offsetAzimuth, offsetAzimuth, durMoving*fs), ...
            linspace(offsetAzimuth, offsetAzimuth, durStatOffset*fs)];        
    end

    %% Generate triangle waves for each portion
    % For stationary portions
    tStatOnset = linspace(0,durStatOnset,durStatOnset*fs); % Time vector from 0 to dur with sampling interval 1/fs
    statOnset = sawtooth(2 * pi * frequency * tStatOnset, 0.5); % 0.5 for a symmetric triangular wave
    tStatOffset = linspace(0,durStatOffset,durStatOffset*fs);
    statOffset  = sawtooth(2 * pi * frequency * tStatOffset, 0.5);

    % For moving portion
    tMoving = linspace(0,durMoving,durMoving*fs);
    moving = sawtooth(2 * pi * frequency * tMoving, 0.5);

    %% Generate targets
    targetITI = 0.1; % 100 ms
    target = sin(2*pi*2400*(0:(1/fs):targetITI-1/fs)); % 100 ms, 2400 Hz
    gap = zeros(targetITI*fs,1);

    % Ramp the target
    rampSamples = round(fs*0.01);
    rampOnset = sin(linspace(0,pi/2,rampSamples));
    rampOffset = cos(linspace(0,pi/2,rampSamples));
    target(1:rampSamples) = target(1:rampSamples)'.*rampOnset';
    target((end-rampSamples+1):end) = target((end-rampSamples+1):end)'.*rampOffset';

    %% Concatenate the sound portions and spatialize the stimuli
    C = [statOnset moving statOffset]; % cue
    win = tukeywin(size(C,2),(0.1*fs)/size(C,2)); % 0.1 s on-offset ramp
    C = C.*win';

    if targetSpat == 0
        if targetTrial(stimNo) == 1
            T = [gap' target; gap' target]; % target
        else
            T = [gap' gap'; gap' gap']; % no target
        end
    else
        T = [gap' target];
    end

    % Spatialize cue
    cue = local_SOFAspat(C',ObjCue,azi,ele,r);

    % Scale cue
    if trajectory < 4
        scaledb = 20 * log10(1 / max(cue,[],"all"));
    else
        scaledb = 20 * log10(0.1 / max(cue,[],"all"));
    end
    scalecue1 = scaletodbspl(cue(:,1)',dbspl(cue(:,1)')+scaledb);
    scalecue2 = scaletodbspl(cue(:,2)',dbspl(cue(:,2)')+scaledb);
    scalecue = [scalecue1; scalecue2];

    switch targetSpat
        case 0              
                out = [scalecue T];

        case 1 % Spatialize target separately, then concatenate
            % rTarget = linspace(mean([EPS,PPS]),mean([EPS,PPS]),(targetITI*fs)*2); % Target exactly between EPS and PPS
            rTarget = linspace(PPS,PPS,(targetITI*fs)*2); % target in PPS

             % Target azimuth is congruent (same as offset azimuth) half the time
             % and incongruent (opposite of offset azimuth) half the time
             if congruence(stimNo) == 1
                 targetAzimuth = offsetAzimuth; 
                 ObjTarget = ObjCue;
             else 
                 targetAzimuth = -offsetAzimuth;
                 switch targetAzimuth
                     case 90
                         targetHRTF = 'HRTF_left.sofa';
                     case -90
                         targetHRTF = 'HRTF_right.sofa';
                 end
                 fullfnT = fullfile(SOFAdbPath, 'database', database, targetHRTF);
                 ObjTarget = SOFAload(fullfnT);
             end

             aziTarget = linspace(targetAzimuth,targetAzimuth,(targetITI*fs)*2);
             targetStim = local_SOFAspat(T',ObjTarget,aziTarget,ele,rTarget);
             % stim = [cue' targetStim'];

             % Scale the target
             scaledb = 20 * log10(1 / max(targetStim,[],"all"));
             scaletarget1 = scaletodbspl(targetStim(:,1)',dbspl(targetStim(:,1)')+scaledb);
             scaletarget2 = scaletodbspl(targetStim(:,2)',dbspl(targetStim(:,2)')+scaledb);
             scaletarget = [scaletarget1; scaletarget2];
             out = [scalecue scaletarget];
             % out = [scalecue'; targetStim]';
    end

    % Save stimulus in .mat file
    digits = ceil(log10(nStimuli + 1));
    stimulusdigits = ceil(log10(stimNo + 1));
    temp = char('');
    for digind = 1:(digits-stimulusdigits)
        temp = strcat(temp, '0');
    end
    filename = strcat(wavDir, '/', wavDir, '-', temp, num2str(stimNo));
    save(strcat(filename,'.mat'), "out");

    % Add parameters to cell array for later saving
    outCsv(stimNo+1, :) = {filename, frequency, totalDur, durStatOnset, durStatOffset,  rMoving(1), ...
        rMoving(2), trajectory, offsetAzimuth, targetTrial(stimNo), congruence(stimNo)};

end % stimulus generation loop

%% Save output
% Convert cell to a table and use first row as variable names
T = cell2table(outCsv(2:end, :), 'VariableNames', outCsv(1, :));

% Write the table to a CSV file, final user message
writetable(T,strcat('./', wavDir, '/', strcat(wavDir, '-', 'StimuliData.csv')));

disp([newline, 'Task done, files and parameters are saved to directory ', wavDir, newline]);

end % function

% Version of SOFAspat that adds an intensity ramp after checking if it's
% possible
function [out, aziActual, eleActual, rActual, idx] = local_SOFAspat(signal,Obj,azi,ele,r)
if length(r)~=length(signal)
    errorText = ['Signal (length: ', num2str(length(signal)), ') and r (length: ',...
        num2str(length(r)),') need to have the same length!'];
    error(errorText)
end
[out, aziActual, eleActual, rActual, idx] = SOFAspat(signal./r(:),Obj,azi,ele,r);
end