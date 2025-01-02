function generateStimuli(nStimuli,trajectory)
% Generate sound stimuli for wp1 pilot in three possible trajectories: looming,
% receding, rotating near, rotating far
%
% Inputs:
% nStimuli      - integer, number of stimuli to generate in the given trajectory
%                 (default = 120)
% trajectory    - 1: looming, 2: receding, 3: rotating near, 4: rotating
%                 far
%
% Output:
% Directory with the generated .wav files
%
% Author: Petra Kovacs

%% Set default parameters
ele = [0 0]; % elevation is always 0

if ~exist('nStimuli', 'var')
    nStimuli = 120;
elseif mod(nStimuli,2) ~= 0
    error('nStimuli has to be an even number.')
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

% Load HRTF dataset for spatialization
if not(exist("SOFAdbPath.m","file"))
    sofaPath = '\\kfs\fileserver\ProjektDaten\CherISH\code\SOFAtoolbox\SOFAtoolbox';
    addpath(sofaPath);
    SOFAstart;
end
database = 'scut';
HRTFfilename = 'SCUT_KEMAR_radius_all.sofa';
fullfn = fullfile(SOFAdbPath, 'database', database, HRTFfilename);
Obj = SOFAload(fullfn);
fs = Obj.Data.SamplingRate;

% Set radii for all stimuli
% Representative distances of each space in meter:
PPS = 0.2;
EPS = 2.0;
v = 1; % velocity in m/s

% Set values which have to be counterbalanced
offsetAzimuthOptions = [...
    repmat(30,1,nStimuli/2), ...
    repmat(-30,1,nStimuli/2)]; % left or right
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
    'offsetAzimuth', ...     % Side of the cue (at offset): 30 - left, -30 - right
    'target', ...            % 1 - target trial, 0 - nontarget trial
    'targetAzimuth', ...     % Side of the target: 30 - left, -30 - right               
    'congruence', ...        % 1 - congruent target, 0 - incongruent target
    };

outCsv = cell(nStimuli+1,length(outFields));
outCsv(1,:) = outFields;

%% Stimulus generation loop
for stimNo = 1:nStimuli
    frequency = (randi(4,1)+4)*100; % 500-800 Hz with round 100 values

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
    durMoving = abs(rMoving(2)-rMoving(1))/v;

    % Calculate total duration
    totalDur = durStatOnset+durMoving+durStatOffset; % in s

    % Set trajectory vector based on the durations
    rMain = [linspace(rOnset(1), rOnset(2), durStatOnset*fs), ...
            linspace(rMoving(1), rMoving(2), durMoving*fs), ...
            linspace(rOffset(1), rOffset(2), durStatOffset*fs)];

    %% Set azimuth
    offsetAzimuth = offsetAzimuthOptions(stimNo);

    if trajectory <= 2 % looming or receding
        % aziMain = [offsetAzimuth offsetAzimuth];
        aziMain = [linspace(offsetAzimuth, offsetAzimuth, durStatOnset*fs), ...
            linspace(offsetAzimuth, offsetAzimuth, durMoving*fs), ...
            linspace(offsetAzimuth, offsetAzimuth, durStatOffset*fs)];
    elseif trajectory >= 3 % rotating near or far
        % aziMain = [0 offsetAzimuth];
        aziMain = [linspace(0, offsetAzimuth, durStatOnset*fs), ...
            linspace(0, offsetAzimuth, durMoving*fs), ...
            linspace(0, offsetAzimuth, durStatOffset*fs)];        
    end

    %% Generate square waves and intensity ramps for each portion
    % Generate square waves for stationary portions
    tStatOnset = linspace(0,durStatOnset,durStatOnset*fs);
    statOnset = square(2*pi*frequency*tStatOnset);
    tStatOffset = linspace(0,durStatOffset,durStatOffset*fs);
    statOffset  = square(2*pi*frequency*tStatOffset);

    % Generate square waves for moving portion
    tMoving = linspace(0,durMoving,durMoving*fs);
    moving = square(2*pi*frequency*tMoving);

    %% Generate targets
    targetITI = 0.1; % 100 ms
    target = sin(2*pi*1200*(0:(1/fs):targetITI-1/fs)); % 100 ms, 1200 Hz
    gap = zeros(targetITI*fs,1);

    % Intensity: On-Off-ramp
    rampT_samples = round(fs*.01);                                   % Ramp in order to avoid startle respoce (> 12 ms)
    rampT = sin(linspace(0, pi/2, rampT_samples));                   % cosine ramp
    rampT_end = cos(linspace(0, pi/2, rampT_samples));               % cosine ramp
    target(1:rampT_samples) = target(1:rampT_samples)'.*rampT';
    target((end-rampT_samples+1):end) = target((end-rampT_samples+1):end)'.*rampT_end';

    % Spatialize target exactly between the two endpoints of the cue
    rTarget = mean([EPS,PPS]);
    rMain = [linspace(rOnset(1),rOnset(2),durStatOnset*fs), ...
        linspace(rMoving(1),rMoving(2),durMoving*fs), ...
        linspace(rOffset(1),rOffset(2),durStatOffset*fs), ...
        linspace(rTarget,rTarget,(targetITI*fs)*2)];

    % Target azimuth is congruent (same as offset azimuth) half the time
    % and incongruent (opposite of offset azimuth) half the time
    if targetTrials(stimNo) == 1 % if it is a target trial
        if congruence(stimNo) == 1
            targetAzimuth = offsetAzimuth; % 30 or -30 deg
        else
            targetAzimuth = -1*offsetAzimuth;
        end
    else
        targetAzimuth = nan;
    end

    aziMain = [aziMain, ...
        linspace(targetAzimuth,targetAzimuth,(targetITI*fs)*2)];

    %% Concatenate the sound portions and spatialize the stimuli
    if targetTrials(stimNo) == 1
        allPortions = [statOnset moving statOffset gap' target];
    else
        allPortions = [statOnset moving statOffset gap' gap'];
    end

    [stim, ~, ~, ~, ~] = local_SOFAspat(allPortions',Obj,aziMain,ele,rMain,EPS,PPS);

    % Save results to wav, add parameters to the cell array later saved out
    % to csv
    digits = ceil(log10(nStimuli + 1));
    stimulusdigits = ceil(log10(stimNo + 1));
    temp = char('');
    for digind = 1:(digits-stimulusdigits)
        temp = strcat(temp, '0');
    end

    filename = strcat(wavDir, '-', temp, num2str(stimNo));
    outCsv(stimNo+1, :) = {filename, frequency, totalDur, durStatOnset, durStatOffset,  rMoving(1), ...
        rMoving(2), trajectory, offsetAzimuth, targetTrials(stimNo), targetAzimuth, congruence(stimNo)};
    audiowrite(strcat('./', wavDir, '/', filename, '.wav'), stim, fs);

end % stimulus generation loop

%% Save output
% Convert cell to a table and use first row as variable names
T = cell2table(outCsv(2:end, :), 'VariableNames', outCsv(1, :));

% Write the table to a CSV file, final user message
writetable(T,strcat('./', wavDir, '/', strcat(wavDir, '-', 'StimuliData.csv')));
disp([newline, 'Task done, files and parameters are saved to directory ', wavDir, newline]);

end % function

function [out, aziActual, eleActual, rActual, idx] = local_SOFAspat(signal,Obj,azi,ele,r,EPS,PPS)
if length(r)~=length(signal)
    errorText = ['Signal (length: ', num2str(length(signal)), ') and r (length: ',...
        num2str(length(r)),') need to have the same length!'];
    error(errorText)
end
signal = db2mag(db(EPS,PPS)+65)*signal; % compensate for distance change from EPS to PPS and additional 65 dB to compensate for the particular type of HRTFs
[out, aziActual, eleActual, rActual, idx] = SOFAspat(signal./r(:),Obj,azi,ele,r);
end
