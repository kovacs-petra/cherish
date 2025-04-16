function generateStimuli(nStimuli,trajectory,offsetAzimuth)
% Generate sound stimuli for CherISH wp1 in four possible trajectories: looming,
% receding, rotating near, rotating far
%
% Inputs:
% nStimuli      - integer, number of stimuli to generate in the given trajectory
%                 (has to be an even number for counterbalancing)
% trajectory    - 1: looming, 2: receding, 3: rotating near, 4: rotating
%                 far
% offsetAzimuth - of the cue, integer, e.g 90 or -90
%                (input rather than randomized: HRTF set is loaded only 1x
%                per run)
%
% Output:
% Directory with the generated .mat files, named after the date and time
% created
%
% #Author: Petra Kovacs, 2025, Acoustics Research Institute, Vienna, Austria

%% Initial user message and input check
disp(['The BRT Renderer App should be open for this. If it isn''t, ', newline,...
    'it''s not too late to press Ctrl + C and rethink your life.',newline]);

if mod(nStimuli,2) ~= 0
    warning(['nStimuli has to be even to allow to counterbalance some parameters. ', ...
        'I''m setting nStimuli to ',num2str(nStimuli+1),'!']);
    nStimuli = nStimuli+1;
elseif trajectory > 4
    disp('trajectory has to be an integer between 1 and 4, please try again.');
    disp(newline, '1 - loom, 2 - recede, 3 - rotate in PPS, 4 - rotate in EPS');
    trajectory = input(newline,'Input new trajectory here (1-4): ');
end

%% Initialize AMT 
if ~exist("amt_start.m","file")
    amtPath = '\\KFS\Fileserver\ProjektDaten\CherISH\code\amtoolbox-full-1.4.0';
    addpath(amtPath);
    amt_start;
end

%% Load an HRTF dataset to determine the right sampling rate
if ~exist("SOFAdbPath.m","file")
    sofaPath = '\\kfs\fileserver\ProjektDaten\CherISH\code\SOFAtoolbox\SOFAtoolbox';
    addpath(sofaPath);
    SOFAstart;
end
database = 'sadie';
HRTFfilename = '3DTI_HRTF_SADIE_II_D2_256s_48000Hz_resampled5.sofa';
% HRTFfilename = '3DTI_HRTF_SADIE_II_D2_256s_48000Hz.sofa';
% database = 'scut';
% HRTFfilename = 'HRTF_EPS.sofa';
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
EPS = 1.0;
% v = 0.5; % velocity in m/s

% Save data about direction (radial or angular)
if trajectory < 3
    direction = 1; % radial
else
    direction = 2; % angular
end

%% Preset some values
% Counterbalance target and no target + congruence-incongruence
targetTrial = repmat([zeros(1),ones(1)],1,nStimuli/2); % 0101
congruence = zeros(1,nStimuli);
congruence(4:4:end) = 1;  % 0001

% Randomize the order of targetTrial and congruence (together)
counterb = [targetTrial;congruence];
counterb = counterb(:,randperm(length(counterb)));
targetTrial = counterb(1,:);
congruence = counterb(2,:);

% Set placeholder for source intensity (will be set upon stimulus
% presentation)
sourceInt = zeros(1,nStimuli);

% Set f0s and durations for stationary onsets to allow unique combos
fmin = 310;
fmax = 2*fmin; % 1 octave
fShort = round(linspace(fmin,fmax,10),-1); % 10 unique f0 values
fOptions = repmat(fShort,1,5); % repeated 5 times for the 5 unique duration values
tOptions = zeros(1,length(fOptions)); % initialize vector for durations
n = 0;
for t = 0.4:0.1:0.8 % in sec
    tOptions(1,(1+n):(length(fShort)+n)) = t; % repeat each duration value for each f0 value
    n = n+length(fShort);
end

%% Create a home for saved parameters
outFields = {
    'filename', ...
    'frequency', ...         % Cue frequency in Hz
    'totalDur', ...          % Cue duration in s
    'durStatOnset', ...      % Duration of stationary onset in the cue
    ...'durStatOffset', ...     % Duration of stationary offset in the cue
    'onsetDistance', ...     % Cue onset distance in m
    'offsetDistance', ...    % Cue offset distance in m
    'direction',...          % 1 - radial, 2 - angular
    'trajectory', ...        % 1 - loom, 2 - rec, 3 - rotate near, 4 - rotate far
    'offsetAzimuth', ...     % Side of the cue (at offset): 90 - left, -90 - right
    'targetTrial', ...       % 1 - target trial, 0 - nontarget trial
    'congruence', ...        % 1 - congruent target, 0 - incongruent target
    'targetAzimuth', ...     % 90 or -90 (deg)
    'sourceInt', ...         % 1 - high source intensity, 0 - low source intensity (placeholder)
    'stimID', ...            % stimulus ID (which unique stimulus)
    'fs', ...                % sampling rate
    };

outCsv = cell(nStimuli+1,length(outFields));
outCsv(1,:) = outFields;

%% Set connection between Matlab and BRT through OSC; set some settings
% Priority(1); % high priority
u = pnet('udpsocket',10017); % Listen port in BRT
oscsend(u, '/control/connect', 'si', 'localhost',10019);
oscsend(u, '/listener/enableSpatialization', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableInterpolation', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableNearFieldEffect', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableITD', 'sB', 'DefaultListener',1);
oscsend(u, '/enableModel', 'sB', 'SDN',1);
% oscsend(u, '/environment/enableDirectPath', 'sB', 'SDN',1);
% oscsend(u, '/environment/enableReverbPath', 'sB', 'SDN',1);
oscsend(u, '/listener/enableParallaxCorrection', 'sB', 'DefaultListener',1);
oscsend(u, '/enableModel', 'sB', 'DirectPath', 1);
oscsend(u, '/enableModel', 'sB', 'ReverbPath', 1);
% oscsend(u, '/environment/enableDistanceAttenuation', 'sB', 'FreeField', 1);
oscsend(u, '/environment/setShoeBoxRoom', 'sfff', 'SDN', 2.8, 4.1, 3); % dimensions of the brown lab

%% Stimulus generation loop
for stimNo = 1:nStimuli

    % Progress update
    clc;
    disp(['Processing stimulus ',num2str(stimNo),'/',num2str(nStimuli),'...']);
    
    % Randomize f0 for cleaner ERPs (picked from previously set unique
    % values)
    frequency = fOptions(stimNo);

    % Set stimulus ID: trajectory, side, and f0 info
    switch offsetAzimuth; case 90; azID = 'L'; case -90; azID = 'R'; end
    stimID = [num2str(trajectory),azID,frequency];

    % Set radii (distance trajectories)
    switch trajectory
        case 1 % Looming
            rOnset = [EPS EPS];  % stationary portion at the onset
            % rOffset = [PPS PPS]; % stationary portion at the offset
            rMoving = [EPS PPS]; % moving portion

        case 2 % Receding
            rOnset = [PPS PPS];  % stationary portion at the onset
            % rOffset = [EPS EPS]; % stationary portion at the offset
            rMoving = [PPS EPS]; % moving portion

        case 3 % Rotating near
            rOnset = [PPS PPS];  % stationary portion at the onset
            % rOffset = rOnset;    % stationary portion at the offset
            rMoving = rOnset;    % moving portion

        case 4 % Rotating far
            rOnset = [EPS EPS];  % stationary portion at the onset
            % rOffset = rOnset;    % stationary portion at the offset
            rMoving = rOnset;    % moving portion
    end

    %% Set durations
    % Stationary portions
    buffer = 1; % 1 s buffer for BRT artefacts
    durStatOnset = tOptions(stimNo)+buffer; 
    % durStatOffset = ((randi(4,1)+4)/10);

    % Moving portion: duration calculated from distance and velocity
    % durMoving = abs(EPS-PPS)/v;
    durMoving = 2;

    % Calculate total duration 
    % totalDur = durStatOnset+durMoving+durStatOffset; % in s
    totalDur = durStatOnset+durMoving;

    % Set rate for updating the sound position in BRT
    updateRate = 2/1000; % every 2 ms

    % Elevation always 0
    % ele = linspace(0,0,round(totalDur*fs));
    eleC = linspace(0,0,round(totalDur/updateRate));

    %% Set distance trajectory vector based on the durations
    rC = [linspace(rOnset(1), rOnset(2), round(durStatOnset/updateRate)), ...
            linspace(rMoving(1), rMoving(2), round(durMoving/updateRate))];

    %% Set azimuth trajectory
    if trajectory <= 2 % looming or receding: same azimuth throughout
        aziC = [linspace(offsetAzimuth, offsetAzimuth, round(durStatOnset/updateRate)), ...
            linspace(offsetAzimuth, offsetAzimuth, round(durMoving/updateRate))];
    elseif trajectory >= 3 % rotating near or far: the laterality changes
        aziC = [linspace(-offsetAzimuth, -offsetAzimuth, round(durStatOnset/updateRate)), ...
            linspace(-offsetAzimuth, offsetAzimuth, round(durMoving/updateRate))];
    end

    %% Generate triangle waves for each portion
    % For stationary portions
    tStatOnset = linspace(0,durStatOnset,durStatOnset*fs); % Time vector from 0 to dur with sampling interval 1/fs
    statOnset = sawtooth(2 * pi * frequency * tStatOnset, 0.5); % 0.5 for a symmetric triangular wave
    % tStatOffset = linspace(0,durStatOffset,durStatOffset*fs);
    % statOffset  = sawtooth(2 * pi * frequency * tStatOffset, 0.5);

    % For moving portion
    tMoving = linspace(0,durMoving,durMoving*fs);
    moving = sawtooth(2 * pi * frequency * tMoving, 0.5);

    %% Generate targets
    targetITI = 0.1; % 100 ms
    target = sin(2*pi*2400*(0:(1/fs):targetITI-1/fs)); % 100 ms, 2400 Hz
    gap = zeros(targetITI*fs,1); % silence

    % Ramp the target
    rampSamples = round(fs*0.01);
    rampOnset = sin(linspace(0,pi/2,rampSamples));
    rampOffset = cos(linspace(0,pi/2,rampSamples));
    target(1:rampSamples) = target(1:rampSamples)'.*rampOnset';
    target((end-rampSamples+1):end) = target((end-rampSamples+1):end)'.*rampOffset';

    %% Concatenate and spatialize the cue
    % C = [statOnset moving statOffset];
    C = [statOnset moving];
    win = tukeywin(size(C,2),(0.1*fs)/size(C,2)); % 0.1 s on-offset ramp
    C = C.*win';

    % Save as wav (BRT needs it)
    digits = ceil(log10(nStimuli + 1));
    stimulusdigits = ceil(log10(stimNo + 1));
    temp = char('');
    for digind = 1:(digits-stimulusdigits)
        temp = strcat(temp, '0');
    end
    filename = fullfile(pwd, strcat(wavDir, '/', wavDir, '-', temp, num2str(stimNo)));
    wavname = strcat(filename,'.wav');   
    audiowrite(wavname, C, fs); 
    savenameC = strcat(filename,'.mat'); 

    % Spatialize cue with BRT
    % WaitSecs(1);
    [cueSpat, cueParams] = BRTspat(wavname,totalDur,aziC,eleC,rC,savenameC,updateRate,u);

    % Cut off the onset and offset artefacts introduced by BRT during recording; 
    % apply Tukey window on the cue

    % Find the first value in the cue which is ~= to start position
    switch trajectory
        case {1,2} % loom, rec
            idx = find(round(rMoving(1),1)==round(cueParams.rActual,1));
        case {3,4} % rotate in PPS, EPS
            idx = find(round(aziC(1),1)==round(rad2deg(cueParams.aziActual),1));
    end
    % idx(end) is the last indx in rActual which is = starting position,
    % i.e. idx(end)+1 is where the sound starts to move
    motionStart = cueParams.timeActual(length(idx)+1); % time in s

    % Count back durStatOnset-buffer ms from this point and cut all before
    % it, then apply window
    cueSpat = cueSpat(:,round((motionStart-durStatOnset+buffer)*fs):end);
    win = tukeywin(size(cueSpat,2),(0.1*fs)/size(cueSpat,2)); % 0.1 s 
    cueSpatWin = cueSpat.*repmat(win',2,1);

    % Update duration data to exclude buffer
    totalDur = length(cueSpatWin)/fs;
    durStatOnset = durStatOnset-buffer;

    % Scale cue
    if trajectory < 4 % sound is in PPS at some point
        scaledb = 20 * log10(1 / max(cueSpatWin,[],"all")); % Scale to 1
    else % sound is in EPS all the way
        scaleto = PPS/EPS;
        scaledb = 20 * log10(scaleto / max(cueSpatWin,[],"all")); % Scale to PPS/EPS (far sound is always softer)
    end
    scalecue1 = scaletodbspl(cueSpatWin(1,:)',dbspl(cueSpatWin(1,:)')+scaledb);
    scalecue2 = scaletodbspl(cueSpatWin(2,:)',dbspl(cueSpatWin(2,:)')+scaledb);
    scalecue = [scalecue1 scalecue2];

    %% Spatialize target separately, then concatenate with the cue
    % Because BRT is not good for very short sounds like the target
    if targetTrial(stimNo) == 1
        T = [gap' target; gap' target]; % target
    else
        T = [gap' gap'; gap' gap']; % no target
    end

    if ~isequal(T,zeros(size(T))) % if there is a target
        % rTarget = linspace(mean([EPS,PPS]),mean([EPS,PPS]),(targetITI*fs)*2); % Target exactly between EPS and PPS
        rTarget = linspace(PPS,PPS,(targetITI*fs)*2); % target in PPS

        % Target azimuth is congruent (same as offset azimuth) half the time
        % and incongruent (opposite of offset azimuth) half the time
        if congruence(stimNo) == 1
            targetAzimuth = offsetAzimuth;
            % ObjTarget = ObjCue; % same HRTF set
        else
            targetAzimuth = -offsetAzimuth;
        end
            
        % switch targetAzimuth
        %     case 90
        %         targetHRTF = 'HRTF_left.sofa';
        %     case -90
        %         targetHRTF = 'HRTF_right.sofa';
        % end

        % Load HRTF set for target
        % fullfnT = fullfile(SOFAdbPath, 'database', database, '3DTI_HRTF_SADIE_II_D2_256s_48000Hz_resampled5.sofa');
        ObjTarget = ObjCue;

        % % Sanity check
        % if isequal(ObjTarget.Data.IR, ObjCue.Data.IR)
        %     error('Wrong HRTF loaded for incongruent target.');
        % end
        % end

        aziTarget = linspace(targetAzimuth,targetAzimuth,(targetITI*fs)*2);
        eleTarget = linspace(0,0,(targetITI*fs)*2);
        targetStim = local_SOFAspat(T(1,:)',ObjTarget,aziTarget,eleTarget,rTarget);
        % stim = [cue' targetStim'];

        % Scale the target
        scaledb = 20 * log10(1 / max(targetStim,[],"all")); % Target is always in PPS, scale to 1
        scaletarget1 = scaletodbspl(targetStim(:,1)',dbspl(targetStim(:,1)')+scaledb); % target, channel 1
        scaletarget2 = scaletodbspl(targetStim(:,2)',dbspl(targetStim(:,2)')+scaledb); % target, channel 2
        scaletarget = [scaletarget1; scaletarget2]; % whole binaural target
        out = [scalecue' scaletarget]; % whole binaural cue and target together
    else
        targetAzimuth = nan;
        out = [scalecue' T]; % cue plus silence
    end % if there is a target

    % Save stimulus in .mat file
    digits = ceil(log10(nStimuli + 1));
    stimulusdigits = ceil(log10(stimNo + 1));
    temp = char('');
    for digind = 1:(digits-stimulusdigits)
        temp = strcat(temp, '0');
    end
    filename = strcat(wavDir, '/', wavDir, '-', temp, num2str(stimNo));
    save(strcat(filename,'.mat'), "out", "fs", "cueParams");

    % Add parameters to cell array for later saving
    outCsv(stimNo+1, :) = {filename, frequency, totalDur, durStatOnset, rMoving(1), ...
        rMoving(2), direction, trajectory, offsetAzimuth, targetTrial(stimNo),...
        congruence(stimNo), targetAzimuth, sourceInt(stimNo), stimID, fs};

    WaitSecs(3); % avoid BRT crashing
    

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