function genStimExt(dur)
% Usage: genStimExt(dur)
%
% Generates 10x3x2 sounds of dur s duration, at 10 different f0s (in one
% octave range), 3 different distances (in head, 20 cm, 100 cm) and 2 sides
% (left, right). The binaural stimuli are generated with the Binaural
% Rendering Toolbox (BRT) in a separate loop than the diotic stimuli.
% 
% Inputs:
% - dur:   duration of each sound in s (default: 2 s)
%
% #Author: Petra Kovacs, 2025, Acoustics Research Institute, Vienna,
% Austria
%
%% Input check and defaults
if ~exist('dur','var')
  dur = 2;
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
fullfn = fullfile(SOFAdbPath, 'database', database, HRTFfilename);
Obj = SOFAload(fullfn);
fs = Obj.Data.SamplingRate;

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

%% Preset values for binaural stimuli
% No. of stimuli and conditons
nF0 = 10; % unique f0 values
nD = 2; % unique distance values for dichotic stimuli (BRT)
nSide = 2; % unique azimuth values
nBinauralStimuli = nF0*nD*nSide;

% Distances (m)
PPS = 0.2;
EPS = 1;
dOptions = repmat([PPS,PPS,EPS,EPS],1,nBinauralStimuli/4);

% Azimuths (deg)
aziOptions = repmat([90,-90],1,nBinauralStimuli/nSide);

% F0
fmin = 310;
fmax = 2*fmin; % 1 octave
fShort = round(linspace(fmin,fmax,nF0),-1);
fOptions = zeros(1,nBinauralStimuli);
j = 1;
for i = 1:4:length(fOptions)
    fOptions(i:i+3) = fShort(j);
    j = j+1;
end

% Source intensity (placeholder)
sourceInt = zeros(1,nBinauralStimuli);

%% Preset values for diotic stimuli
% No. of stimuli and conditons
nDdi = 1; % unique distance values for diotic stimuli
nDioticStimuli = nF0*nDdi*nSide;

% Distances (m)
dOptionsDi = zeros(1,nDioticStimuli);

% Dummy azimuths (deg)
aziOptionsDi = repmat([90,-90],1,nDioticStimuli/nSide);

% F0
fOptionsDi = zeros(1,nDioticStimuli);
j = 1;
for i = 1:2:length(fOptionsDi)
    fOptionsDi(i:i+1) = fShort(j);
    j = j+1;
end

% Source intensity (placeholder)
sourceIntDi = zeros(1,nDioticStimuli);

%% Create a home for saved parameters
outFields = {
    'filename', ...
    'f0', ...                % F0 in Hz
    'dur', ...               % Duration in s
    'distance', ...          % Distance in m
    'azimuth', ...           % Side: 90 - left, -90 - right
    'sourceInt', ...         % 1 - high source intensity, 0 - low source intensity (placeholder)
    'fs', ...                % sampling rate
    };

outCsv = cell(nBinauralStimuli+nDioticStimuli+1,length(outFields));
outCsv(1,:) = outFields;

%% Set connection between Matlab and BRT through OSC; set some settings
u = pnet('udpsocket',10017); % Listen port in BRT
oscsend(u, '/control/connect', 'si', 'localhost',10019);
oscsend(u, '/listener/enableSpatialization', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableInterpolation', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableNearFieldEffect', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableITD', 'sB', 'DefaultListener',1);
oscsend(u, '/enableModel', 'sB', 'FreeField',1);
oscsend(u, '/environment/enableDirectPath', 'sB', 'FreeField',1);
oscsend(u, '/listener/enableParallaxCorrection', 'sB', 'DefaultListener',1);
oscsend(u, '/enableModel', 'sB', 'DirectPath', 1);
oscsend(u, '/environment/enableDistanceAttenuation', 'sB', 'FreeField', 1);

%% Stimulus generation loop for dichotic sounds (BRT)
for stimNo = 1:nBinauralStimuli

    % Progress update
    clc;
    disp(['Processing stimulus ',num2str(stimNo),'/',num2str(nBinauralStimuli+nDioticStimuli),'...']);
    
    % Randomize f0 (picked from previously set unique values)
    f0 = fOptions(stimNo);

    %% Timing variables
    buffer = 1; % 1 s buffer for BRT artefacts

    %% Position variables
    azi = aziOptions(stimNo);
    ele = 0; % Elevation always 0
    r = dOptions(stimNo);

    %% Generate triangle wave
    t = linspace(0,dur,dur*fs); % Time vector from 0 to dur with sampling interval 1/fs
    triwave = sawtooth(2 * pi * f0 * t, 0.5); % 0.5 for a symmetric triangular wave

    %% Prepare for spatialization
    win = tukeywin(size(triwave,2),(0.1*fs)/size(triwave,2)); % 0.1 s on-offset ramp
    triwave = triwave.*win';

    % Save as wav (BRT needs it)
    digits = ceil(log10(nBinauralStimuli + 1));
    stimulusdigits = ceil(log10(stimNo + 1));
    temp = char('');
    for digind = 1:(digits-stimulusdigits)
        temp = strcat(temp, '0');
    end
    filename = fullfile(pwd, strcat(wavDir, '/', wavDir, '-', temp, num2str(stimNo)));
    wavname = strcat(filename,'.wav');   
    audiowrite(wavname, triwave, fs); 
    savename = strcat(filename,'.mat'); 

    %% Spatialize with BRT
    oscsend(u,'/removeAllSources', 's', '');
    oscsend(u, '/source/loadSource', 'sss', 'source', wavname, 'DirectivityModel'); 
    [x,y,z] = sph2cart(deg2rad(azi),deg2rad(ele),r);
    oscsend(u, '/source/location', 'sfff', 'source', x, y, z); WaitSecs(2);
    oscsend(u, '/source/playAndRecord', 'sssf', 'source', savename, 'mat', dur+buffer); WaitSecs(3);
    oscsend(u, '/stop', 's', ''); WaitSecs(1);

    %% Create output
    params.sofa = load(savename);
    spat = params.sofa.Data.Receiver;
    params.rGoal = r;
    params.aziGoal = azi;
    params.timeGoal = linspace(0,dur,length(r));
    [aziActual,~,rActual] = cart2sph(params.sofa.EmitterPosition(1,1,:), params.sofa.EmitterPosition(1,2,:),params.sofa.EmitterPosition(1,3,:));
    params.rActual = squeeze(rActual);
    params.aziActual = squeeze(aziActual);
    params.timeActual = params.sofa.M;

    % Cut off the onset and offset artefacts introduced by BRT during recording; 
    % apply Tukey window on the sound
    spat = spat(:,1:round(dur*fs));
    win = tukeywin(size(spat,2),(0.1*fs)/size(spat,2)); % 0.1 s 
    spatWin = spat.*repmat(win',2,1);

    % Update duration data to exclude buffer
    totalDur = length(spatWin)/fs;

    % Scale sound to insure intensity differences
    switch r
        case PPS % Scale to intern/PPS
            scaleto = 1;
        case EPS % Scale to intern/EPS
            scaleto = PPS/EPS;
    end
    
    scaledb = 20 * log10(scaleto / max(spatWin,[],"all"));

    scaled1 = scaletodbspl(spatWin(1,:)',dbspl(spatWin(1,:)')+scaledb);
    scaled2 = scaletodbspl(spatWin(2,:)',dbspl(spatWin(2,:)')+scaledb);
    out = [scaled1 scaled2];

    % Save stimulus in .mat file
    digits = ceil(log10(nBinauralStimuli + 1));
    stimulusdigits = ceil(log10(stimNo + 1));
    temp = char('');
    for digind = 1:(digits-stimulusdigits)
        temp = strcat(temp, '0');
    end
    filename = strcat(wavDir, '/', wavDir, '-', temp, num2str(stimNo));
    save(strcat(filename,'.mat'), "out", "fs", "params");

    % Add parameters to cell array for later saving
    outCsv(stimNo+1, :) = {filename, f0, totalDur, r, azi, sourceInt(stimNo), fs};

    WaitSecs(2); % avoid BRT crashing

end % stimulus generation loop for BRT

%% Stimulus generation for diotic sounds not spatialized
% These sounds exemplify an in-the-head percept. Still, a dummy azimuth
% value is specified in order not to confuse the stimulus presentation
% script.

for stimDi = 1:nBinauralStimuli/2
% Progress update
    clc;
    disp(['Processing stimulus ',num2str(nBinauralStimuli+stimDi),'/',num2str(nBinauralStimuli+nDioticStimuli),'...']);
    
    % Randomize f0 (picked from previously set unique values)
    f0 = fOptionsDi(stimDi);

    %% Position variables
    azi = aziOptionsDi(stimDi);
    r = dOptionsDi(stimDi);

    %% Generate triangle wave
    t = linspace(0,dur,dur*fs); % Time vector from 0 to dur with sampling interval 1/fs
    triwave = sawtooth(2 * pi * f0 * t, 0.5); % 0.5 for a symmetric triangular wave

    %% Window and stereo
    win = tukeywin(size(triwave,2),(0.1*fs)/size(triwave,2)); % 0.1 s on-offset ramp
    triwave = repmat(triwave.*win',2,1);

    % Scale sound to insure intensity differences
    scaleto = 1;    
    scaledb = 20 * log10(scaleto / max(triwave,[],"all"));

    scaled1 = scaletodbspl(triwave(1,:)',dbspl(triwave(1,:)')+scaledb);
    scaled2 = scaletodbspl(triwave(2,:)',dbspl(triwave(2,:)')+scaledb);
    out = [scaled1 scaled2];

    % Save stimulus in .mat file
    digits = ceil(log10(nDioticStimuli + 1));
    stimulusdigits = ceil(log10(stimNo + 1));
    temp = char('');
    for digind = 1:(digits-stimulusdigits)
        temp = strcat(temp, '0');
    end
    filename = strcat(wavDir, '/', wavDir, '-', temp, num2str(nBinauralStimuli+stimDi));
    save(strcat(filename,'.mat'), "out", "fs", "params");

    % Add parameters to cell array for later saving
    outCsv(nBinauralStimuli+stimDi+1, :) = {filename, f0, totalDur, r, azi, sourceIntDi(stimDi), fs};

end

%% Save output
% Convert cell to a table and use first row as variable names
T = cell2table(outCsv(2:end, :), 'VariableNames', outCsv(1, :));

% Write the table to a CSV file, final user message
writetable(T,strcat('./', wavDir, '/', strcat(wavDir, '-', 'StimuliData.csv')));

disp([newline, 'Task done, files and parameters are saved to directory ', wavDir, newline]);
end
