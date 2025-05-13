function stimStruct = genStimExt(dur,f0)
% Usage: genStimExt(dur)
%
% Generates 10x3x2 sounds of dur s duration, at 10 different f0s (in one
% octave range), 3 different distances (in head, 20 cm, 100 cm) and 2 sides
% (left, right). The binaural stimuli are generated with the Binaural
% Rendering Toolbox (BRT) in a separate loop than the diotic stimuli.
%
% Inputs:
% - dur:   duration of each sound in s (default: 2 s)
% - f0:    set constant f0 for all stimuli (default: varies in one octave range)
%
% #Author: Petra Kovacs, 2025, Acoustics Research Institute, Vienna,
% Austria
%
%% Preset values for binaural stimuli
% No. of stimuli and conditons
nF0 = 8; % unique f0 values
nD = 7; % unique distance values
nSide = 2; % unique azimuth values
nStimuli = nF0*nD*nSide;
stimStruct.stim = cell(nF0,nD,nSide);

% Distances (m)
in = 0;
PPS = 0.2;
EPS = 2;
dOptions = [in,PPS,EPS,round(linspace(PPS+0.2,EPS-0.2,nD-3),1)];  
                                              % 1st three columns: the fixed
                                              % distances we want to use in the
                                              % experiment.
                                              % Then, the distances inbetween.

% Azimuths (deg)
% aziOptions = repmat([90,-90],1,nBinauralStimuli/nSide);
aziOptions = [90,-90];

%% Preset values for diotic stimuli
% No. of stimuli and conditons
nDdi = 1; % unique distance values for diotic stimuli
nDioticStimuli = nF0*nDdi*nSide;

%% Input check and defaults
% Duration
if ~exist('dur','var')
    dur = 2;
end

% F0
if ~exist('f0', 'var')
    fmin = 310;
    fmax = 2*fmin; % 1 octave
    fOptions = round(linspace(fmin,fmax,nF0),-1);
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

if exist('f0','var') % if constant f0 is requested
    wavDir = strcat(wavDir,'_',num2str(f0));
end
mkdir(wavDir);


%% Create a home for saved parameters
outFields = {
    'filename', ...
    'f0', ...                % F0 in Hz
    'dur', ...               % Duration in s
    'distance', ...          % Distance in m
    'azimuth', ...           % Side: 90 - left, -90 - right
    'fs', ...                % sampling rate
    };

outCsv = cell(nStimuli+nDioticStimuli+1,length(outFields));
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

stimNo = 1;
%% Stimulus generation loop for dichotic sounds (BRT)
for ff = 1:length(fOptions)
    % Randomize f0 (picked from previously set unique values)
    f0 = fOptions(ff);

    %% Timing variables
    buffer = 1; % 1 s buffer for BRT artefacts

    %% Position variables
    for dd = 1:length(dOptions)
        for aa = 1:length(aziOptions)

            % Progress update
            clc;
            disp(['Processing stimulus ',num2str(stimNo),'/',num2str(nStimuli),'...']);

            azi = aziOptions(aa);
            ele = 0; % Elevation always 0
            r = dOptions(dd);

            %% Generate triangle wave
            t = linspace(0,dur,dur*fs); % Time vector from 0 to dur with sampling interval 1/fs
            triwave = sawtooth(2 * pi * f0 * t, 0.5); % 0.5 for a symmetric triangular wave

            if r > in % if spatialized stimuli
                    %% Prepare for spatialization
                    win = tukeywin(size(triwave,2),(0.1*fs)/size(triwave,2)); % 0.1 s on-offset ramp
                    triwave = triwave.*win';

                    % Save as wav (BRT needs it)
                    digits = ceil(log10(nStimuli + 1));
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
                    oscsend(u, '/source/loadSource', 'sss', 'source', wavname, 'DirectivityModel'); WaitSecs(1);
                    [x,y,z] = sph2cart(deg2rad(azi),deg2rad(ele),r);
                    oscsend(u, '/source/location', 'sfff', 'source', x, y, z); WaitSecs(5);
                    oscsend(u, '/source/playAndRecord', 'sssf', 'source', savename, 'mat', dur+buffer); WaitSecs(3);
                    oscsend(u, '/stop', 's', ''); WaitSecs(3);

                    %% Create output
                    params.sofa = load(savename);
                    spat = params.sofa.Data.Receiver;

                    % Cut off the onset and offset artefacts introduced by BRT during recording;
                    % apply Tukey window on the sound
                    spat = spat(:,1:round(dur*fs));
                    win = tukeywin(size(spat,2),(0.1*fs)/size(spat,2)); % 0.1 s
                    spatWin = spat.*repmat(win',2,1);

                    % Update duration data to exclude buffer
                    totalDur = length(spatWin)/fs;

                    % Scale sound to insure intensity differences
                    scaleto = PPS/r;
                    scaledb = 20 * log10(scaleto / max(spatWin,[],"all"));

                    scaled1 = scaletodbspl(spatWin(1,:)',dbspl(spatWin(1,:)')+scaledb);
                    scaled2 = scaletodbspl(spatWin(2,:)',dbspl(spatWin(2,:)')+scaledb);
                    out = [scaled1 scaled2];

                    % Save stimulus in .mat file
                    digits = ceil(log10(nStimuli + 1));
                    stimulusdigits = ceil(log10(stimNo + 1));
                    temp = char('');
                    for digind = 1:(digits-stimulusdigits)
                        temp = strcat(temp, '0');
                    end
                    filename = strcat(wavDir, '/', wavDir, '-', temp, num2str(stimNo));
                    save(strcat(filename,'.mat'), "out", "fs");
                    WaitSecs(4); % avoid BRT crashing

            elseif r == in
                    %% Window and stereo
                    win = tukeywin(size(triwave,2),(0.1*fs)/size(triwave,2)); % 0.1 s on-offset ramp
                    triwave = repmat(triwave.*win',2,1);

                    % Scale sound 
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
                    filename = strcat(wavDir, '/', wavDir, '-', temp, num2str(stimNo));
                    save(strcat(filename,'.mat'), "out", "fs");
                    totalDur = length(out)/fs;
            end

            % Add parameters to cell array for later saving
            outCsv(stimNo+1, :) = {filename, f0, totalDur, r, azi, fs};
            stimStruct.fs = fs;
            stimStruct.stim{ff,dd,aa} = {out, r};

            stimNo = stimNo+1;
        end % azi options
    end % distance options
end % f0 options // stimulus generation loop for BRT

%% Save output
% Convert cell to a table and use first row as variable names
T = cell2table(outCsv(2:end, :), 'VariableNames', outCsv(1, :));

% Write the table to a CSV file, final user message
writetable(T,strcat('./', wavDir, '/', strcat(wavDir, '-', 'StimuliData.csv')));
save("stimStruct.mat","stimStruct");

disp([newline, 'Task done, files and parameters are saved to directory ', wavDir, newline]);
end
