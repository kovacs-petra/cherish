function stimStruct = genStimExt(dur,f0)
% Usage: genStimExt(dur)
%
% Generates sounds of dur s duration, at different f0s (in one
% octave range), distances, and azimuths. The stimuli are generated
% with the Binaural Rendering Toolbox (BRT).
%
% Inputs:
% - dur:   duration of each sound in s (default: 2 s)
% - f0:    set constant f0 for all stimuli (default: varies in one octave range)
%
% Requirements:
% - BRT Renderer App (needs to be open)
% - oscsend.m (insert path in the script)
%
% #Author: Petra Kovacs, 2025, Acoustics Research Institute, Vienna,
% Austria
%

%% Input check and defaults
% Duration
if ~exist('dur','var')
    dur = 2;
end

% F0
if ~exist('f0', 'var')
    nF0 = 8;
    fmin = 310;
    fmax = 2*fmin; % 1 octave
    fOptions = round(linspace(fmin,fmax,nF0),-1);
else
    fOptions = f0;
    nF0 = 1;
end

nD = 6; % unique distance values
nSide = 2; % unique azimuth values
nStimCond = 2; % stimuli type conditions
nStimuli = nF0*nD*nSide*nStimCond;
stimStruct.stim = cell(nF0,nD,nSide);

% Distances (m)
% in = 0;
PPS = 0.2;
EPS = 2;
dOptions = [PPS,EPS,round(linspace(PPS+0.2,EPS-0.2,nD-2),1)];
                                                            % 1st two columns: the fixed
                                                            % distances we want to use in the
                                                            % experiment.
                                                            % Then, the distances inbetween.
% Azimuths (deg)
aziOptions = [90,-90];

% Stimuli type conditions
stimTypeOptions = [1,2]; % triangle wave, white noise

% Path of oscsend needed
oscsend_path = 'C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusGeneration';
addpath(oscsend_path);

%% Initialize AMT
if ~exist("amt_start.m","file")
    amtPath = '\\KFS\Fileserver\ProjektDaten\CherISH\code\amtoolbox-full-1.4.0';
    addpath(amtPath);
    amt_start;
end

%% Load an HRTF dataset to determine the right sampling rate
% if ~exist("SOFAdbPath.m","file")
%     sofaPath = '\\kfs\fileserver\ProjektDaten\CherISH\code\SOFAtoolbox\SOFAtoolbox';
%     addpath(sofaPath);
%     SOFAstart;
% end
% database = 'sadie';
% HRTFfilename = '3DTI_HRTF_SADIE_II_D2_256s_48000Hz_resampled5.sofa';
% fullfn = fullfile(SOFAdbPath, 'database', database, HRTFfilename);
% Obj = SOFAload(fullfn);
% fs = Obj.Data.SamplingRate;
fs = 48e3;

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
% outFields = {
%     'filename', ...
%     'f0', ...                % F0 in Hz
%     'dur', ...               % Duration in s
%     'distance', ...          % Distance in m
%     'azimuth', ...           % Side: 90 - left, -90 - right
%     'stimType',...           % Stimulus type: 1 - triangle wave, 2 - white noise
%     'fs', ...                % sampling rate
%     };
% 
% outCsv = cell(nStimuli+1,length(outFields));
% outCsv(1,:) = outFields;

%% Set connection between Matlab and BRT through OSC; set some settings
u = pnet('udpsocket',10017); % Listen port in BRT
oscsend(u, '/control/connect', 'si', 'localhost',10019);
oscsend(u, '/listener/enableSpatialization', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableInterpolation', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableNearFieldEffect', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableITD', 'sB', 'DefaultListener',1);
oscsend(u, '/enableModel', 'sB', 'FreeField',0);
oscsend(u, '/enableModel', 'sB', 'ReverbPath',1);
oscsend(u, '/environment/enableDirectPath', 'sB', 'SDN',1);
oscsend(u, '/environment/enableReverbPath', 'sB', 'SDN',1);
oscsend(u, '/listener/enableParallaxCorrection', 'sB', 'DefaultListener',1);
oscsend(u, '/enableModel', 'sB', 'DirectPath', 1);
% oscsend(u, '/environment/enableDistanceAttenuation', 'sB', 'FreeField', 1);
oscsend(u, '/environment/enableDistanceAttenuation', 'sB', 'SDN', 1);

buffer = 1; % buffer size in s for BRT artefacts
stimNo = 1;
aa = 1;
ele = 0; % Elevation always 0

%% Stimulus generation loop
for tt = 1:length(stimTypeOptions) % Stimulus type
    stimType = stimTypeOptions(tt);

    for ff = 1:length(fOptions) % F0
        f0 = fOptions(ff);

        for dd = 1:length(dOptions) % Distance
            r = dOptions(dd);

            azi = aziOptions(aa);


            % Progress update
            clc; disp(['Processing stimulus ',num2str(stimNo),'/',num2str(nStimuli/2),'...',...
                newline,'f0 = ',num2str(f0),'; r = ',num2str(r),'; stimType = ', ...
                num2str(stimType)]);

            switch stimType
                case 1 % Generate triangle wave with on-offset ramp
                    t = linspace(0,dur+buffer,(dur+buffer)*fs); % Time vector from 0 to dur with sampling interval 1/fs
                    signal = sawtooth(2 * pi * f0 * t, 0.5); % 0.5 for a symmetric triangular wave
                    % win = tukeywin(size(triwave,2),(0.1*fs)/size(triwave,2)); % 0.1 s on-offset ramp
                    % signal = triwave.*win';

                case 2 % Generate white noise with on-offset ramp
                    signal = (randn(fs*(dur+buffer),1)+1)';
                    % win = tukeywin(size(wn,2),(0.1*fs)/size(wn,2)); % 0.1 s on-offset ramp
                    % signal = wn.*win';
            end

            % Save signal as wav (BRT needs it)
            digits = ceil(log10(nStimuli + 1));
            stimulusdigits = ceil(log10(stimNo + 1));
            temp = char('');
            for digind = 1:(digits-stimulusdigits)
                temp = strcat(temp, '0');
            end
            filename = fullfile(pwd, strcat(wavDir, '/', wavDir, '-', temp, num2str(stimNo)));
            wavname = strcat(filename,'.wav');
            audiowrite(wavname, signal, fs);
            savename = strcat(filename,'.mat');

            %% Spatialize with BRT
            oscsend(u,'/removeAllSources', 's', ''); WaitSecs(2);
            oscsend(u, '/source/loadSource', 'sss', 'source', wavname, 'DirectivityModel'); WaitSecs(3);
            [x,y,z] = sph2cart(deg2rad(azi),deg2rad(ele),r);
            oscsend(u, '/source/location', 'sfff', 'source', x, y, z); WaitSecs(3);
            oscsend(u, '/source/playAndRecord', 'sssf', 'source', savename, 'mat', dur+buffer*3); WaitSecs(6);
            oscsend(u, '/stop', 's', ''); WaitSecs(1);

            %% Create output
            params.sofa = load(savename);
            spat = params.sofa.Data.Receiver;

            % Cut off the onset and offset artefacts introduced by BRT during recording;
            % apply Tukey window on the sound
            % spat = spat(:,buffer*fs:round(dur*fs));
            spat = spat(:,(buffer/2*fs)+1:(dur+buffer/2)*fs);
            win = tukeywin(size(spat,2),(0.1*fs)/size(spat,2)); % 0.1 s
            spatWin = spat.*repmat(win',2,1);

            % Update duration data to exclude buffer
            disp(['totalDur = ', num2str(length(spatWin)/fs)]);

            % Scale sound to ensure intensity differences
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
            WaitSecs(2); % avoid BRT crashing

            % Add parameters to cell array for later saving
            stimStruct.fs = fs;
            stimStruct.stim{ff,dd,aa,tt} = {out,r};
            stimStruct.stim{ff,dd,aa+1,tt} = {flip(out,2),r}; % flip symmetrically for other side
            % outCsv(stimNo+1, :) = {filename, f0, totalDur, r, azi, stimType, fs};

            stimNo = stimNo+1;

            figure; plot(out); drawnow; WaitSecs(1); % Plot figure for checking
            save("stimStruct.mat","stimStruct");

        end % distance options
    end % f0 options
end % stimType options // stimulus generation loop

%% Save output
% Convert cell to a table and use first row as variable names
% T = cell2table(outCsv(2:end, :), 'VariableNames', outCsv(1, :));

% Write the table to a CSV file, final user message
% writetable(T,strcat('./', wavDir, '/', strcat(wavDir, '-', 'StimuliData.csv')));
% save("stimStruct.mat","stimStruct");

disp([newline, 'Task done, files and parameters are saved to directory ', wavDir, newline]);
end
