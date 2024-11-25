function generateStimuli(nBuildingBlocks,whichBlock,trajectory)
% Generate stimuli for wp1 pilot: looming and receding sounds in duplets
% Address complaints to: Petra Kovacs

% Inputs:
%% trajectory
% trajectory = 1; % loom-rec (block 1-2) // 0 to 30 deg (block 3)
% trajectory = 2; % rec-loom (block 1-2) // 30 to 0 deg (block 3)
% trajectory = 3; % loom-loom (block 2)
% trajectory = 4; % rec-rec (block 2)

%% no. of stimuli in each distance x direction pair (termed a bulding block):
% Note that the total no. of stimuli per block will be 18*nBuildingBlocks in block 1 and
% 36*nStimuli in block 2. % See also nStimuli.xlsx for an illustration.
% nBuildingBlocks = 1;
stimuliPerTrajectory = nBuildingBlocks*9;

%% block (1:3):
% whichBlock = 1; % BIG block (task block)

%% Parameters and to-do independent of block number

ele = [0 0]; % elevation is always 0

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
PPS = [0.05 0.1 0.15];
ARS = [0.8 0.9 1.0];
EPS = [2.5 3.0 3.5];
allSpace = [PPS ARS EPS];
v = 1; % m/s

%% Create a home for saved parameters
if whichBlock ~= 3
    outFields = {
        'filename', ...         % File data
        'frequency',...
        'durStatStart',...      % Duration data
        'durMov1',...
        'durStatMiddle',...
        'durMov2',...
        'durStatEnd',...
        'totalDur',...
        'startDistMov1',...     % Space data - distance in m
        'stopDistMov1',...
        'startDistMov2',...
        'stopDistMov2',...
        'startSpaceMov1',...    % Space data - space category (PPS,ARS,EPS)
        'stopSpaceMov1',...
        'startSpaceMov2',...
        'stopSpaceMov2',...
        'trajectory'            % Space data - direction (loom,rec)
        };

    outCsv = cell(stimuliPerTrajectory+1,length(outFields));
    outCsv(1,:) = outFields;

else
    outFields = {
        'filename', ...    % File data
        'frequency',...
        'durStatStart',... % Duration data
        'durMov1',...
        'durStatMiddle',...
        'durMov2',...
        'durStatEnd',...
        'totalDur',...
        'distance',...     % Space data - distance in m
        'whichSpace',...   % Space data - space category (PPS,ARS,EPS)
        'trajectory'       % Space data - direction (0-30deg or umgekehrt)
        };

    outCsv = cell(stimuliPerTrajectory+1,length(outFields));
    outCsv(1,:) = outFields;
end

%% Preset position of stationary portions for all stimuli based on trajectory
% Blocks 1 and 2: x, y and z denote distances (radii)
if whichBlock ~= 3
    switch trajectory
        case 1  % loom-rec: x > y, y < z, z ~= 0.1
            z = [...
                repmat(PPS(2),1,nBuildingBlocks),...
                repmat(PPS(3),1,2*nBuildingBlocks),...
                repmat(ARS(1),1,nBuildingBlocks),...
                repmat(ARS(2),1,nBuildingBlocks),...
                repmat(ARS(3),1,nBuildingBlocks),...
                repmat(EPS(1),1,nBuildingBlocks),...
                repmat(EPS(2),1,nBuildingBlocks),...
                repmat(EPS(3),1,nBuildingBlocks),...
                ];
            y = zeros(1,length(z));
            x = zeros(1,length(z));

            for i = 1: length(z)
                possibleSpaceY = allSpace(allSpace < z(i));
                idx = randi(length(possibleSpaceY),1);
                y(i) = possibleSpaceY(idx);
            end

            for i = 1:length(y)
                possibleSpaceX = allSpace(allSpace > y(i));
                idx = randi(length(possibleSpaceX),1);
                x(i) = possibleSpaceX(idx);
            end

        case 2 % rec-loom: x < y, y > z, z ~= 2.1
            z = [...
                repmat(PPS(1),1,nBuildingBlocks),...
                repmat(PPS(2),1,nBuildingBlocks),...
                repmat(PPS(3),1,nBuildingBlocks),...
                repmat(ARS(1),1,nBuildingBlocks),...
                repmat(ARS(2),1,nBuildingBlocks),...
                repmat(ARS(3),1,nBuildingBlocks),...
                repmat(EPS(1),1,2*nBuildingBlocks),...
                repmat(EPS(2),1,nBuildingBlocks),...
                ];
            y = zeros(1,length(z)); % has to be greater than z
            x = zeros(1,length(z)); % has to be smaller than y

            for i = 1:length(z)
                possibleSpaceY = allSpace(allSpace > z(i));
                idx = randi(length(possibleSpaceY),1);
                y(i) = possibleSpaceY(idx);
            end

            for i = 1:length(y)
                possibleSpaceX = allSpace(allSpace < y(i));
                idx = randi(length(possibleSpaceX),1);
                x(i) = possibleSpaceX(idx);
            end

        case 3 % loom-loom: x > y, y > z, z ~= 2.1
            z = [...
                repmat(PPS(1),1,nBuildingBlocks),...
                repmat(PPS(2),1,nBuildingBlocks),...
                repmat(PPS(3),1,nBuildingBlocks),...
                repmat(ARS(1),1,nBuildingBlocks),...
                repmat(ARS(2),1,nBuildingBlocks),...
                repmat(ARS(3),1,nBuildingBlocks),...
                repmat(EPS(1),1,3*nBuildingBlocks),...
                ];
            y = zeros(1,length(z));
            x = zeros(1,length(z));

            for i = 1:length(z)
                possibleSpaceY = allSpace(allSpace > z(i));
                possibleSpaceY = possibleSpaceY(possibleSpaceY < EPS(3));
                idx = randi(length(possibleSpaceY),1);
                y(i) = possibleSpaceY(idx);
            end

            for i = 1:length(y)
                possibleSpaceX = allSpace(allSpace > y(i));
                idx = randi(length(possibleSpaceX),1);
                x(i) = possibleSpaceX(idx);
            end

        case 4 % rec-rec
            z = [...
                repmat(PPS(3),1,3*nBuildingBlocks),...
                repmat(ARS(1),1,nBuildingBlocks),...
                repmat(ARS(2),1,nBuildingBlocks),...
                repmat(ARS(3),1,nBuildingBlocks),...
                repmat(EPS(1),1,nBuildingBlocks),...
                repmat(EPS(2),1,nBuildingBlocks),...
                repmat(EPS(3),1,nBuildingBlocks),...
                ];
            y = zeros(1,length(z));
            x = zeros(1,length(z));

            for i = 1:length(z)
                possibleSpaceY = allSpace(allSpace < z(i));
                possibleSpaceY = possibleSpaceY(possibleSpaceY > PPS(1));
                idx = randi(length(possibleSpaceY),1);
                y(i) = possibleSpaceY(idx);
            end

            for i = 1:length(y)
                possibleSpaceX = allSpace(allSpace < y(i));
                idx = randi(length(possibleSpaceX),1);
                x(i) = possibleSpaceX(idx);
            end

    end % switch trajectory blocks 1-2

elseif whichBlock == 3 % x, y and z denote azimuth (set later). R is set here:
    rBlock3 = [...
        repmat(PPS(1),1,nBuildingBlocks),...
        repmat(PPS(2),1,nBuildingBlocks),...
        repmat(PPS(3),1,nBuildingBlocks),...
        repmat(ARS(1),1,nBuildingBlocks),...
        repmat(ARS(2),1,nBuildingBlocks),...
        repmat(ARS(3),1,nBuildingBlocks),...
        repmat(EPS(1),1,nBuildingBlocks),...
        repmat(EPS(2),1,nBuildingBlocks),...
        repmat(EPS(3),1,nBuildingBlocks)...
        ];
   
end % if not block 3

%% Stimulus generation loop
if whichBlock ~= 3
    azi = [30 30]; % azi is fixed in the first two blocks

    for stimNo = 1:stimuliPerTrajectory
        frequency = (randi(4,1)+4)*100; % 500-800 Hz with round 100 values

        
        %% Radii for this stimulus
        rStatStart = [x(stimNo) x(stimNo)];
        rStatMiddle = [y(stimNo) y(stimNo)];
        rStatEnd = [z(stimNo) z(stimNo)];

        rMov1 = [x(stimNo) y(stimNo)];
        rMov2 = [y(stimNo) z(stimNo)];

        % rMain = [rStatStart,rMov1,rStatMiddle,rMov2,rStatEnd];

        % Save which space category the moving portions begin and end in
        % Start of moving portion 1
        if ismember(rMov1(1),PPS)
            startSpaceMov1 = 1;
        elseif ismember(rMov1(1),ARS)
            startSpaceMov1 = 2;
        else
            startSpaceMov1 = 3;
        end

        % Stop of moving portion 1
        if ismember(rMov1(2),PPS)
            stopSpaceMov1 = 1;
        elseif ismember(rMov1(2),ARS)
            stopSpaceMov1 = 2;
        else
            stopSpaceMov1 = 3;
        end

        % Start of moving portion 2
        if ismember(rMov2(1),PPS)
            startSpaceMov2 = 1;
        elseif ismember(rMov2(1),ARS)
            startSpaceMov2 = 2;
        else
            startSpaceMov2 = 3;
        end

        % Stop of moving portion 2
        if ismember(rMov2(2),PPS)
            stopSpaceMov2 = 1;
        elseif ismember(rMov2(2),ARS)
            stopSpaceMov2 = 2;
        else
            stopSpaceMov2 = 3;
        end

        %% Set durations
        % Static portions
        durStatStart = (randi(10,1)+10)/10; % 1000-2000 ms with round 100 values
        durStatMiddle = (randi(10,1)+10)/10;
        durStatEnd = (randi(10,1)+10)/10;

        % Moving portions
        durMov1 = abs(rMov1(2)-rMov1(1))/v;
        durMov2 = abs(rMov2(2)-rMov2(1))/v;

        durations = struct( ...
            'durStatStart',durStatStart, ...
            'durMov1',durMov1, ...
            'durStatMiddle',durStatMiddle, ...
            'durMov2',durMov2, ...
            'durStatEnd', durStatEnd);

        rMain = [linspace(rStatStart(1),rStatStart(2),durStatStart*fs),...
            linspace(rMov1(1),rMov1(2),durMov1*fs),...
            linspace(rStatMiddle(1),rStatMiddle(2),durStatMiddle*fs),...
            linspace(rMov2(1),rMov2(2),durMov2*fs),...
            linspace(rStatEnd(1),rStatEnd(2),durStatEnd*fs)];


        %% Generate square waves and intensity ramps for each portion

        % Generate square waves for stationary portions
        tStatStart = 0:1/fs:durations.durStatStart;
        statStart = square(2*pi*frequency*tStatStart);

        tStatMiddle = 0:1/fs:durations.durStatMiddle;
        statMiddle = square(2*pi*frequency*tStatMiddle);

        tStatEnd = 0:1/fs:durations.durStatEnd;
        statEnd  = square(2*pi*frequency*tStatEnd);

        % Create corresponding intensity ramps
        intensityRampStart = 1./linspace(rStatStart(1),rStatStart(2),length(statStart));
        intensityRampMiddle = 1./linspace(rStatMiddle(1),rStatMiddle(2),length(statMiddle));
        intensityRampEnd = 1./linspace(rStatEnd(1),rStatEnd(2),length(statEnd));

        %
        % Generate square waves for moving portions
        tMov1 = 0:1/fs:durations.durMov1;
        mov1 = square(2*pi*frequency*tMov1);

        tMov2 = 0:1/fs:durations.durMov2;
        mov2 = square(2*pi*frequency*tMov2);

        % Create corresponding intensity ramps
        intensityRampMov1 = 1./linspace(rMov1(1),rMov1(2),length(mov1));
        intensityRampMov2 = 1./linspace(rMov2(1),rMov2(2),length(mov2));

        %% Concatenate and spatialize
        % intensityRampMain = [intensityRampStart intensityRampMov1 intensityRampMiddle intensityRampMov2 intensityRampEnd];
        allPortions = [statStart mov1 statMiddle mov2 statEnd];
        [stim, ~, ~, rActual, idx] = local_SOFAspat(allPortions',Obj,azi,ele,rMain);
        stim = rescale(stim,-0.8,0.8);

        totalDur = durMov1+durMov2+durStatStart+durStatMiddle+durStatEnd; % in s

        % Save results to wav, add parameters to the cell array later saved out
        % to csv
        digits = ceil(log10(stimuliPerTrajectory + 1));
        stimulusdigits = ceil(log10(stimNo + 1));
        temp = char('');
        for digind = 1:(digits-stimulusdigits)
            temp = strcat(temp, '0');
        end

        filename = strcat(wavDir, '-', temp, num2str(stimNo));
        outCsv(stimNo+1, :) = {filename, frequency, durations.durStatStart, durations.durMov1, durations.durStatMiddle, ...
            durations.durMov2, durations.durStatEnd, totalDur, rMov1(1), ...
            rMov1(2),rMov2(1),rMov2(2), startSpaceMov1, stopSpaceMov1, startSpaceMov2,...
            stopSpaceMov2, trajectory};
        audiowrite(strcat('./', wavDir, '/', filename, '.wav'), stim, fs);

    end % for stimNo

elseif whichBlock == 3 % Azimuth block
    for stimNo = 1:stimuliPerTrajectory
        frequency = (randi(4,1)+4)*100; % 500-800 Hz with round 100 values

        %% Set azimuth:
        switch trajectory
            case 1
                x = 0;
                y = 30;
                z = 0;
            case 2
                x = 30;
                y = 0;
                z = 30;
        end % switch trajectory block 3

        aziStatStart = [x x];
        aziStatMiddle = [y y];
        aziStatEnd = [z z];

        aziMov1 = [x y];
        aziMov2 = [y z];

        aziMain = [aziStatStart,aziMov1,aziStatMiddle,aziMov2,aziStatEnd];

        %% Set durations
        % Static portions: 800-1100 ms, with round 100 ms values, in secs
        durStatStart = (randi(4,1)+7)/10;
        durStatMiddle = (randi(4,1)+7)/10;
        durStatEnd = (randi(4,1)+7)/10;

        % Moving portions: 1000-2000 ms, with round 100 ms values, in secs
        durMov1 = (randi(10,1)+10)/10;
        durMov2 = (randi(10,1)+10)/10;

        durations = struct( ...
            'durStatStart',durStatStart, ...
            'durMov1',durMov1, ...
            'durStatMiddle',durStatMiddle, ...
            'durMov2',durMov2, ...
            'durStatEnd', durStatEnd);

        %% Generate square waves and intensity ramps for each portion

        % Generate square waves for stationary portions
        tStatStart = 0:1/fs:durations.durStatStart;
        statStart = square(2*pi*frequency*tStatStart);

        tStatMiddle = 0:1/fs:durations.durStatMiddle;
        statMiddle = square(2*pi*frequency*tStatMiddle);

        tStatEnd = 0:1/fs:durations.durStatEnd;
        statEnd  = square(2*pi*frequency*tStatEnd);

        %
        % Generate square waves for moving portions
        tMov1 = 0:1/fs:durations.durMov1;
        mov1 = square(2*pi*frequency*tMov1);

        tMov2 = 0:1/fs:durations.durMov2;
        mov2 = square(2*pi*frequency*tMov2);

        %% Concatenate and spatialize
        allPortions = [statStart mov1 statMiddle mov2 statEnd];
        stim = local_SOFAspat(allPortions',Obj,aziMain,ele,rBlock3(stimNo));

        stim = rescale(stim,-1,1);

        % Save which space category this stim belongs to
        if ismember(rBlock3(stimNo),PPS)
            whichSpace = 1;
        elseif ismember(rBlock3(stimNo),ARS)
            whichSpace = 2;
        else
            whichSpace = 3;
        end
        
        totalDur = durMov1+durMov2+durStatStart+durStatMiddle+durStatEnd; % in s

        % Save results to wav, add parameters to the cell array later saved out
        % to csv
        digits = ceil(log10(stimuliPerTrajectory + 1));
        stimulusdigits = ceil(log10(stimNo + 1));
        temp = char('');
        for digind = 1:(digits-stimulusdigits)
            temp = strcat(temp, '0');
        end

        filename = strcat(wavDir, '-', temp, num2str(stimNo));
        outCsv(stimNo+1, :) = {filename, frequency, durations.durStatStart, durations.durMov1, durations.durStatMiddle, ...
            durations.durMov2, durations.durStatEnd, totalDur, rBlock3(stimNo), whichSpace, trajectory};
        audiowrite(strcat('./', wavDir, '/', filename, '.wav'), stim, fs);

    end % stimulus generation loop for block 3


end % if blockNo

% Convert cell to a table and use first row as variable names
T = cell2table(outCsv(2:end, :), 'VariableNames', outCsv(1, :));

% Write the table to a CSV file, final user message
writetable(T,strcat('./', wavDir, '/', strcat(wavDir, '-', 'StimuliData.csv')));
disp([newline, 'Task done, files and parameters are saved to directory ', wavDir, newline]);

end % function

function [out, aziActual, eleActual, rActual, idx] = local_SOFAspat(signal,Obj,azi,ele,r)
if length(r)~=length(signal)
    error('Error: signal and r trajectory need to have the same length')
end

[out, aziActual, eleActual, rActual, idx] = SOFAspat(signal./r,Obj,azi,ele,r);
end
