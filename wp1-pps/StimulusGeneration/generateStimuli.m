function generateStimuli(nBuildingBlocks,whichBlock,trajectory)
% Generate stimuli for wp1 pilot: looming and receding sounds in duplets
% Address complaints to: Petra Kovacs

% Inputs:
%% nBuildingBlocks: no. of stimuli in each distance x direction pair:
% Note that the total no. of stimuli per block will be 18*nBuildingBlocks in block 1 and
% 36*nBuildingBlocks in block 2. See also nStimuli.xlsx for an illustration.
stimuliPerTrajectory = nBuildingBlocks*9;

% Sanity check: nBuildingBlocks has to be an even no. so that half of the
% trials can have a target
if mod(nBuildingBlocks,2) ~= 0
    error('nBuildingBlocks has to be an even number');
end

%% whichBlock (1:3):
% Big block (task block)
% 1 - loom-rec and rec-loom stimuli only
% 2 - besides the conditions in 1, there's also loom-loom and rec-rec
% conditions
% 3 - sound is fixed in distance, moves left or right

%% trajectory
% 1 - loom-rec (block 1-2) // 0 to 30 deg i.e. left (block 3)
% 2 - rec-loom (block 1-2) // 0 to -30 deg i.e. right (block 3)
% 3 - loom-loom (block 2)
% 4 - rec-rec (block 2)

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
v = 0.5; % velocity in m/s

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
        'trajectory',...        % Space data - direction (loom, rec)
        'azimuth',...
        'target'                % Target or nontarget trial
        };

    outCsv = cell(stimuliPerTrajectory+1,length(outFields));
    outCsv(1,:) = outFields;

else
    outFields = {
        'filename', ...    % File data
        'frequency',...
        'durStatStart',... % Duration data
        'durMov1',...
        'durStatEnd',...
        'totalDur',...
        'distance',...     % Space data - distance in m
        'whichSpace',...   % Space data - space category (PPS,ARS,EPS)
        'trajectory'       % Space data - direction (0-30 deg or umgekehrt)
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

            % Allocate target and nontarget trials:
            % Half of the unique end distances in z have to be target, half nontarget
            % Follow the logic of variable z, but alternate 0s and 1s at half the
            % frequency
            target = [...
                ones(1,nBuildingBlocks/2),...   % PPS(2) - 1st half of trials
                zeros(1,nBuildingBlocks/2),...  % PPS(2) - 2nd half of trials
                ones(1,nBuildingBlocks),...     % PPS(3)
                zeros(1,nBuildingBlocks),...    % PPS(3)
                ones(1,nBuildingBlocks/2),...   % ARS(1)
                zeros(1,nBuildingBlocks/2),...  % ARS(1)
                ones(1,nBuildingBlocks/2),...   % ARS(2)
                zeros(1,nBuildingBlocks/2),...  % ARS(2)
                ones(1,nBuildingBlocks/2),...   % ARS(3)
                zeros(1,nBuildingBlocks/2),...  % ARS(3)
                ones(1,nBuildingBlocks/2),...   % EPS(1)
                zeros(1,nBuildingBlocks/2),...  % EPS(1)
                ones(1,nBuildingBlocks/2),...   % EPS(2)
                zeros(1,nBuildingBlocks/2),...  % EPS(2)
                ones(1,nBuildingBlocks/2),...   % EPS(3)
                zeros(1,nBuildingBlocks/2)...   % EPS(3)
                ];

            % x and y coordinates:
            y = zeros(1,length(z)); % has to be smaller than z
            x = zeros(1,length(z)); % has to be greater than y
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

            % Allocate target and nontarget trials:
            target = [...
                ones(1,nBuildingBlocks/2),...   % PPS(1) - 1st half of trials
                zeros(1,nBuildingBlocks/2),...  % PPS(1) - 2nd half of trials
                ones(1,nBuildingBlocks/2),...   % PPS(2)
                zeros(1,nBuildingBlocks/2),...  % PPS(2)
                ones(1,nBuildingBlocks/2),...   % PPS(3)
                zeros(1,nBuildingBlocks/2),...  % PPS(3)
                ones(1,nBuildingBlocks/2),...   % ARS(1)
                zeros(1,nBuildingBlocks/2),...  % ARS(1)
                ones(1,nBuildingBlocks/2),...   % ARS(2)
                zeros(1,nBuildingBlocks/2),...  % ARS(2)
                ones(1,nBuildingBlocks/2),...   % ARS(3)
                zeros(1,nBuildingBlocks/2),...  % ARS(3)
                ones(1,nBuildingBlocks),...     % EPS(1)
                zeros(1,nBuildingBlocks),...    % EPS(1)
                ones(1,nBuildingBlocks/2),...   % EPS(2)
                zeros(1,nBuildingBlocks/2)...   % EPS(2)
                ];

            % x and y coordinates:
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

             % Allocate target and nontarget trials:
             target = [...
                ones(1,nBuildingBlocks/2),...   % PPS(1) - 1st half of trials
                zeros(1,nBuildingBlocks/2),...  % PPS(1) - 2nd half of trials
                ones(1,nBuildingBlocks/2),...   % PPS(2)
                zeros(1,nBuildingBlocks/2),...  % PPS(2)
                ones(1,nBuildingBlocks/2),...   % PPS(3)
                zeros(1,nBuildingBlocks/2),...  % PPS(3)
                ones(1,nBuildingBlocks/2),...   % ARS(1)
                zeros(1,nBuildingBlocks/2),...  % ARS(1)
                ones(1,nBuildingBlocks/2),...   % ARS(2)
                zeros(1,nBuildingBlocks/2),...  % ARS(2)
                ones(1,nBuildingBlocks/2),...   % ARS(3)
                zeros(1,nBuildingBlocks/2),...  % ARS(3)
                ones(1,3*nBuildingBlocks/2),... % EPS(1)
                zeros(1,3*nBuildingBlocks/2),...% EPS(1)
                ];

            % x and y coordinates:
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

            % Allocate target and nontarget trials:
             target = [...
                ones(1,3*nBuildingBlocks/2),... % PPS(3) - 1st half of trials
                zeros(1,3*nBuildingBlocks/2),...% PPS(3) - 2nd half of trials
                ones(1,nBuildingBlocks/2),...   % ARS(1)
                zeros(1,nBuildingBlocks/2),...  % ARS(1)
                ones(1,nBuildingBlocks/2),...   % ARS(2)
                zeros(1,nBuildingBlocks/2),...  % ARS(2)
                ones(1,nBuildingBlocks/2),...   % ARS(3)
                zeros(1,nBuildingBlocks/2),...  % ARS(3)
                ones(1,nBuildingBlocks/2),...   % EPS(1)
                zeros(1,nBuildingBlocks/2),...  % EPS(1)
                ones(1,nBuildingBlocks/2),...   % EPS(2)
                zeros(1,nBuildingBlocks/2),...  % EPS(2)
                ones(1,nBuildingBlocks/2),...   % EPS(3)
                zeros(1,nBuildingBlocks/2),...  % EPS(3)
                ];

            % x and y coordinates: 
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

elseif whichBlock == 3 % x and y denote azimuth (set later). R is set here:
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
    azi = repmat([-30 30],1,stimuliPerTrajectory/2);

    for stimNo = 1:stimuliPerTrajectory
        frequency = (randi(4,1)+4)*100; % 500-800 Hz with round 100 values

        
        %% Radii for this stimulus
        rStatStart = [x(stimNo) x(stimNo)];
        rStatMiddle = [y(stimNo) y(stimNo)];
        rStatEnd = [z(stimNo) z(stimNo)];

        rMov1 = [x(stimNo) y(stimNo)];
        rMov2 = [y(stimNo) z(stimNo)];

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
        durStatStart = (randi(4,1)+5)/10; % 600-900 ms with round 100 values
        durStatMiddle = (randi(4,1)+5)/10;
        durStatEnd = (randi(4,1)+5)/10;

        % Moving portions: duration calculated from distance and velocity
        durMov1 = abs(rMov1(2)-rMov1(1))/v;
        durMov2 = abs(rMov2(2)-rMov2(1))/v;

        rMain = [linspace(rStatStart(1),rStatStart(2),durStatStart*fs),...
            linspace(rMov1(1),rMov1(2),durMov1*fs),...
            linspace(rStatMiddle(1),rStatMiddle(2),durStatMiddle*fs),...
            linspace(rMov2(1),rMov2(2),durMov2*fs),...
            linspace(rStatEnd(1),rStatEnd(2),durStatEnd*fs)];


        %% Generate square waves and intensity ramps for each portion

        % Generate square waves for stationary portions
        tStatStart = linspace(0,durStatStart,durStatStart*fs);
        statStart = square(2*pi*frequency*tStatStart);

        tStatMiddle = linspace(0,durStatMiddle,durStatMiddle*fs);
        statMiddle = square(2*pi*frequency*tStatMiddle);

        tStatEnd = linspace(0,durStatEnd,durStatEnd*fs);
        statEnd  = square(2*pi*frequency*tStatEnd);

        %
        % Generate square waves for moving portions
        tMov1 = linspace(0,durMov1,durMov1*fs);
        mov1 = square(2*pi*frequency*tMov1);

        tMov2 = linspace(0,durMov2,durMov2*fs);
        mov2 = square(2*pi*frequency*tMov2);

        %% Concatenate and spatialize
        % intensityRampMain = [intensityRampStart intensityRampMov1 intensityRampMiddle intensityRampMov2 intensityRampEnd];
        allPortions = [statStart mov1 statMiddle mov2 statEnd];
        [stim, ~, ~, rActual, ~] = local_SOFAspat(allPortions',Obj,azi(stimNo),ele,rMain);
        
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
        outCsv(stimNo+1, :) = {filename, frequency, durStatStart, durMov1, durStatMiddle, ...
            durMov2, durStatEnd, totalDur, rMov1(1), ...
            rMov1(2),rMov2(1),rMov2(2), startSpaceMov1, stopSpaceMov1, startSpaceMov2,...
            stopSpaceMov2, trajectory, azimuth(stimNo), target(stimNo)};
        audiowrite(strcat('./', wavDir, '/', filename, '.wav'), stim, fs);

    end % for stimNo

elseif whichBlock == 3 % Azimuth block
    for stimNo = 1:stimuliPerTrajectory
        frequency = (randi(4,1)+4)*100; % 500-800 Hz with round 100 values

        %% Set azimuth:
        switch trajectory
            case 1
                x = 0;
                y = 30; % left
            case 2
                x = 0;
                y = -30; % right
        end % switch trajectory block 3

        aziStatStart = [x x];
        aziStatEnd = [y y];

        aziMov1 = [x y];

        aziMain = [aziStatStart,aziMov1,aziStatEnd];

        %% Set durations
        % Static portions: 600-900 ms, with round 100 ms values, in secs
        durStatStart = (randi(4,1)+5)/10;
        durStatEnd = (randi(4,1)+5)/10;

        % % Moving portion: duration calculated from velocity
        % % arcLength = r*phi where phi is radian
        % % WARNING - duration is too short in this way
        % arcLength = rBlock3(stimNo)*deg2rad(30);
        % durMov1 = arcLength/v;
        durMov1 = (randi(4,1)+5)/10;
        
        totalDur = durStatStart+durMov1+durStatEnd; % in s

        %% Generate square waves and intensity ramps for each portion

        % Generate square waves for stationary portions
        tStatStart = linspace(0,durStatStart,durStatStart*fs);
        statStart = square(2*pi*frequency*tStatStart);

        tStatEnd = linspace(0,durStatEnd,durStatEnd*fs);
        statEnd  = square(2*pi*frequency*tStatEnd);

        %
        % Generate square waves for moving portions
        tMov1 = linspace(0,durMov1,durMov1*fs);
        mov1 = square(2*pi*frequency*tMov1);

        %% Concatenate and spatialize
        allPortions = [statStart mov1 statEnd];
        rMain = linspace(rBlock3(stimNo),rBlock3(stimNo),length(allPortions));
        stim = local_SOFAspat(allPortions',Obj,aziMain,ele,rMain);

        % Save which space category this stim belongs to
        if ismember(rBlock3(stimNo),PPS)
            whichSpace = 1;
        elseif ismember(rBlock3(stimNo),ARS)
            whichSpace = 2;
        else
            whichSpace = 3;
        end
        
        % Save results to wav, add parameters to the cell array later saved out
        % to csv
        digits = ceil(log10(stimuliPerTrajectory + 1));
        stimulusdigits = ceil(log10(stimNo + 1));
        temp = char('');
        for digind = 1:(digits-stimulusdigits)
            temp = strcat(temp, '0');
        end

        filename = strcat(wavDir, '-', temp, num2str(stimNo));
        outCsv(stimNo+1, :) = {filename, frequency, durStatStart, durMov1, ...
            durStatEnd, totalDur, rBlock3(stimNo), whichSpace, trajectory};
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
    errorText = ['Signal (length: ', num2str(length(signal)), ') and r (length: ',...
        num2str(length(r)),') need to have the same length!'];
    error(errorText)
end
signal = db2mag(-100)*signal; % -26 dB to compensate for distance change from 1 to 0.05 and additional 65 dB to compensate for the particular type of HRTFs
[out, aziActual, eleActual, rActual, idx] = SOFAspat(signal./r(:),Obj,azi,ele,r);
end
