function runExp(subNum, noBlocks)
% Experimental script for CherISH wp1
%
% % Inputs:
% subNum    - participant ID (integer)
% noBlocks  - number of experimental blocks (even number)

% % Triggers:
% 100       - Block start
% 101-110   - 100 + number of block
% 150-166   - Motion onset EEG triggers (see below)
% 170-176   - Motion onset GSR triggers (only no-target trials; see below)
% 199       - Target onset
% 200       - Trial start (before ITI)
% 210       - Playback start
% 220       - Response (button press)

% % % Stimulus type triggers at motion onset:
% EEG   GSR     Trajectory  Target  Offset-azi
% --------------------------------------------
% 151   171     loom        no      -90
% 152   172     loom        no       90
% 153           loom        yes     -90
% 154           loom        yes      90
% 155   175     rec         no      -90
% 156   176     rec         no       90
% 157           rec         yes     -90
% 158           rec         yes      90
% 159   179     PPS         no      -90
% 160   180     PPS         no       90
% 161           PPS         yes     -90
% 162           PPS         yes      90
% 163   183     EPS         no      -90
% 164   184     EPS         no       90
% 165           EPS         yes     -90
% 166           EPS         yes      90

% try
diary(strcat('diaries/s',num2str(subNum),'_',string(datetime("today"))));
%% Input check
if mod(noBlocks,2) ~=0
    warning(['Number of experimental blocks must be an even number, so I''m setting noBlocks to ',...
        num2str(noBlocks+1),'!']);
    noBlocks = noBlocks+1;
end

% Set volume 
system('\\KFS\Fileserver\Users\pkovacs\SetVol40.bat.lnk');

%% Flags
flag.headphoneCheck = 0;
flag.triggers = 1;
breakBlock = 5; % after which block to take a break

%% Set order of source intensity blocks
sourceInt = [zeros(1,noBlocks/2), ones(1,noBlocks/2)];
if mod(subNum,2) == 0
    sourceInt = flip(sourceInt);
end

%% Load/set params, stimuli, check for conflicts
% Look for stimArray file in the subject's folder
stimArrayName = 'stimArray.mat';
stimArrayFileStruct = dir(stimArrayName);

% if there was no stimArray file or there were multiple
if isempty(stimArrayFileStruct) || length(stimArrayFileStruct)~=1
    error(['Either found too many or no stimArrayFile at ./subject', num2str(subNum), '/stimArray*.mat !!!']);
else
    stimArrayFile = [stimArrayFileStruct.folder, '/', stimArrayFileStruct.name];
end

% user message
disp([char(10), 'Loading params and stimuli, checking ',...
    'for existing files for the subject']);

% a function handles all stimulus sorting to blocks and potential conflicts
% with earlier recordings for same subject
[stimArray, ~, ~,...
    startBlockNo, blockIdx, trialIdx,...
    logVar, subLogF, returnFlag,...
    logHeader, stimTypes] = handleParams(subNum, stimArrayFile, noBlocks);

% if there was early return from expParamsHandler.m, abort
if returnFlag
    return
end

% Collect column numbers from the resulting stimArray
%%%%%%% HARD-CODED VALUES %%%%%%%
targetColumn = 10;
audioColumn = 15;
stimTypeColumn = 16;
fsColumn = 14;
durColumn = 3;
durStatColumn = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = stimArray{1,fsColumn};

%% Stimulus features for triggers + logging
totalDurCue = cell2mat(logVar(2:end, strcmp(logHeader, 'totalDur'))); % duration of cue (w/o target or gap)
totalDurWithTarget = totalDurCue + 0.2; % in s
durStatOnset = cell2mat(logVar(2:end, strcmp(logHeader, 'durStatOnset'))); % required for EEG triggering
% durStatOffset = cell2mat(logVar(2:end, strcmp(logHeader, 'durStatOffset'))); % for EEG triggering

% user message
disp([char(10), 'Extracted stimulus features']);

%% Initialize PsychPortAudio and Screen Number
addpath('C:\Users\experimentator.KFS\Documents\cherish\wp1-pps\StimulusPresentation\externalization_pilot');

%%% Thanks to David Meijer for providing code

InitializePsychSound(1); % 1 pushes for low latency

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select sound card %%% OptionsDialog.m needed
%%%%%%%%%%%%%%%%%%%%%%%%%

computer_IDs = {};
sound_cards = {};
screen_numbers = {};

%On what computer are we?
S.computer_ID = getenv('computername');
computer_matches = cellfun(@(x) strcmp(x,S.computer_ID),computer_IDs);
if any(computer_matches)
    computer_nr = find(computer_matches);
else
    computer_nr = numel(computer_IDs)+1;
    computer_IDs{end+1,1} = S.computer_ID;
end


%Retrieve all available sound devices
sound_devices = PsychPortAudio('GetDevices');
for i=1:numel(sound_devices)
    full_name_audioDevices{i} = [num2str(sound_devices(i).DeviceIndex), '. ', sound_devices(i).HostAudioAPIName, ' - ', ...
        sound_devices(i).DeviceName, ' - NrOfOutputChan = ', num2str(sound_devices(i).NrOutputChannels)];
    full_name_audioDevices{i} = regexprep(full_name_audioDevices{i},'[\n\r]+','');  %Erase newlines and carriage returns
end

%Check if the saved audiocard for this computer still exists, if so set as default
sound_card_default = 0;
if (size(sound_cards,1) >= computer_nr) && ~isempty(sound_cards{computer_nr,1})
    if any(cellfun(@(x) strcmp(x,sound_cards{computer_nr,1}),full_name_audioDevices))
        sound_card_default = find(cellfun(@(x) strcmp(x,sound_cards{computer_nr,1}),full_name_audioDevices));
    end
end

%If no known sound card was found, default to a Windows WASAPI sound card if available
if ~sound_card_default
    WASAPI_devices = strfind(full_name_audioDevices, 'WASAPI');             %cell-array with indices of characters of where 'WASAPI' was found in full_name_audioDevices
    Loudspeaker_devices = strfind(full_name_audioDevices, 'Lautsprecher'); % keep only multi-channel (termed Lautsprecher)
    first_WASAPI_Lsp_device = find(and(cellfun(@(x) ~isempty(x),WASAPI_devices),cellfun(@(x) ~isempty(x),Loudspeaker_devices)),1); %first device with WASAPI and Lautsprecher in the name
    if ~isempty(first_WASAPI_Lsp_device)
        sound_card_default = first_WASAPI_Lsp_device;
    else
        sound_card_default = 1;                                             %Default to first sound card if WASAPI was not found
    end
end

%Let the user select the appropriate sound card
[selected_option,continueBool] = OptionsDialog(full_name_audioDevices,sound_card_default,'Sound Card','Select the sound card');
if ~continueBool
    error('Program is aborted');                        %dialog window was closed
end
S.sound_card = full_name_audioDevices{selected_option};
S.sound_device_idx = sound_devices(selected_option).DeviceIndex;

if strcmp(getenv('computername'),'LAB-BROWN')  % Manually overwrite the default sampling rate on the lab computer (default is 44100, but that leads to auditory delays relative to visual)
    S.sound_sample_rate = 48000;
else
    S.sound_sample_rate = sound_devices(selected_option).DefaultSampleRate;
end

%If this sound_card was not known yet, then remember it for this computer
if (size(sound_cards,1) < computer_nr) || isempty(sound_cards{computer_nr,1}) || ~strcmp(S.sound_card,sound_cards{computer_nr,1})
    sound_cards{computer_nr,1} = S.sound_card;
end

%%%%%%%%%%%%%%%%%%%%%
%%% Select screen %%% OptionsDialog.m needed
%%%%%%%%%%%%%%%%%%%%%

available_screens = Screen('Screens');

%Check if the saved screen_number is available to set as default choice
screen_number_default = NaN;
if (size(screen_numbers,1) >= computer_nr) && ~isempty(screen_numbers{computer_nr,1})
    if ismember(screen_numbers{computer_nr,1},available_screens)
        screen_number_default = find(screen_numbers{computer_nr,1} == available_screens);
    end
end
if isnan(screen_number_default)
    screen_number_default = numel(available_screens);
end

%Let the user select the appropriate screen number
[selected_option,continueBool] = OptionsDialog(num2cell(available_screens),screen_number_default,'Screen Number','Select the screen number');
if ~continueBool
    error('Program is aborted');                                            %dialog window was closed
end
S.screen_number = available_screens(selected_option);

%If this screen_number was not known yet, then remember it for this computer
if (size(screen_numbers,1) < computer_nr) || isempty(screen_numbers{computer_nr,1}) || ~strcmp(S.screen_number,screen_numbers{computer_nr,1})
    screen_numbers{computer_nr,1} = S.screen_number;
end

fs = 48000;
pahandle = PsychPortAudio('Open', S.sound_device_idx, 1, 1, fs, 4);
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1);

% user message
disp([char(10), 'Set audio parameters']);

% Define the specific keys we use
KbName('UnifyKeyNames');
keys = struct;
keys.abort = KbName('ESCAPE');
keys.go = KbName('SPACE');
keys.respTarget = KbName('RETURN');
% keys.left = KbName('LeftArrow');
% keys.right = KbName('RightArrow');
% keys.up = KbName('UpArrow');

% restrict keys to the ones we use
keysFields = fieldnames(keys);
keysVector = zeros(1, length(keysFields));
for f = 1:length(keysFields)
    keysVector(f) = keys.(keysFields{f});
end
RestrictKeysForKbCheck(keysVector);

% Force costly mex functions into memory to avoid latency later on
GetSecs; WaitSecs(0.1); KbCheck();

% screen params, screen selection
backGroundColor = [0 0 0];
textColor = [255 255 255];

% open stimulus window
[win, rect] = Screen('OpenWindow', S.screen_number, backGroundColor);

% query frame duration for window
ifi = Screen('GetFlipInterval', win);
% set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% Setup the text type for the window
Screen('TextFont', win, 'Ariel');
Screen('TextSize', win, 30);

% set up a central fixation cross into a texture / offscreen window
% get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(rect);
% Here we set the size of the arms of our fixation cross
fixCrossDimPix = 40;
% set the coordinates
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];
% set the line width for our fixation cross
lineWidthPix = 4;
% command to draw the fixation cross
fixCrossWin = Screen('OpenOffscreenWindow', win, backGroundColor, rect);
Screen('BlendFunction', fixCrossWin, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('DrawLines', fixCrossWin, allCoords,...
    lineWidthPix, textColor, [xCenter yCenter], 2);

% initial start & stop of audio device to avoid potential initial latencies
tmpSound = zeros(4, fs/10);  % silence
tmpBuffer = PsychPortAudio('CreateBuffer', [], tmpSound);  % create buffer
PsychPortAudio('FillBuffer', pahandle, tmpBuffer);  % fill the buffer of audio device with silence
PsychPortAudio('Start', pahandle, 1);  % start immediately
PsychPortAudio('Stop', pahandle, 1);  % stop when playback is over


% user message
disp([char(10), 'Ready to start the experiment']);

%% Triggers
% basic triggers for trial start, sound onset and response
trig = struct;
trig.trialStart = 200;
trig.playbackStart = 210;
trig.resp = 220;
trig.blockStart = 100;
trig.target = 199;

% triggers for stimulus types, based on the number of unique stimulus types
% we assume that stimTypes is a cell array (with headers) that contains the
% unique stimulus feature combinations, with an index for each combination
% in the last column
uniqueStimTypes = cell2mat(stimTypes(2:end,end));
if length(uniqueStimTypes) > 49
    error('Too many stimulus types for properly triggering them');
end
% triggers for stimulus types are integers in the range 151-199
trigTypes = uniqueStimTypes+150;
% add trigger info to stimTypes cell array as an extra column
stimTypes = [stimTypes, [{'trigger'}; num2cell(trigTypes)]];

% create triggers for stimulus types, for all trials
trig.stimType = cell2mat(stimArray(:, stimTypeColumn))+150;
trig.gsr  = cell2mat(stimArray(:, stimTypeColumn))+170;

% add triggers to logging / results variable
logVar(2:end, strcmp(logHeader, 'trigger')) = num2cell(trig.stimType);

% user message
disp([char(10), 'Set up triggers']);

% set random ITI between 2600-2800 ms, with round 100 ms values
iti = (randi([26,28], [size(stimArray, 1) 1]))/10;  % in secs

% response time interval
respInt = 1;

% response variables preallocation
respTime = nan(size(stimArray, 1), 1);
accuracy = respTime;

% set flag for aborting experiment
abortFlag = 0;
% hide mouse
HideCursor(S.screen_number);
% suppress keyboard input to command window
ListenChar(-1);
% realtime priority
Priority(1);

if flag.triggers
    % init serial port control
    % if ~exist(TB, 'var')
        TB = IOPort('OpenSerialPort', 'COM3'); % replace with number of COM port % test
    % end

    % Read data from the TriggerBox
    Available = IOPort('BytesAvailable', TB);
    if(Available > 0)
        disp(IOPort('Read', TB, 0, Available));
    end

    % Set the port to zero state 0
    IOPort('Write', TB, uint8(0), 0);
    pause(0.01);
end

% user message
disp([char(10), 'Initialized psychtoolbox basics, opened window, ',...
    'started PsychPortAudio device']);

%% Instructions phase
% instructions text
instrText = ['Sie werden jetzt Töne hören, die sich in verschiedenen Richtungen bewegen. \n',...
    'Nach einigen dieser Töne hören Sie einen kurzen Zielton. \n\n',...
    'Erkenne den Zielton so rasch und pünktlich wie möglich: Drücken Sie die Eingabetaste, wenn Sie ihn gehört haben. \n',...
    'Drücken Sie nichts, wenn Sie den Zielton nicht gehört haben. \n\n',...
    '<Drücken Sie die Leertaste, um fortzufahren>'];

% write instructions to text
Screen('FillRect', win, backGroundColor);
DrawFormattedText(win, instrText, 'center', 'center', textColor);
Screen('Flip', win);

% user message
disp([char(10), 'Showing the instructions text right now...']);

% wait for key press to start
while 1
    [keyIsDown, ~, keyCode] = KbCheck;
    if keyIsDown
        % if subject is ready to start
        if find(keyCode) == keys.go
            break;
        elseif find(keyCode) == keys.abort
            % if abort was requested
            abortFlag = 1;
            break;
        end
    end
end

if abortFlag
    if flag.triggers
        % ppdev_mex('Close', 1);
        IOport('Close', TB);
    end
    ListenChar(0);
    Priority(0);
    RestrictKeysForKbCheck([]);
    PsychPortAudio('Close');
    Screen('CloseAll');
    ShowCursor(S.screen_number);
    diary off
    return;
end

% user message
disp([char(10), 'Subject signalled she/he is ready, we go ahead with the task']);

%% Blocks loop

% start from the block specified by the parameters/settings parts
for block = startBlockNo:noBlocks

    % If low source intensity block, scale down
    if sourceInt(block) == 0
        SPL = db(0.2/2);
        disp('Low intensity block');
    else
        SPL = db(1);
        disp('High intensity block');
    end

    % uniform background
    Screen('FillRect', win, backGroundColor);
    Screen('Flip', win);

    % wait for releasing keys before going on
    releaseStart = GetSecs;
    KbReleaseWait([], releaseStart+2);

    % fill a dynamic buffer with data for whole block
    % get trial index list for current block
    trialList = trialIdx(blockIdx==block);
    buffer = [];
    for trial = min(trialList):max(trialList)
        audioData = stimArray{trial, audioColumn};
        audioData = 10^(SPL/20)*audioData; % Control intensity

        % Set up audio triggers for StimTrak
        fsd=1000; % sampling rate of the trigger channel
        % SW=[ones(40,1); -ones(40,1)]; % square wave
        stimVec = [ones(40,1); -ones(40,1)];
        siglen = length(audioData); % signal length in samples
        % statlen = fs*stimArray{trial,durStatColumn}; % length of stationary onset
        % targlen = fs*0.1; % 100 ms target

        % SW=resample(SW, fs, fsd); % upsample to the audio sampling rate
        stimVec = resample(stimVec, fs, fsd);
        stimVec=[stimVec;zeros(siglen-length(stimVec),1)]; 
        % if stimArray{trial,targetColumn} == 1
        %     stimVec = [SW;zeros(statlen-length(SW),1);SW;zeros(siglen-statlen-length(SW)-2*targlen,1)];
        %     stimVec = [stimVec;zeros(targlen,1);SW;zeros(targlen-length(SW),1)]; % target trigger, if there is one
        %      % onset trigger, silence, motion trigger, silence
        % else
        %     stimVec = [SW;zeros(statlen-length(SW),1);SW;zeros(siglen-statlen-length(SW),1)]; % onset trigger, silence, motion trigger, silence
        % 
        % end
        % stimVec=[stimVec zeros(length(stimVec),1)]; % make it two channels
        stimVec = [stimVec zeros(siglen,1)];
        % stimVec = [zeros(siglen,1) stimVec];

        % Fill the buffer with the audio and the StimTrak trigger
        buffer(end+1) = PsychPortAudio('CreateBuffer', [], [audioData;  stimVec']);
    end

    % counter for trials in given block
    trialCounterForBlock = 0;

    % user message
    disp([char(10), 'Buffered all stimuli for block ', num2str(block),...
        ', showing block start message']);

    % block starting text
    if sourceInt(block) == 0
        intensityInfo = 'In diesem Block sind alle Töne relativ leise.';
    else
        intensityInfo = 'In diesem Block sind alle Töne relativ laut.';
    end

    blockStartText = ['Block ', num2str(block), ' beginnt. \n\n',...
        intensityInfo, '\n',...
        'Es gibt ', num2str(length(trialList)), ' Proben im Block.\n',...
        'Vergessen Sie nicht ENTER zu drücken, wenn Sie den kurzen Zielton hören. \n\n\n',...
        'Bereit? Drücken Sie die Leertaste um zu beginnen.'];

    % uniform background
    Screen('FillRect', win, backGroundColor);
    % draw block-starting text
    DrawFormattedText(win, blockStartText, 'center', 'center', textColor);
    Screen('Flip', win);
    % wait for key press to start
    while 1
        [keyIsDown, ~, keyCode] = KbCheck;

        if keyIsDown
            % if subject is ready to start
            if find(keyCode) == keys.go
                break;
            elseif find(keyCode) == keys.abort
                abortFlag = 1;
                break;
            end
        end
    end

    if abortFlag
        if flag.triggers
            % ppdev_mex('Close', 1);
            IOport('Close',TB);
        end
        ListenChar(0);
        Priority(0);
        RestrictKeysForKbCheck([]);
        PsychPortAudio('Close');
        Screen('CloseAll');
        ShowCursor(S.screen_number);
        diary off
        return;
    end

    if flag.triggers
        % block start trigger + block number trigger
        IOPort('Write',TB,uint8(trig.blockStart),0);
        pause(0.01);
        IOPort('Write', TB, uint8(0), 0);
        pause(0.01);

        IOPort('Write', TB, uint8(trig.blockStart+block),0);
        pause(0.01);
        IOPort('Write', TB, uint8(0), 0);
        pause(0.01);
    end


    %% Trials loop

    % trial loop (over the trials for given block)
    for trial = min(trialList):max(trialList)

        % relative trial number (trial in given block)
        trialCounterForBlock = trialCounterForBlock+1;

        % background with fixation cross, get trial start timestamp
        Screen('CopyWindow', fixCrossWin, win);
        Screen('DrawingFinished', win);
        trialStart = Screen('Flip', win);

        if flag.triggers
            % trial start trigger
            IOPort('Write',TB,uint8(trig.trialStart),0);
            pause(0.01);
            IOPort('Write', TB, uint8(0), 0);
            pause(0.01);
        end

        % user message
        disp([char(10), 'Starting trial ', num2str(trialCounterForBlock)]);

        % fill audio buffer with next stimuli
        PsychPortAudio('FillBuffer', pahandle, buffer(trialCounterForBlock));
        % PsychPortAudio('FillBuffer', pahandle, [sigpair stimVec]');

        % wait till we are 100 ms from the start of the playback
        while GetSecs-trialStart <= iti(trial)-100
            WaitSecs(0.001);
        end

        % blocking playback start for precision
        startTime = PsychPortAudio('Start', pahandle, 1, trialStart+iti(trial), 1);

        if flag.triggers
            % playback start trigger
            IOPort('Write', TB, uint8(trig.playbackStart),0);
            pause(0.01);
            IOPort('Write', TB, uint8(0),0);
            pause(0.01);

            % wait for motion onset
            delay = (startTime + durStatOnset(trial)) - GetSecs;
            if delay > 0
                WaitSecs(delay);
            end

            % stimulus type trigger at motion onset (EEG)
            IOPort('Write',TB,uint8(trig.stimType(trial)),0);
            pause(0.01);
            IOPort('Write', TB, uint8(0), 0);
            pause(0.01);

            % GSR stimulus type trigger at motion onset (only no-target trials)
            if cell2mat(stimArray(trial,targetColumn)) == 0 % if it's a no-target trial
                IOPort('Write',TB,uint8(trig.gsr(trial)),0);
                pause(0.01);
                IOPort('Write', TB, uint8(0), 0);
                pause(0.01);
            end

            % wait for target if there is one
            if cell2mat(stimArray(trial,targetColumn)) == 1
                targ_delay = (startTime + totalDurWithTarget(trial)-durStatOnset(trial) - 0.01) - GetSecs;
                if targ_delay > 0
                    WaitSecs(targ_delay);
                end

                % send target trigger
                IOPort('Write',TB,uint8(trig.target),0);
                pause(0.01);
                IOPort('Write', TB, uint8(0), 0);
                pause(0.01);
            end
        end

        % user message
        disp(['Audio started at ', num2str(startTime-trialStart), ' secs after trial start']);
        disp(['(Target ITI was ', num2str(iti(trial)), ' secs)']);

        % Wait for response (target detection)
        respFlag = 0;

        % while GetSecs >= startTime+totalDurCue(trial) % cue has ended
        while GetSecs <= startTime+totalDurWithTarget(trial)+respInt % response interval hasn't ended
            [keyIsDown, respSecs, keyCode] = KbCheck;
            if keyIsDown % if subject responded
                if find(keyCode) == keys.respTarget % and it was the detection button
                    respFlag = 1;
                    if cell2mat(stimArray(trial,targetColumn)) == 1 % and there is a target in the trial
                        accuracy(trial) = 1;
                    else
                        accuracy(trial) = 0;
                    end
                    % response trigger
                    if flag.triggers
                        IOPort('Write',TB,uint8(trig.resp),0);
                        pause(0.01);
                        IOPort('Write', TB, uint8(0),0);
                        pause(0.01);
                    end
                    break;
                elseif find(keyCode) == keys.abort % if it was the escape button
                    abortFlag = 1;
                    break;
                end
            else % if subject did not respond
                if cell2mat(stimArray(trial,targetColumn)) == 0 % and there is no target in the trial
                    accuracy(trial) = 1;
                else
                    accuracy(trial) = 0;
                end
            end
        end
        % end

        % if abort was requested, quit
        if abortFlag
            if flag.triggers
                IOPort('Close', TB);
            end
            ListenChar(0);
            Priority(0);
            RestrictKeysForKbCheck([]);
            PsychPortAudio('Close');
            Screen('CloseAll');
            ShowCursor(S.screen_number);
            diary off
            return;
        end

        % switch visual right when the audio finishes
        % Screen('Flip', win, startTime+totalDurWithTarget(trial)-0.5*ifi);

        % response time into results variable
        if respFlag
            respTime(trial) = 1000*(respSecs-(startTime+totalDurWithTarget(trial)));
        end

        % user messages
        if cell2mat(stimArray(trial,targetColumn)) == 1 % if it's a target trial
            if isnan(respTime(trial))
                disp('Subject did not respond in time');
            else
                disp(['RT: ', num2str(respTime(trial))]);
            end
        end

        % Cumulative accuracy and RT in block
        blockAcc = sum(accuracy(trial-trialCounterForBlock+1:trial), 'omitnan')/trialCounterForBlock*100;
        disp(['Overall accuracy in block so far is ', num2str(blockAcc), '%']);
        blockRT = mean(respTime(trial-trialCounterForBlock+1:trial),'omitnan');
        disp(['Average RT in block so far is ', num2str(blockRT), ' ms']);

        % accumulating results in logging / results variable
        logVar(trial+1,strcmp(logHeader, 'sourceInt')) = {sourceInt(block)};
        logVar(trial+1,strcmp(logHeader, 'accuracy')) = {accuracy(trial)};
        logVar(trial+1,strcmp(logHeader, 'respTime')) = {respTime(trial)};
        logVar(trial+1,strcmp(logHeader, 'promptness')) = {1/respTime(trial)};

        % save logging/results variable
        save(subLogF, 'logVar');

    end  % trial for loop

    % user messages
    disp([char(10), char(10), 'Block no. ', num2str(block), ' has ended,'...
        'showing block-ending text to participant']);
    disp([char(10), 'Overall accuracy in block was ', num2str(blockAcc),'%']);

    %% Feedback to subject at the end of block
    % if not last block and not a break
    if (block ~= noBlocks) && (block ~= breakBlock)
        % block ending text
        blockEndText = ['Ende von Block ', num2str(block), '! \n\n\n',...
            'Sie hatten ',num2str(blockAcc),'% richtige Antworten in diesem Block. \n',...
            'Ihre durchschnittliche Reaktionszeit in diesem Block war ',num2str(round(blockRT, 2)), ' ms. \n\n\n',...
            'Drücken Sie die Leertaste, um fortzufahren.'];
        % uniform background
        Screen('FillRect', win, backGroundColor);
        % draw block-starting text
        DrawFormattedText(win, blockEndText, 'center', 'center', textColor);
        Screen('Flip', win);
        % wait for key press to start
        while 1
            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown
                % if subject is ready to start
                if find(keyCode) == keys.go
                    break;
                elseif find(keyCode) == keys.abort
                    abortFlag = 1;
                    break;
                end
            end
        end
        if abortFlag
            if flag.triggers
                % ppdev_mex('Close', 1);
                IOPort('Close',TB);
            end
            ListenChar(0);
            Priority(0);
            RestrictKeysForKbCheck([]);
            PsychPortAudio('Close');
            Screen('CloseAll');
            ShowCursor(S.screen_number);
            diary off
            return;
        end

    % if not last block and there is a break
    elseif (block ~= noBlocks) && (block == breakBlock)
        % user message
        disp([char(10), 'There is a BREAK now!']);
        disp('Only the experimenter can start the next block - press "SPACE" when ready');
        
        % block ending text
        blockEndText = ['Ende von Block ', num2str(block), '. \n\n\n',... 
                'In diesem Block hatten Sie ', num2str(round(blockAcc, 2)), '% richtige Antworten.\n\n\n',... 
                'Jetzt machen wir eine kurze Pause.'];
        % uniform background
        Screen('FillRect', win, backGroundColor);
        % draw block-starting text
        DrawFormattedText(win, blockEndText, 'center', 'center', textColor);   
        Screen('Flip', win);

        while 1
            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown
                % if subject is ready to start
                if find(keyCode) == keys.go
                    break;
                elseif find(keyCode) == keys.abort
                    abortFlag = 1;
                    break;
                end
            end
        end
        if abortFlag
            if flag.triggers
                % ppdev_mex('Close', 1);
                IOPort('Close',TB);
            end
            ListenChar(0);
            Priority(0);
            RestrictKeysForKbCheck([]);
            PsychPortAudio('Close');
            Screen('CloseAll');
            ShowCursor(S.screen_number);
            diary off
            return;
        end

    % if this was the last block
    elseif block == noBlocks 
        % user message
        disp([char(10), 'The task has ended.']);
        blockEndText = ['Ende vom letzten Block. \n',...
            'Sie hatten ',num2str(blockAcc),'% richtige Antworten in diesem Block. \n',...
            'Ihre durchschnittliche Reaktionszeit im letzten Block war ',num2str(round(blockRT, 2)), ' ms. \n\n', ...
            'Das Experiment ist abgeschlossen. Dieses Fenster schließt sich bald...'];
        % uniform background
        Screen('FillRect', win, backGroundColor);
        % draw block-starting text
        DrawFormattedText(win, blockEndText, 'center', 'center', textColor);
        Screen('Flip', win);
        WaitSecs(5);

    end  % if last block

end  % block for loop
diary off
% catch
%     if flag.triggers
%     % ppdev_mex('Close', 1);
%     IOPort('Close',TB);
%     end
%     ListenChar(0);
%     Priority(0);
%     RestrictKeysForKbCheck([]);
%     PsychPortAudio('Close');
%     Screen('CloseAll');
%     ShowCursor;

%% Ending, cleaning up
if flag.triggers
    % ppdev_mex('Close', 1);
    IOPort('Close',TB);
end
ListenChar(0);
Priority(0);
RestrictKeysForKbCheck([]);
PsychPortAudio('Close');
Screen('CloseAll');
ShowCursor;
% Set volume back
system('C:\Users\experimentator.KFS\Desktop\SetVol100.bat.lnk');
return
end