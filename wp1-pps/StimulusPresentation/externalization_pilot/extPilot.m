function extPilot(subNum,round,noBlocks)
% Usage: extPilot(subNum)
%
% Externalization pilot for CherISH DC6 WP1
%
% Input:
%   - subNum:   subject number, integer
%   - round:    first or second round, integer
%   - noBlocks: number of experimental blocks, integer
%
% Sounds are played in 2x2 blocks: 
%   - 2 variable f0 blocks (60 trials: 20 for each distance) 
%   - 2 constant f0 blocks (60 trials: 20 for each distance)
% Crossed with:
%   - 2 low intensity blocks
%   - 2 high intensity blocks
%
% 3AFC distance judgement task (in my head, at 20 cm, at 1 m)
%
% One trial consists of:
%   - Stimulus
%   - Prompt showing the three distances and the corresponding button
%   - Unspeeded response
%
% Responses recorded:
%   - Keypress (which distance was perceived)
%   - Accuracy (1 or 0)

%% Set order of source intensity blocks
sourceInt = [zeros(1,noBlocks/2), ones(1,noBlocks/2)]; % high or low source intensity
f0Condition = [zeros(1,noBlocks/2), ones(1,noBlocks/2)]; % constant or variable f0
sourceInt = sourceInt(randperm(length(sourceInt)));

%% Load/set params, stimuli, check for conflicts
% Open stimArray file
if mod(subNum,2) == 0 % even subjects start with constant f0 condition
    if round == 1 % open constant f0 stimArray corresponding to the subject
        stimArrayName = strcat('stimArrayExt',num2str(subNum),'.mat');
    else % open variable f0 stimArray in the 2nd round
        stimArrayName = 'stimArrayExt.mat';
    end
else % odd subjects start with variable f0 condition
    if round == 1 % open variable f0 stimArray in the 1st round
        stimArrayName = 'stimArrayExt.mat';
    else % open constant f0 condition corresponding to the subject
        stimArrayName = strcat('stimArrayExt',num2str(subNum),'.mat');
    end
end
stimArrayFileStruct = dir(stimArrayName);

% if there was no stimArray file or there were multiple
if isempty(stimArrayFileStruct) || length(stimArrayFileStruct)~=1
    error('Either found too many or no stimulus array files!');
else
    stimArrayFile = [stimArrayFileStruct.folder, '/', stimArrayFileStruct.name];
end

% user message
disp([char(10), 'Loading params and stimuli, checking ',...
    'for existing files for the subject']);

% A function handles all stimulus sorting to blocks and potential conflicts
% with earlier recordings for same subject
[stimArray, ~, ~,...
    startBlockNo, blockIdx, trialIdx,...
    logVar, subLogF, returnFlag,...
    logHeader, stimTypes] = handleParams(subNum, stimArrayFile, noBlocks);

% if there was early return from handleParams.m, abort
if returnFlag
    return
end

% Collect column numbers from the resulting stimArray
%%%%%%% HARD-CODED VALUES %%%%%%%
audioColumn = 8;
distanceColumn = 4;
stimTypeColumn = 9;
fsColumn = 7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Stimulus features for triggers + logging
fs = cell2mat(stimArray(1,fsColumn));
totalDur = cell2mat(logVar(2:end, strcmp(logHeader, 'dur'))); % sound duration

% user message
disp([char(10), 'Extracted stimulus features']);

%% Initialize Psychtoolbox (screen and audio)
PsychDefaultSetup(1);
Screen('Preference', 'SkipSyncTests', 1);
InitializePsychSound(1);

% Select audio device
device = []; 
% tmpDevices = PsychPortAudio('GetDevices');
% for i = 1:numel(tmpDevices)
%     % if strcmp(tmpDevices(i).DeviceName, 'Headphones (Conexant HD Audio headphone)')
%     % if strcmp(tmpDevices(i).DeviceName, 'Speakers/Headphones (Realtek(R) Audio)')
%         % device = tmpDevices(i).DeviceIndex;
%     end
% end
% mode is simple playback
mode = 1;
% reqlatencyclass is set to low-latency
reqLatencyClass = 2;
% 2 channels output
nrChannels = 2;

% user message
disp([char(10), 'Set audio parameters']);

% Define the specific keys we use
KbName('UnifyKeyNames');
keys = struct;
keys.abort = KbName('ESCAPE');
keys.go = KbName('SPACE');
keys.intern = KbName('1');
keys.near = KbName('2');
keys.far = KbName('3');
respKeys = [keys.intern,keys.near,keys.far];

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
screens=Screen('Screens');
screenNumber=max(screens);  % look into XOrgConfCreator and XOrgConfSelector
% screenNumber = 1;

% open stimulus window
[win, rect] = Screen('OpenWindow', screenNumber, backGroundColor);

% query frame duration for window
ifi = Screen('GetFlipInterval', win);
% set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% Setup the text type for the window
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

% open PsychPortAudio device for playback
pahandle = PsychPortAudio('Open', device, mode, reqLatencyClass, fs, nrChannels);
% get and display device status
pahandleStatus = PsychPortAudio('GetStatus', pahandle);
disp([char(10), 'PsychPortAudio device status: ']);
disp(pahandleStatus);

% initial start & stop of audio device to avoid potential initial latencies
tmpSound = zeros(2, fs/10);  % silence
tmpBuffer = PsychPortAudio('CreateBuffer', [], tmpSound);  % create buffer
PsychPortAudio('FillBuffer', pahandle, tmpBuffer);  % fill the buffer of audio device with silence
PsychPortAudio('Start', pahandle, 1);  % start immediately
PsychPortAudio('Stop', pahandle, 1);  % stop when playback is over

% user message
disp([char(10), 'Ready to start the experiment']);

%% Triggers
% basic triggers for trial start, sound onset and response
trig = struct;
trig.resp = 220;

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

% add triggers to logging / results variable
logVar(2:end, strcmp(logHeader, 'trigger')) = num2cell(trig.stimType);

% user message
disp([char(10), 'Set up triggers']);

% Set iti in sec
iti = 1;

% response time interval
respInt = 10; 

% response variables preallocation
accuracy = nan(size(stimArray, 1), 1);
response = accuracy;

% set flag for aborting experiment
abortFlag = 0;
% hide mouse
HideCursor(screenNumber);
% suppress keyboard input to command window
ListenChar(-1);
% realtime priority
Priority(1);

% user message
disp([char(10), 'Initialized psychtoolbox basics, opened window, ',...
    'started PsychPortAudio device']);

%% Instructions phase
% instructions text
instrText = ['You will hear sounds at different distances: at 1 m, at 20 cm, or in your head. \n',...
    'After each sound, indicate how far you think the sound was from you. \n',...
    'The experiment will be divided into blocks with some rest in between. \n',...
    'Some blocks will have soft sounds, while other blocks will be louder. \n\n', ...
    '<Press SPACE to continue>'];

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
    ListenChar(0);
    Priority(0);
    RestrictKeysForKbCheck([]);
    PsychPortAudio('Close');
    Screen('CloseAll');
    ShowCursor(screenNumber);
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
        SPL = db(0.2); 
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
        % Intensity roving:
        audioData = 10^(SPL/20)*audioData;
        buffer(end+1) = PsychPortAudio('CreateBuffer', [], audioData');
    end

    % counter for trials in given block
    trialCounterForBlock = 0;

    % user message
    disp([char(10), 'Buffered all stimuli for block ', num2str(block),...
        ', showing block start message']);

    % block starting text
    if sourceInt(block) == 0
        intensityInfo = 'This block will have relatively soft sounds.';
    else
        intensityInfo = 'This block will have relatively loud sounds.';
    end

    blockStartText = ['Beginning block ', num2str(block), '. \n\n',...
        intensityInfo, '\n',...
        'There will be ', num2str(length(trialList)), ' trials in the block.\n\n\n',...
        'Ready? Press SPACE to begin.'];

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
        ShowCursor(screenNumber);
        diary off
        return;
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

        % user message
        disp([char(10), 'Starting trial ', num2str(trialCounterForBlock)]);

        % fill audio buffer with next stimuli
        PsychPortAudio('FillBuffer', pahandle, buffer(trialCounterForBlock));

        % wait till we are 100 ms from the start of the playback
        while GetSecs-trialStart <= iti-100
            WaitSecs(0.001);
        end

        % blocking playback start for precision
        startTime = PsychPortAudio('Start', pahandle, 1, trialStart+iti, 1);

        % user message
        disp(['Audio started at ', num2str(startTime-trialStart), ' secs after trial start']);
        disp(['(Target ITI was ', num2str(iti), ' secs)']);

        % Wait for sound to end, then display 3AFC prompt
        WaitSecs(totalDur(trial));
        promptText = ['How far was the sound? \n\n',...
                      '1 - in my head \n 2 - close to my ear \n 3 - farther from my ear'];
        Screen('FillRect', win, backGroundColor);
        DrawFormattedText(win, promptText, 'center', 'center', textColor);
        Screen('Flip', win);

        % Wait for response (target detection)
        respFlag = 0;
            while GetSecs <= startTime+totalDur(trial)+respInt % response interval hasn't ended
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyIsDown % if subject responded
                    if ismember(find(keyCode),respKeys) % and it was a response button
                        respFlag = 1;
                        response(trial) = str2double(KbName(find(keyCode)));
                        if find(keyCode) == keys.intern % if response was in-the-head
                            if cell2mat(stimArray(trial,distanceColumn)) == 0 % and distance is 0
                                accuracy(trial) = 1;
                            else
                                accuracy(trial) = 0;
                            end
                        elseif find(keyCode) == keys.near % if response was 20 cm
                            if cell2mat(stimArray(trial,distanceColumn)) == 0.2 % and distance is 20 cm
                                accuracy(trial) = 1;
                            else
                                accuracy(trial) = 0;
                            end
                        elseif find(keyCode) == keys.far % if response is 1 m
                            if cell2mat(stimArray(trial,distanceColumn)) == 1 % and distance is 1 m
                                accuracy(trial) = 1;
                            else
                                accuracy(trial) = 0;
                            end
                        end
                        break;
                    elseif find(keyCode) == keys.abort % if it was the escape button
                        abortFlag = 1;
                        break;
                    end
                end
            end

        % if abort was requested, quit
        if abortFlag
            ListenChar(0);
            Priority(0);
            RestrictKeysForKbCheck([]);
            PsychPortAudio('Close');
            Screen('CloseAll');
            ShowCursor(screenNumber);
            diary off
            return;
        end

        % user messages
            if isnan(accuracy(trial))
                disp('Subject did not respond in time');
            else
                disp(['Distance: ', num2str(cell2mat(stimArray(trial,distanceColumn))),' m']);
                disp(['Response: ', num2str(response(trial))]);
                disp(['Accuracy: ', num2str(accuracy(trial))]);
            end

        % Cumulative accuracy and RT in block
        blockAcc = sum(accuracy(trial-trialCounterForBlock+1:trial), 'omitnan')/trialCounterForBlock*100;
        disp(['Overall accuracy in block so far is ', num2str(blockAcc), '%']);

        % accumulating results in logging / results variable
        logVar(trial+1,strcmp(logHeader, 'sourceInt')) = {sourceInt(block)};
        logVar(trial+1,strcmp(logHeader, 'accuracy')) = {accuracy(trial)};
        logVar(trial+1,strcmp(logHeader, 'response')) = {response(trial)};

        % save logging/results variable
        save(subLogF, 'logVar');

    end  % trial for loop

    % user messages
    disp([char(10), char(10), 'Block no. ', num2str(block), ' has ended,'...
        'showing block-ending text to participant']);
    disp([char(10), 'Overall accuracy in block was ', num2str(blockAcc),'%']);

    %% Feedback to subject at the end of block
    % if not last block
    if block ~= noBlocks
        % block ending text
        blockEndText = ['End of block ', num2str(block), '! \n\n\n',...
            'You had ',num2str(blockAcc),'% correct responses in this block. \n',...
            'Press SPACE to begin next block.'];
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
            ListenChar(0);
            Priority(0);
            RestrictKeysForKbCheck([]);
            PsychPortAudio('Close');
            Screen('CloseAll');
            ShowCursor(screenNumber);
            diary off
            return;
        end
    elseif block == noBlocks % if this was the last block
        % user message
        disp([char(10), 'The task has ended.']);
            blockEndText = ['End of last block. \n',...
                'You had ',num2str(blockAcc),'% correct responses in this block. \n',...
                'This is the end of the experiment. This window closes soon...'];
        % uniform background
        Screen('FillRect', win, backGroundColor);
        % draw block-starting text
        DrawFormattedText(win, blockEndText, 'center', 'center', textColor);
        Screen('Flip', win);
        WaitSecs(5);

    end  % if last block

end  % block for loop
diary off

%% Ending, cleaning up
ListenChar(0);
Priority(0);
RestrictKeysForKbCheck([]);
PsychPortAudio('Close');
Screen('CloseAll');
ShowCursor;  
return
end