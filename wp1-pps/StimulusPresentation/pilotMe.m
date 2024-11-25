function pilotMe(subNum, bigBlock, noBlocks, triggers)
% Experimental script for CherISH wp-1 pilot
% Adapted (-ing) from SFGmain.m

%% Input checks

switch bigBlock
    % Look for stimArray file in the subject's folder corresponding to the
    % block
    case 1
        stimArrayFileStruct = dir('stimArrayBlock1.mat');
    case 2
        stimArrayFileStruct = dir('stimArrayBlock2.mat');
    case 3
        stimArrayFileStruct = dir('stimArrayBlock3.mat');
end

% if there was no stimArray file or there were multiple
if isempty(stimArrayFileStruct) || length(stimArrayFileStruct)~=1
    error(['Either found too many or no stimArrayFile at ./subject', num2str(subNum), '/stimArray*.mat !!!']);
else
    stimArrayFile = [stimArrayFileStruct.folder, '/', stimArrayFileStruct.name];
end

clc;

%% Load/set params, stimuli, check for conflicts

% user message
disp([char(10), 'Loading params and stimuli, checking ',...
    'for existing files for the subject']);

% a function handles all stimulus sorting to blocks and potential conflicts
% with earlier recordings for same subject
[stimArray, ~, ~,...
    startBlockNo, blockIdx, trialIdx,...
    logVar, subLogF, returnFlag,...
    logHeader, stimTypes] = handleParams(subNum, stimArrayFile, noBlocks, bigBlock);

% if there was early return from expParamsHandler.m, abort
if returnFlag
    return
end

% user message
disp([char(10), 'Ready to start the experiment']);



%% Stimulus features for triggers + logging
endSpace = cell2mat(logVar(2:end, strcmp(logHeader, 'space')));
trajectory = cell2mat(logVar(2:end, strcmp(logHeader, 'trajectory')));
stimLength = cell2mat(logVar(2:end, strcmp(logHeader, 'totalDur')));
lastPortionLength = cell2mat(stimArray(:,7)); % required for defining response interval

% user message
disp([char(10), 'Extracted stimulus features']);

%% Triggers

% basic triggers for trial start, sound onset and response
trig = struct;
trig.trialStart = 200;
trig.playbackStart = 210;
trig.resp = 220;
trig.blockStart = 100;

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
if bigBlock < 3
    trig.stimType = cell2mat(stimArray(:, 19))+150;
else
    trig.stimType = cell2mat(stimArray(:, 13))+150;
end

% add triggers to logging / results variable
logVar(2:end, strcmp(logHeader, 'trigger')) = num2cell(trig.stimType);

% user message
disp([char(10), 'Set up triggers']);


%% Psychtoolbox initialization

% General init (AssertOpenGL, 'UnifyKeyNames')
PsychDefaultSetup(1);

% init PsychPortAudio with pushing for lowest possible latency
InitializePsychSound(1);
%% Audio parameters for PsychPortAudio

% sampling rate is derived from the HRTF set, 44100 Hz in this case
fs = 44100;

% Select audio device: this is where I could use o_ptb
% device = [];  % system default is our default as well
%
tmpDevices = PsychPortAudio('GetDevices');
for i = 1:numel(tmpDevices)
    if strcmp(tmpDevices(i).DeviceName, 'Lautsprecher (RME Fireface UCX ')
        device = tmpDevices(i).DeviceIndex;
    end
end

% mode is simple playback
mode = 1;
% reqlatencyclass is set to low-latency
reqLatencyClass = 2;
% 2 channels output
nrChannels = 2;

% user message
disp([char(10), 'Set audio parameters']);

% Define the specific keys we use
keys = struct;
keys.abort = KbName('ESCAPE');
keys.go = KbName('SPACE');
keys.respStop = KbName('RETURN');
keys.respDistancePPS = KbName('1');
keys.respDistanceARS = KbName('2');
keys.respDistanceEPS = KbName('3');
keys.respLeft = KbName('LeftArrow');
keys.respRight = KbName('RightArrow');
distanceResps = [keys.respDistancePPS keys.respDistanceARS keys.respDistanceEPS];

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

% open stimulus window
[win, rect] = Screen('OpenWindow', screenNumber, backGroundColor);

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

% set random ITI between 500-800 ms, with round 100 ms values
iti = (randi(4, [size(stimArray, 1) 1])+4)/10;  % in secs

% response time interval
respInt1 = 5; 
respInt2 = 20; 

% response variables preallocation
respTime = nan(size(stimArray, 1), 1);
accDistance = respTime;
accDirection = respTime;
respDistance = respTime;
respDirection = respTime;

% set flag for aborting experiment
abortFlag = 0;
% hide mouse
HideCursor(screenNumber);
% suppress keyboard input to command window
ListenChar(-1);
% realtime priority
Priority(1);

if triggers
    % init serial port control
    TB = IOPort('OpenSerialPort', 'COM3'); % replace with number of COM port % test

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
if bigBlock < 3
    instrText = ['Now you will hear pairs of sounds, each moving either closer to you \n',...
        'or farther. In the beginning, middle and end, the sound stops for a few moments. \n\n',...
        'Detect when the sound stopped for the LAST time: press ENTER at the end \n',...
        'of the second movement. \n\n',...
        'Then, you will have more time to answer how far the sound stopped from you: \n',...
        '1: uncomfortably close to my face \n',...
        '2: not extremely close, but still within arm''s reach \n',...
        '3: quite far, out of reach \n \n',...
        'Press SPACE to continue.'];
else
    instrText = ['Now you will hear sounds moving either to the left or to the right. \n',...
        'Detect as fast and accurately as possible which way the sound moved: \n',...
        'press LEFT arrow or RIGHT arrow. \n\n',...
        'Then, you will have more time to answer how far the sound was from you: \n',...
        '1: uncomfortably close to my face \n',...
        '2: not extremely close, but still within arm''s reach \n',...
        '3: quite far, out of reach \n \n',...
        'Press SPACE to continue.'];
end

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
    if triggers
        % ppdev_mex('Close', 1);
        IOport('Close', TB);
    end
    ListenChar(0);
    Priority(0);
    RestrictKeysForKbCheck([]);
    PsychPortAudio('Close');
    Screen('CloseAll');
    ShowCursor(screenNumber);
    return;
end

% user message
disp([char(10), 'Subject signalled she/he is ready, we go ahead with the task']);


%% Blocks loop

SPL = -30; % changes in the next block; values: SPL + [0 3 6]

% start from the block specified by the parameters/settings parts
for block = startBlockNo:noBlocks

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
        if bigBlock < 3
            audioData = stimArray{trial, 18};
            % Intensity roving:
            audioData = 10^(SPL/20)*audioData;
            buffer(end+1) = PsychPortAudio('CreateBuffer', [], audioData');
        else
            audioData = stimArray{trial, 12};
            % Intensity roving:
            audioData = 10^(SPL/20)*audioData;
            buffer(end+1) = PsychPortAudio('CreateBuffer', [], audioData');
        end
    end

    % counter for trials in given block
    trialCounterForBlock = 0;

    % user message
    disp([char(10), 'Buffered all stimuli for block ', num2str(block),...
        ', showing block start message']);

    % block starting text
    blockStartText = ['Beginning block ', num2str(block), '. \n\n',...
        'As always in this round, there will be ', num2str(length(trialList)), ' trials in the block.\n\n\n',...
        'Press SPACE to begin!'];

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
        if triggers
            % ppdev_mex('Close', 1);
            IOport('Close',TB);
        end
        ListenChar(0);
        Priority(0);
        RestrictKeysForKbCheck([]);
        PsychPortAudio('Close');
        Screen('CloseAll');
        ShowCursor(screenNumber);
        return;
    end

    if triggers
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

        if triggers
            % trial start trigger + trial type trigger
            IOPort('Write',TB,uint8(trig.trialStart),0);
            pause(0.01);
            IOPort('Write', TB, uint8(0), 0);
            pause(0.01);

            IOPort('Write',TB,uint8(trig.stimType),0);
            pause(0.01);
            IOPort('Write', TB, uint8(0), 0);
            pause(0.01);
        end

        % user message
        disp([char(10), 'Starting trial ', num2str(trialCounterForBlock)]);

        % fill audio buffer with next stimuli
        PsychPortAudio('FillBuffer', pahandle, buffer(trialCounterForBlock));

        % wait till we are 100 ms from the start of the playback
        while GetSecs-trialStart <= iti(trial)-100
            WaitSecs(0.001);
        end

        % blocking playback start for precision
        startTime = PsychPortAudio('Start', pahandle, 1, trialStart+iti(trial), 1);

        % playback start trigger
        if triggers
            IOPort('Write', TB, uint8(trig.playbackStart),0);
            pause(0.01);
            IOPort('Write', TB, uint8(0),0);
            pause(0.01);
        end

        % user message
        disp(['Audio started at ', num2str(startTime-trialStart), ' secs after trial start']);
        disp(['(Target ITI was ', num2str(iti(trial)), ' secs)']);

        % % prepare screen change for response period
        % Screen('CopyWindow', qMarkWin, win);
        % Screen('DrawingFinished', win);
        % Wait for response 1 (when did the movement stop?)
        respFlag1 = 0;
        
        while GetSecs <= startTime+stimLength(trial)+respInt1
            [keyIsDown, respSecs, keyCode] = KbCheck;
            if keyIsDown
                % if subject responded to movement stop (blocks 1-2)
                if bigBlock < 3
                    if find(keyCode) == keys.respStop
                        respFlag1 = 1;
                        % response trigger
                        if triggers
                            IOPort('Write',TB,uint8(trig.resp),0);
                            pause(0.01);
                            IOPort('Write', TB, uint8(0),0);
                            pause(0.01);
                        end
                        break;
                    elseif find(keyCode) == keys.abort
                        abortFlag = 1;
                        break;
                    end
                elseif bigBlock == 3 
                    if trajectory(trial) == 1 % left, 30 deg
                        if find(keyCode) == keys.respLeft
                            respFlag1 = 1;
                            accDirection(trial) = 1;
                            respDirection(trial) = 1;
                        elseif find(keyCode) == keys.respRight
                            respFlag1 = 1;
                            accDirection(trial) = 0;
                            respDirection(trial) = 2;
                        elseif find(keyCode) == keys.abort
                            abortFlag = 1;
                            break;
                        end
                    elseif trajectory(trial) == 2 % right, -30deg
                        if find(keyCode) == keys.respLeft
                            respFlag1 = 1;
                            accDirection(trial) = 0;
                            respDirection(trial) = 1;
                        elseif find(keyCode) == keys.respRight
                            respFlag1 = 1;
                            accDirection(trial) = 1;
                            respDirection(trial) = 2;
                        elseif find(keyCode) == keys.abort
                            abortFlag = 1;
                            break;
                        end
                    end               
                end
            end
        end

        % if abort was requested, quit
        if abortFlag
            if triggers
                IOPort('Close', TB);
            end
            ListenChar(0);
            Priority(0);
            RestrictKeysForKbCheck([]);
            PsychPortAudio('Close');
            Screen('CloseAll');
            ShowCursor(screenNumber);
            return;
        end

        % switch visual right when the audio finishes
        respStart = Screen('Flip', win, startTime+stimLength(trial)-0.5*ifi);
        
        % response time into results variable
        if respFlag1
            respTime(trial) = 1000*(respSecs-respStart-lastPortionLength(trial)); 
        end

        % user messages
        if isnan(respTime(trial))
            disp('Subject did not respond in time');
        end

        distanceText = ['Distance? \n\n 1: uncomfortably close to my face \n',...
    '2: not extremely close, but still within arm''s reach \n',...
    '3: quite far, out of reach \n \n'];

        % write instructions to text
        Screen('FillRect', win, backGroundColor);
        DrawFormattedText(win, distanceText, 'center', 'center', textColor);
        Screen('Flip', win);

        % user message
        disp(['Visual flip for response period start was ', num2str(respStart-startTime),...
            ' secs after audio start (should equal ', num2str(stimLength(trial)), ')']);

        % Wait for response 2 (how far was the sound at the end point?)
        respFlag2 = 0;
        % while GetSecs <= startTime+stimLength(trial)+respInt2
        while GetSecs-(startTime+stimLength(trial)) <= respInt2
            [keyIsDown, respSecs, keyCode] = KbCheck;
            if keyIsDown
                % if subject responded to distance
                if ismember(find(keyCode), distanceResps)
                    respFlag2 = 1;
                    % response trigger
                    if triggers
                        IOPort('Write',TB,uint8(trig.resp),0);
                        pause(0.01);
                        IOPort('Write',TB,uint8(0),0);
                        pause(0.01);
                    end
                    break;
                elseif find(keyCode) == keys.abort
                    abortFlag = 1;
                    break;
                end
            end
        end

        % if abort was requested, quit
        if abortFlag
            if triggers
                % ppdev_mex('Close', 1);
                IOPort('Close',TB);
            end
            ListenChar(0);
            Priority(0);
            RestrictKeysForKbCheck([]);
            PsychPortAudio('Close');
            Screen('CloseAll');
            ShowCursor(screenNumber);
            return;
        end

        % accuracy and logging of distance response
        if respFlag2
            if endSpace(trial) == 1
                if find(keyCode) == keys.respDistancePPS
                    disp('Subject''s response was accurate (PPS)');
                    accDistance(trial) = 1;
                    respDistance(trial) = 1;
                elseif find(keyCode) == keys.respDistanceARS
                    disp('Subject made an error (ARS instead of PPS)');
                    accDistance(trial) = 0;
                    respDistance(trial) = 2;
                elseif find(keyCode) == keys.respDistanceEPS
                    disp('Subject made an error (EPS instead of PPS)');
                    accDistance(trial) = 0;
                    respDistance(trial) = 3;
                end
            elseif endSpace(trial) == 2
                if find(keyCode) == keys.respDistanceARS
                    disp('Subject''s response was accurate (ARS)');
                    accDistance(trial) = 1;
                    respDistance(trial) = 2;
                elseif find(keyCode) == keys.respDistancePPS
                    disp('Subject made an error (PPS instead of ARS)');
                    accDistance(trial) = 0;
                    respDistance(trial) = 1;
                elseif find(keyCode) == keys.respDistanceEPS
                    disp('Subject made an error (EPS instead of ARS)');
                    accDistance(trial) = 0;
                    respDistance(trial) = 3;
                end
            elseif endSpace(trial) == 3
                if find(keyCode) == keys.respDistanceEPS
                    disp('Subject''s response was accurate (EPS)');
                    accDistance(trial) = 1;
                    respDistance(trial) = 3;
                elseif find(keyCode) == keys.respDistancePPS
                    disp('Subject made an error (PPS instead of EPS)');
                    accDistance(trial) = 0;
                    respDistance(trial) = 1;
                elseif find(keyCode) == keys.respDistanceARS
                    disp('Subject made an error (ARS instead of EPS)');
                    accDistance(trial) = 0;
                    respDistance(trial) = 2;
                end
            end
        end
       
        % cumulative accuracy in block
        blockAcc = sum(accDistance(trial-trialCounterForBlock+1:trial), 'omitnan')/trialCounterForBlock*100;
        disp(['Overall accuracy in block so far is ', num2str(blockAcc), '%']);

        % accumulating all results in logging / results variable
        logVar(trial+1,strcmp(logHeader, 'SPL')) = {SPL};
        logVar(trial+1,strcmp(logHeader, 'accDistance')) = {accDistance(trial)};
        logVar(trial+1,strcmp(logHeader, 'accDirection')) = {accDirection(trial)};
        logVar(trial+1,strcmp(logHeader, 'respTime')) = {respTime(trial)};
        logVar(trial+1,strcmp(logHeader, 'respSpace')) = {respDistance(trial)};
        logVar(trial+1,strcmp(logHeader, 'respDirection')) = {respDirection(trial)};
        
        % logVar(trial+1, 10:end-1) = {acc(trial), ...
        %     respTime(trial), iti(trial),...
        %     trialStart, startTime-trialStart,...
        %     respStart-startTime};

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
            'In this block, you had ', num2str(round(blockAcc, 2)), '% correct responses regarding distance.\n\n\n',...
            'Press SPACE to begin next block!'];
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
            if triggers
                % ppdev_mex('Close', 1);
                IOPort('Close',TB);
            end
            ListenChar(0);
            Priority(0);
            RestrictKeysForKbCheck([]);
            PsychPortAudio('Close');
            Screen('CloseAll');
            ShowCursor(screenNumber);
            return;
        end

        % if not last block and there is a break
        % elseif (block ~= noBlocks) && ismembertol(block, breakBlocks)
        %
        %     % user message
        %     disp([char(10), 'There is a BREAK now!']);
        %     disp('Only the experimenter can start the next block - press "SPACE" when ready');
        %
        %     % block ending text
        %     blockEndText = ['Vége a(z) ', num2str(block), '. blokknak!\n\n\n',...
        %             'Ebben a blokkban a próbák ', num2str(round(blockAcc, 2)), '%-ra adott helyes választ.\n\n\n',...
        %             'Most tartunk egy rövid szünetet, a kísérletvezető hamarosan beszél Önnel.'];
        %     % uniform background
        %     Screen('FillRect', win, backGroundColor);
        %     % draw block-starting text
        %     DrawFormattedText(win, blockEndText, 'center', 'center', textColor);
        %     Screen('Flip', win);
        %
        %     % approximate wait time
        %
        %     % wait for key press to start
        %     while 1
        %         [keyIsDownExp, ~, keyCodeExp] = KbCheck(KbIdxExp);
        %         % experimenter key down
        %         if keyIsDownExp
        %             % if subject is ready to start
        %             if find(keyCodeExp) == keys.go
        %                 break;
        %             % if abort was requested
        %             elseif find(keyCodeExp) == keys.abort
        %                 abortFlag = 1;
        %                 break;
        %             end
        %         end
        %     end
        %     if abortFlag
        %         if triggers
        %             % ppdev_mex('Close', 1);
        %             IOPort('Close',TB);
        %         end
        %         ListenChar(0);
        %         Priority(0);
        %         RestrictKeysForKbCheck([]);
        %         PsychPortAudio('Close');
        %         Screen('CloseAll');
        %         ShowCursor(screenNumber);
        %         return;
        %     end
        %
        %
        % if last block ended now
    elseif block == noBlocks

        % user message
        disp([char(10), 'The task has ended!']);

        % block ending text
        blockEndText = ['End of task!\n',...
            'In the last block you had ', num2str(round(blockAcc, 2)), '% correct responses.\n\n',...
            'Great effort! This window closes soon...'];
        % uniform background
        Screen('FillRect', win, backGroundColor);
        % draw block-starting text
        DrawFormattedText(win, blockEndText, 'center', 'center', textColor);
        Screen('Flip', win);

        WaitSecs(5);

    end  % if block

    SPL = SPL + 3;
    disp([char(10), 'SPL is set to ', num2str(SPL)]);

end  % block for loop


%% Ending, cleaning up

if triggers
    % ppdev_mex('Close', 1);
    IOPort('Close',TB);
end
ListenChar(0);
Priority(0);
RestrictKeysForKbCheck([]);
PsychPortAudio('Close');
Screen('CloseAll');
ShowCursor(screenNumber);

return
end