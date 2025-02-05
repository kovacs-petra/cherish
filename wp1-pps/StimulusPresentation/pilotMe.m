function pilotMe(subNum, noBlocks, triggers)
% Experimental script for CherISH wp-1 pilot
% Adapted (-ing) from SFGmain.m

% try
%% Input checks
% Determine task from subject number (even Ss start with task 1, odds with
% task 2)
if mod(subNum,2) == 0
    taskOrder = [1,2];
else
    taskOrder = [2,1];
end

%% Initialize Psychtoolbox (screen and audio)
PsychDefaultSetup(1);
Screen('Preference', 'SkipSyncTests', 1);
InitializePsychSound(1);

% Get sampling rate from the hrtf file that was used for stimulus
% generation
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

% Select audio device
% device = [];  % system default is our default as well
tmpDevices = PsychPortAudio('GetDevices');
for i = 1:numel(tmpDevices)
    % if strcmp(tmpDevices(i).DeviceName, 'Lautsprecher (RME Fireface UCX ')
    if strcmp(tmpDevices(i).DeviceName, 'Speakers/Headphones (Realtek(R) Audio)')
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
KbName('UnifyKeyNames');
keys = struct;
keys.abort = KbName('ESCAPE');
keys.go = KbName('SPACE');
keys.respTarget = KbName('RETURN');
keys.left = KbName('LeftArrow');
keys.right = KbName('RightArrow');
keys.up = KbName('UpArrow');

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

%% Headphone check
% Play a sound on the left, right, middle, and have the listener respond
% which side they heard it on (left arrow, right arrow, up arrow,
% respectively)
addpath('..\StimulusGeneration\sig_triwave.m');

% Generate the sounds (triangle waves)
testDur = 2;
signal = sig_triwave(400,fs,testDur);
left = [signal; zeros(1,length(signal))];
right = [zeros(1,length(signal)); signal];
middle = [signal; signal];
testSounds = {left; right; middle}; % fixed order
testResp = [keys.left, keys.right, keys.up]; % correct responses

% Fill buffer with the test sounds
testBuffer = [];
for sounds = 1:3
    testBuffer(end+1) = PsychPortAudio('CreateBuffer', [], testSounds{sounds,1});
end

% Write instructions on the screen
testInstr = ['Headphone check. \n We will play 3 sounds one by one. Please indicate ',...
    'which side you heard each sound on: \n left side - left arrow \n right side - right arrow',...
    'middle - up arrow. \n\n Press SPACE to start the headphone check.'];
Screen('FillRect', win, backGroundColor);
DrawFormattedText(win, testInstr, 'center', 'center', textColor);
Screen('Flip', win);
% user message
disp([char(10), 'Showing headphone check instructions right now...']);

% Wait for SPACE then play a sound
abortFlag = 0;
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
    return;
end

% user message
disp([char(10), 'Subject signalled she/he is ready, we start headphone check']);

Screen('CopyWindow', fixCrossWin, win);
Screen('DrawingFinished', win);
respInt = 4;
success = 0;
for check = 1:3
    PsychPortAudio('FillBuffer', pahandle, testBuffer(check)); % bufferdata: 2xN
    testStart = PsychPortAudio('Start', pahandle, 1, [], 1);

    % Record a response and write it in a user msg
    while GetSecs <= testStart+testDur+respInt
        [keyIsDown, respSecs, keyCode] = KbCheck;
        if keyIsDown % if subject responded
            if find(keyCode) == testResp(check) % and it was the correct button
                disp(strcat('Side ',num2str(check), ': Correct response.'));
                success = success+1;
                break;
            else
                disp(strcat('Side ',num2str(check), ': Incorrect response.'));
                break;
            end
        end
    end
end % headphone check

if success == 3
    % display msg to listener: successful headphone check, we go on
    Screen('FillRect', win, backGroundColor);
    DrawFormattedText(win, 'Headphone check successful.', 'center', 'center', textColor);
    Screen('Flip', win);
    WaitSecs(3);
else
    % display blank window to listener
    % user msg: go and check what's wrong, press return to continue
    Screen('FillRect', win, backGroundColor);
    Screen('Flip', win);
    disp('There was a problem with the headphone check, investigate!');
    disp('Press ENTER to continue or ESC to abort.');
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
        return;
    end

end % action after headphone check


    %% Task loop
    for taskIdx = 1:2
        task = taskOrder(taskIdx);
        % Look for stimArray file in the subject's folder corresponding to task
        stimArrayName = ['stimArray', num2str(task), '.mat'];
        stimArrayFileStruct = dir(stimArrayName);

        % if there was no stimArray file or there were multiple
        if isempty(stimArrayFileStruct) || length(stimArrayFileStruct)~=1
            error(['Either found too many or no stimArrayFile at ./subject', num2str(subNum), '/stimArray*.mat !!!']);
        else
            stimArrayFile = [stimArrayFileStruct.folder, '/', stimArrayFileStruct.name];
        end

        %% Load/set params, stimuli, check for conflicts
        % user message
        disp([char(10), 'Loading params and stimuli, checking ',...
            'for existing files for the subject']);

        % a function handles all stimulus sorting to blocks and potential conflicts
        % with earlier recordings for same subject
        [stimArray, ~, ~,...
            startBlockNo, blockIdx, trialIdx,...
            logVar, subLogF, returnFlag,...
            logHeader, stimTypes] = handleParams(subNum, stimArrayFile, noBlocks, task);

        % if there was early return from expParamsHandler.m, abort
        if returnFlag
            return
        end

        % Collect column numbers from the resulting stimArray
        %%%%%%% HARD-CODED VALUES %%%%%%%
        targetColumn = 10;
        audioColumn = 12;
        stimTypeColumn = 13;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % user message
        disp([char(10), 'Ready to start the experiment']);

        %% Stimulus features for triggers + logging
        totalDurCue = cell2mat(logVar(2:end, strcmp(logHeader, 'totalDur'))); % duration of cue (w/o target or gap)
        totalDurWithTarget = totalDurCue + 2/1000;
        durStatOnset = cell2mat(logVar(2:end, strcmp(logHeader, 'durStatOnset'))); % required for EEG triggering
        durStatOffset = cell2mat(logVar(2:end, strcmp(logHeader, 'durStatOffset'))); % for EEG triggering

        % user message
        disp([char(10), 'Extracted stimulus features']);

        %% Triggers
        % basic triggers for trial start, sound onset and response
        trig = struct;
        trig.trialStart = 200;
        trig.playbackStart = 210;
        trig.resp = 220;
        trig.blockStart = 100;
        trig.motionOnset = 280;
        trig.motionOffset = 290;

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


        % set ITI
        iti = 2;  % in secs

        % response time interval
        respInt = 1;

        % response variables preallocation
        respTime = nan(size(stimArray, 1), 1);
        accuracy = respTime;

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
            if ~exist(TB, 'var')
                TB = IOPort('OpenSerialPort', 'COM3'); % replace with number of COM port % test
            end

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
        instrText = ['Now you will hear sounds moving in different directions. \n',...
            'After some of these sounds, you will sometimes hear a brief target sound. \n\n',...
            'Detect the target sound as quickly and accurately as possible: press ENTER if you heard it. \n',...
            'Don''t press anything if you didn''t hear the target sound. \n\n',...
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

        SPLoptions = [-20 -33 -30 -30];
        SPLidx = 1; % changes in the next block
        SPL = SPLoptions(SPLidx);

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
                audioData = stimArray{trial, audioColumn};
                % Intensity roving:
                audioData = 10^(SPL/20)*audioData;
                buffer(end+1) = PsychPortAudio('CreateBuffer', [], audioData);
            end

            % counter for trials in given block
            trialCounterForBlock = 0;

            % user message
            disp([char(10), 'Buffered all stimuli for block ', num2str(block),...
                ', showing block start message']);

            % block starting text
            blockStartText = ['Beginning block ', num2str(block), '. \n\n',...
                'As always in this round, there will be ', num2str(length(trialList)), ' trials in the block.\n\n\n',...
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

                % wait till we are 100 ms from the start of the playback
                while GetSecs-trialStart <= iti-100
                    WaitSecs(0.001);
                end

                % blocking playback start for precision
                startTime = PsychPortAudio('Start', pahandle, 1, trialStart+iti, 1);

                if triggers
                    % playback start trigger
                    IOPort('Write', TB, uint8(trig.playbackStart),0);
                    pause(0.01);
                    IOPort('Write', TB, uint8(0),0);
                    pause(0.01);

                    % stimulus type trigger
                    IOPort('Write',TB,uint8(trig.stimType),0);
                    pause(0.01);
                    IOPort('Write', TB, uint8(0), 0);
                    pause(0.01);

                    % motion onset trigger
                    while 1
                        if GetSecs == startTime + durStatOnset(trial)
                            IOPort('Write', TB, uint8(trig.motionOnset),0);
                            pause(0.01);
                            IOPort('Write', TB, uint8(0),0);
                            pause(0.01);
                            break
                        end
                    end

                    % motion offset trigger
                    while 1
                        if GetSecs == startTime + totalDurCue(trial) - durStatOffset(trial)
                            IOPort('Write', TB, uint8(trig.motionOffset),0);
                            pause(0.01);
                            IOPort('Write', TB, uint8(0),0);
                            pause(0.01);
                            break
                        end
                    end
                end

                % user message
                disp(['Audio started at ', num2str(startTime-trialStart), ' secs after trial start']);
                disp(['(Target ITI was ', num2str(iti), ' secs)']);

                % Wait for response 1 (target detection)
                respFlag = 0;

                while GetSecs <= startTime+totalDurWithTarget(trial)+respInt
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
                            if triggers
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
                logVar(trial+1,strcmp(logHeader, 'SPL')) = {SPL};
                logVar(trial+1,strcmp(logHeader, 'respTime')) = {respTime(trial)};
                logVar(trial+1,strcmp(logHeader, 'accuracy')) = {accuracy(trial)};

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
                    'Your average reaction time in this block was ',num2str(round(blockRT, 2)), ' ms. \n\n\n',...
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
            elseif block == noBlocks % if this was the last block in the task
                % user message
                disp([char(10), 'The task has ended.']);
                if taskIdx == 1 % if this is the 1st task
                    blockEndText = ['End of round 1. \n',...
                        'Your average reaction time in the last block was ',num2str(round(blockRT, 2)), ' ms. \n\n', ...
                        'Let''s take a break!'];
                else % if it's the 2nd task
                    blockEndText = ['End of round 2. \n',...
                        'Your average reaction time in the last block was ',num2str(round(blockRT, 2)), ' ms. \n\n', ...
                        'This is the end of the experiment. This window closes soon...'];
                end
                % uniform background
                Screen('FillRect', win, backGroundColor);
                % draw block-starting text
                DrawFormattedText(win, blockEndText, 'center', 'center', textColor);
                Screen('Flip', win);

                WaitSecs(5);

            end  % if block

            SPLidx = SPLidx+1; % changes in the next block
            SPL = SPLoptions(SPLidx);
            disp([char(10), 'SPL is set to ', num2str(SPL)]); %#ok<*CHARTEN>

        end  % block for loop

    end % task for loop

    % catch
    %     if triggers
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
    if triggers
        % ppdev_mex('Close', 1);
        IOPort('Close',TB);
    end
    ListenChar(0);
    Priority(0);
    RestrictKeysForKbCheck([]);
    PsychPortAudio('Close');
    Screen('CloseAll');
    ShowCursor;

    return
end