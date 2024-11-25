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

%% Audio parameters for PsychPortAudio

% sampling rate is derived from the HRTF set, 44100 Hz in this case
fs = 44100;

% Select audio device: this is where I could use o_ptb
device = [];  % system default is our default as well % test

% mode is simple playback
mode = 1;
% reqlatencyclass is set to low-latency
reqLatencyClass = 2;
% 2 channels output
nrChannels = 2;

% user message
disp([char(10), 'Set audio parameters']);

%% Stimulus features for triggers + logging
endSpace = cell2mat(logVar(2:end, strcmp(logHeader, 'space')));
% trajectory = cell2mat(logVar(2:end, strcmp(logHeader, 'trajectory')));
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

% Keyboard params - names %%%%%%%%%%% o_ptb probably better %%%%%%%%%%
KbNameSub = 'Logitech USB Keyboard'; % set
KbNameExp = 'CASUE USB KB'; % set
% detect attached devices
[keyboardIndices, productNames, ~] = GetKeyboardIndices; % test
% define subject's and experimenter keyboards
KbIdxSub = keyboardIndices(ismember(productNames, KbNameSub));
KbIdxExp = keyboardIndices(ismember(productNames, KbNameExp));

% Define the specific keys we use
keys = struct;
keys.abort = KbName('ESCAPE');
keys.go = KbName('SPACE');
keys.respStop = KbName('RETURN');
keys.respDistance = KbName('1','2','3'); % test

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

% set up the question mark (stimulus marking response period) into a
% texture / offscreen window
qMarkWin = Screen('OpenOffscreenWindow', win, backGroundColor, rect);
Screen('BlendFunction', qMarkWin, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('TextSize', qMarkWin, 50);
Screen('DrawText', qMarkWin, ['Distance? \n1: very close \n2: within arm''s reach \n3: out of reach'],...
    xCenter-15, yCenter-15, textColor);

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
respInt1 = 2; % 2 secs from the end of the last moving portion (includes last static portion)
respInt2 = 5; % 5 secs from the end of playback

% response variables preallocation
respTime = nan(size(stimArray, 1), 1);
acc = respTime;

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
instrText = ['A feladat ugyanaz lesz, mint az előző blokkok során - \n',... 
    'jelezze, ha a hangmintában hall egy emelkedő hangsort.\n\n',...
    'Nyomja meg a SPACE billentyűt ha készen áll!'];

% write instructions to text
Screen('FillRect', win, backGroundColor);
DrawFormattedText(win, instrText, 'center', 'center', textColor);
Screen('Flip', win);

% user message
disp([char(10), 'Showing the instructions text right now...']);

% wait for key press to start
while 1
    [keyIsDownSub, ~, keyCodeSub] = KbCheck(KbIdxSub);
    [keyIsDownExp, ~, keyCodeExp] = KbCheck(KbIdxExp);
    % subject key down
    if keyIsDownSub 
        % if subject is ready to start
        if find(keyCodeSub) == keys.go
            break;
        end
    % experimenter key down    
    elseif keyIsDownExp
        % if abort was requested    
        if find(keyCodeExp) == keys.abort
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

SPL = -3; % changes in the next block; values: [-3 0 3]

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
        else
            audioData = stimArray{trial, 12};
        end
    end

    % Intensity roving:
    audioData = 10^(SPL/20)*audioData;

    buffer(end+1) = PsychPortAudio('CreateBuffer', [], audioData');

    % counter for trials in given block
    trialCounterForBlock = 0;    
    
    % user message
    disp([char(10), 'Buffered all stimuli for block ', num2str(block),... 
        ', showing block start message']);    
     
    % block starting text
    blockStartText = ['Kezdhetjük a(z) ', num2str(block), '. blokkot,\n\n',... 
            'mint mindegyik, ez is ', num2str(length(trialList)), ' próbából fog állni.\n\n\n',... 
            'Nyomja meg a SPACE billentyűt ha készen áll!'];
    
    % uniform background
    Screen('FillRect', win, backGroundColor);
    % draw block-starting text
    DrawFormattedText(win, blockStartText, 'center', 'center', textColor);
    Screen('Flip', win);
    % wait for key press to start
    while 1
        [keyIsDownSub, ~, keyCodeSub] = KbCheck(KbIdxSub);
        [keyIsDownExp, ~, keyCodeExp] = KbCheck(KbIdxExp);
        % subject key down
        if keyIsDownSub 
            % if subject is ready to start
            if find(keyCodeSub) == keys.go
                break;
            end
        % experimenter key down    
        elseif keyIsDownExp
            % if abort was requested    
            if find(keyCodeExp) == keys.abort
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
        
        % prepare screen change for response period       
        Screen('CopyWindow', qMarkWin, win);
        Screen('DrawingFinished', win);
        
        % switch visual right when the audio finishes
        respStart = Screen('Flip', win, startTime+stimLength(trial)-0.5*ifi);
        
        % user message
        disp(['Visual flip for response period start was ', num2str(respStart-startTime),... 
            ' secs after audio start (should equal ', num2str(stimLength(trial)), ')']);
        
        % Wait for response 1 (when did the movement stop?)
        respFlag1 = 0;
        while GetSecs-(startTime+stimLength(trial)-lastPortionLength(trial)) <= respInt1
            [keyIsDownSub, respSecs, keyCodeSub] = KbCheck(KbIdxSub);
            [keyIsDownExp, ~, keyCodeExp] = KbCheck(KbIdxExp);
            % subject key down
            if keyIsDownSub
                % if subject responded figure presence/absence
                if find(keyCodeSub) == keys.respStop
                    respFlag1 = 1;
                    % response trigger
                    if triggers
                        IOPort('Write',TB,uint8(trig.resp),0);
                        pause(0.01);
                        IOPort('Write', TB, uint8(0),0);
                        pause(0.01);
                    end
                    break;
                end
            % experimenter key down    
            elseif keyIsDownExp
                % if abort was requested    
                if find(keyCodeExp) == keys.abort
                    abortFlag = 1;
                    break;
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
        
        % response time into results variable
        if respFlag1
            respTime(trial) = 1000*(respSecs-respStart-lastPortionLength(trial));
        end
        
        % user messages
        if isnan(respTime(trial))
            disp('Subject did not respond in time');
        end

        % Wait for response 2 (how far was the sound at the end point?)
        respFlag2 = 0;
        while GetSecs-(startTime+stimLength(trial)) <= respInt2
            [keyIsDownSub, respSecs, keyCodeSub] = KbCheck(KbIdxSub);
            [keyIsDownExp, ~, keyCodeExp] = KbCheck(KbIdxExp);
            % subject key down
            if keyIsDownSub
                % if subject responded 
                if ismember(find(keyCodeSub), keys.respDistance)
                    respFlag2 = 1;
                    % response trigger
                    if triggers
                        IOPort('Write',TB,uint8(trig.resp),0);
                        pause(0.01);
                        IOPort('Write',TB,uint8(0),0);
                        pause(0.01);
                    end
                    break;
                end
            % experimenter key down    
            elseif keyIsDownExp
                % if abort was requested    
                if find(keyCodeExp) == keys.abort
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

        % accuracy
        if respFlag2
            if find(keyCodeSub) == endSpace(trial)
                disp('Subject''s response was accurate');
                acc(trial) = 1;
            else
                disp('Subject made an error');
                acc(trial) = 0;
            end
        end

        % cumulative accuracy in block
        blockAcc = sum(acc(trial-trialCounterForBlock+1:trial), 'omitnan')/trialCounterForBlock*100;
        disp(['Overall accuracy in block so far is ', num2str(blockAcc), '%']);
        
        % accumulating all results in logging / results variable
        logVar(trial+1,strcmp(logHeader, 'SPL')) = SPL;
        logVar(trial+1,strcmp(logHeader, 'accuracy')) = acc(trial);
        logVar(trial+1,strcmp(logHeader, 'respTime')) = respTime(trial);
        logVar(trial+1,strcmp(logHeader, 'iti')) = iti(trial);
        logVar(trial+1,strcmp(logHeader, 'trialStart')) = trialStart;
        logVar(trial+1,strcmp(logHeader, 'soundOnset')) = startTime-trialStart;

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
        blockEndText = ['Vége a(z) ', num2str(block), '. blokknak!\n\n\n',... 
                'Ebben a blokkban a próbák ', num2str(round(blockAcc, 2)), '%-ra adott helyes választ.\n\n\n',... 
                'Nyomja meg a SPACE billentyűt ha készen áll a következő blokkra!'];
        % uniform background
        Screen('FillRect', win, backGroundColor);
        % draw block-starting text
        DrawFormattedText(win, blockEndText, 'center', 'center', textColor);   
        Screen('Flip', win);
        % wait for key press to start
        while 1
            [keyIsDownSub, ~, keyCodeSub] = KbCheck(KbIdxSub);
            [keyIsDownExp, ~, keyCodeExp] = KbCheck(KbIdxExp);
            % subject key down
            if keyIsDownSub 
                % if subject is ready to start
                if find(keyCodeSub) == keys.go
                    break;
                end
            % experimenter key down    
            elseif keyIsDownExp
                % if abort was requested    
                if find(keyCodeExp) == keys.abort
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
        blockEndText = ['Vége a feladatnak!\n',...
            'Az utolsó blokkban a próbák ', num2str(round(blockAcc, 2)), '%-ra adott helyes választ.\n\n',...
            'Köszönjük a részvételt!'];       
        % uniform background
        Screen('FillRect', win, backGroundColor);
        % draw block-starting text
        DrawFormattedText(win, blockEndText, 'center', 'center', textColor);  
        Screen('Flip', win);
        
        WaitSecs(5);
        
    end  % if block       
         
SPL = SPL + 3;

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