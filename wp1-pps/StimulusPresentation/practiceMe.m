function practiceMe(noTrials)
% Familiarization phase for CherISH wp1 pilot
% After a headphone check, the listener practices target detection.

% Feedback: accurate or inaccurate response; RT of target
% detection

%% Headphone check
flag.headphoneCheck = 0;

if flag.headphoneCheck
    headphoneCheck;
end

%% Find and load stimArray file
% stimArrayFileStruct = dir('stimArrayPractice.mat');
stimArrayFileStruct = dir('stimArray.mat');

% if there was no stimArray file or there were multiple
if isempty(stimArrayFileStruct) || length(stimArrayFileStruct)~=1
    error(['Either found too many or no stimArrayFile at ./subject', num2str(subNum), '/stimArray*.mat !!!']);
else
    stimArrayFile = [stimArrayFileStruct.folder, '/', stimArrayFileStruct.name];
end

load(stimArrayFile, 'stimArray');
audioColumn = 17;
totalDurColumn = 3;
targetColumn = 11;

% get number of trials from stimuli array
trialNo = size(stimArray, 1);

%% Randomize trial order
trialIdx = 1:trialNo;
trialIdx = trialIdx(randperm(length(trialIdx)));

stimArray = [stimArray, num2cell(trialIdx')];
% sort into trial order
[stimArray, ~] = sortrows(stimArray, size(stimArray, 2));

% user message
disp([char(10), 'Sorted stimuli into final stimArray']);

%% Psychtoolbox initialization

% General init (AssertOpenGL, 'UnifyKeyNames')
PsychDefaultSetup(1);
Screen('Preference', 'SkipSyncTests', 1);

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

%% Keyboard and screen
% Define the specific keys we use
keys = struct;
keys.abort = KbName('ESCAPE');
keys.go = KbName('SPACE');
keys.respTarget = KbName('RETURN');
% keys.respDistancePPS = KbName('1');
% keys.respDistanceARS = KbName('2');
% keys.respDistanceEPS = KbName('3');
% distanceResps = [keys.respDistancePPS keys.respDistanceARS keys.respDistanceEPS];

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
% ifi = Screen('GetFlipInterval', win);
% set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% Setup the text type for the window
Screen('TextFont', win, 'Arial');
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

%% open PsychPortAudio device for playback
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

% set ITI and response interval  and feedback text duration (in secs)
% iti = 3; 
respInt = 1;
fbDur = 2; 

% set flag for aborting experiment
abortFlag = 0;
% hide mouse
HideCursor(screenNumber);
% suppress keyboard input to command window
ListenChar(-1);
% realtime priority
Priority(1);

% Log accuracy and RT
accuracy = zeros(1,noTrials);
respTime = accuracy;

% SPL = -30; % volume ctrl

% user message
disp([char(10), 'Initialized psychtoolbox basics, opened window, ',...
    'started PsychPortAudio device']);

%% Instructions phase

% instructions text
instrText = ['You will hear sounds moving in different directions. \n',...
    'After some of these sounds, you will sometimes hear a brief target sound. \n\n',...
    'Detect the target sound as quickly and accurately as possible: press ENTER if you heard it. \n',...
    'Don''t press anything if you didn''t hear the target sound. \n\n',...
    'Respond within 1 second after the target. \n\n', ...
    '<Press SPACE to continue>'];

% write instructions on the screen
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
            abortFlag = 1;
            break;
        end
    end
end

if abortFlag
    disp('Terminating at user''s or subject''s request');
    ListenChar(0);
    Priority(0);
    RestrictKeysForKbCheck([]);
    PsychPortAudio('Close');
    Screen('CloseAll');
    ShowCursor(screenNumber);
    return;
end

% user message
disp([char(10), 'Subject is ready to practice, we are starting...']);

%% Loop of playing requested stimuli

for trial = 1:noTrials  

    % display fixation cross
    % background with fixation cross, get trial start timestamp
    Screen('CopyWindow', fixCrossWin, win);
    Screen('DrawingFinished', win);
    trialStart = Screen('Flip', win);

    % load a random but corresponding stimulus from the stimArray into
    % buffer
    randIdx = randi(length(stimArray),1);
    
    soundOutput = cell2mat(stimArray(randIdx,audioColumn));
    % soundOutput = soundOutput*10^(SPL/20);
    totalDur = cell2mat(stimArray(randIdx,totalDurColumn));
    buffer = PsychPortAudio('CreateBuffer', pahandle, soundOutput);
    PsychPortAudio('FillBuffer', pahandle, buffer);
    % play stimulus - blocking start
    startTime = PsychPortAudio('Start', pahandle, 1, trialStart+respInt, 1);

            % Wait for response (target detection)
        respFlag = 0;

        while GetSecs <= startTime+totalDur+respInt % response interval hasn't ended
            [keyIsDown, respSecs, keyCode] = KbCheck;
            if keyIsDown % if subject responded
                if find(keyCode) == keys.respTarget % and it was the detection button
                    respFlag = 1;
                    if cell2mat(stimArray(randIdx,targetColumn)) == 1 % and there is a target in the trial
                        accuracy(trial) = 1;
                    else
                        accuracy(trial) = 0;
                    end
                    break;
                elseif find(keyCode) == keys.abort % if it was the escape button
                    abortFlag = 1;
                    break;
                end
            else % if subject did not respond
                if cell2mat(stimArray(randIdx,targetColumn)) == 0 % and there is no target in the trial
                    accuracy(trial) = 1;
                else
                    accuracy(trial) = 0;
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
            return;
        end

        % switch visual right when the audio finishes
        % Screen('Flip', win, startTime+totalDurWithTarget(trial)-0.5*ifi);

        % response time into results variable
        if respFlag
            respTime(trial) = 1000*(respSecs-(startTime+totalDur));
        end

        % user messages
        if cell2mat(stimArray(randIdx,targetColumn)) == 1 % if it's a target trial
            if isnan(respTime(trial))
                disp('Subject did not respond in time');
            else
                disp(['RT: ', num2str(respTime(trial))]);
            end
        end

        % participant feedback
        if accuracy(trial) == 1
            if cell2mat(stimArray(randIdx,targetColumn)) == 1 % there was a target
                accMsg = 'Correct detection!';
                rtMsg = ['Reaction time: ', num2str(respTime(trial)), ' ms'];
                Screen('FillRect', win, backGroundColor);
                DrawFormattedText(win, [accMsg, '\n', rtMsg], 'center', 'center', textColor);
                Screen('Flip', win);
                WaitSecs(fbDur);
            else % if there was no target
                accMsg = 'Correct, no target!';
                Screen('FillRect', win, backGroundColor);
                DrawFormattedText(win, accMsg, 'center', 'center', textColor);
                Screen('Flip', win);
                WaitSecs(fbDur);
            end

        elseif accuracy(trial) == 0
            if cell2mat(stimArray(randIdx,targetColumn)) == 1 % there was a target
                accMsg = 'Incorrect, there was a target!';
                Screen('FillRect', win, backGroundColor);
                DrawFormattedText(win, accMsg, 'center', 'center', textColor);
                Screen('Flip', win);
                WaitSecs(fbDur);
            else % if there was no target
                accMsg = 'Incorrect, there was no target!';
                Screen('FillRect', win, backGroundColor);
                DrawFormattedText(win, accMsg, 'center', 'center', textColor);
                Screen('Flip', win);
                WaitSecs(fbDur);
            end

        end % feedback loop

end % trial loop

DrawFormattedText(win, 'End of practice. This window closes soon...', 'center', 'center', textColor);
Screen('Flip', win);
WaitSecs(3);

ListenChar(0);
Priority(0);
RestrictKeysForKbCheck([]);
PsychPortAudio('Close');
Screen('CloseAll');


