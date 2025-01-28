function practiceMe
% Familiarization phase for CherISH wp1 pilot
% The listener practices target detection, and especially distance judgements.

% Feedback: accurate or inaccurate distance judgement; RT of target
% detection
% Should happen at the beginning of each run (when sounds go back from loud
% to soft), not in a separate script

%% Find and load stimArray file
stimArrayFileStruct = dir('stimArrayPractice.mat');

% if there was no stimArray file or there were multiple
if isempty(stimArrayFileStruct) || length(stimArrayFileStruct)~=1
    error(['Either found too many or no stimArrayFile at ./subject', num2str(subNum), '/stimArray*.mat !!!']);
else
    stimArrayFile = [stimArrayFileStruct.folder, '/', stimArrayFileStruct.name];
end

load(stimArrayFile, 'stimArray');

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

%% Keyboard and screen
% Define the specific keys we use
keys = struct;
keys.abort = KbName('ESCAPE');
keys.go = KbName('SPACE');
keys.respTarget = KbName('RETURN');
keys.respDistancePPS = KbName('1');
keys.respDistanceARS = KbName('2');
keys.respDistanceEPS = KbName('3');
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

% set ITI
iti = 0.6;  % in secs

% set flag for aborting experiment
abortFlag = 0;
% hide mouse
HideCursor(screenNumber);
% suppress keyboard input to command window
ListenChar(-1);
% realtime priority
Priority(1);

SPL = -30; % volume ctrl

% user message
disp([char(10), 'Initialized psychtoolbox basics, opened window, ',...
    'started PsychPortAudio device']);

%% Instructions phase

% instructions text
instrText = ['For practice, you can listen to what tones sound like when they end at different distances. \n',...
    'To request a sound ending at a specific distance category, press one of these keys \n',...
    'on the number pad: \n\n',...
    '1: uncomfortably close to my face \n',...
    '2: not extremely close, but still within arm''s reach \n',...
    '3: quite far, out of reach \n \n',...
    'Press a key (1, 2 or 3) to listen.'];

taskText = ['Play me a sound that ends up: \n', '1: uncomfortably close to my face \n',...
    '2: not extremely close, but still within arm''s reach \n',...
    '3: quite far, out of reach \n \n',...
    'Press a key (1, 2 or 3) to listen. \n\n',...
    'If you''re done practicing, press ESCAPE.'];

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
        if find(keyCode) == keys.respDistancePPS
            nextTrial = 1;
            break;
        elseif find(keyCode) == keys.respDistanceARS
            % if abort was requested
            nextTrial = 2;
            break;
        elseif find(keyCode) == keys.respDistanceEPS
            nextTrial = 3;
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
if nextTrial==1
    trialMessage = 'PPS';
elseif nextTrial==2
    trialMessage = 'ARS';
elseif nextTrial==3
    trialMessage = 'EPS';
end
disp([char(10), 'Subject requested a trial in ', trialMessage, ', we are starting...']);

%% Loop of playing requested stimuli

while 1  % until abort is requested

    % display fixation cross
    % background with fixation cross, get trial start timestamp
    Screen('CopyWindow', fixCrossWin, win);
    Screen('DrawingFinished', win);
    trialStart = Screen('Flip', win);

    % load a random but corresponding stimulus from the stimArray into
    % buffer
    randIdx = randi(length(stimArray),1);

    while cell2mat(stimArray(randIdx,16)) ~= nextTrial
        randIdx = randi(length(stimArray),1);
    end
    
    soundOutput = cell2mat(stimArray(randIdx,18));
    soundOutput = soundOutput*10^(SPL/20);
    totalDur = cell2mat(stimArray(randIdx,8));
    buffer = PsychPortAudio('CreateBuffer', pahandle, soundOutput');
    PsychPortAudio('FillBuffer', pahandle, buffer);
    % play stimulus - blocking start
    startTime = PsychPortAudio('Start', pahandle, 1, trialStart+iti, 1);

    % wait till playback is over
    WaitSecs('UntilTime', startTime+totalDur);

    % display text for subject about requesting next stimulus
    Screen('FillRect', win, backGroundColor);
    DrawFormattedText(win, taskText, 'center', 'center', textColor);
    Screen('Flip', win);

    % wait for key press to play next stimulus
    while 1
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            % if subject is ready to start
            if find(keyCode) == keys.respDistancePPS
                nextTrial = 1;
                break;
            elseif find(keyCode) == keys.respDistanceARS
                % if abort was requested
                nextTrial = 2;
                break;
            elseif find(keyCode) == keys.respDistanceEPS
                nextTrial = 3;
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
    if nextTrial==1
        trialMessage = 'PPS';
    elseif nextTrial==2
        trialMessage = 'ARS';
    elseif nextTrial == 3
        trialMessage = 'EPS';
    end
    disp([char(10), 'Subject requested a trial in ', trialMessage, ', we are starting...']);

    % only go on to next stimulus when keys are released
    KbReleaseWait;
end


