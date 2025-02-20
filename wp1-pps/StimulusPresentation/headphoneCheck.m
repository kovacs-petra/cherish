function headphoneCheck
% Headphone check
% Play a sound on the left, right, middle, and have the listener respond
% which side they heard it on (left arrow, right arrow, up arrow,
% respectively)

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

%% Generate the sounds (triangle waves)
addpath('..\StimulusGeneration\sig_triwave.m');
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
    'which side you heard each sound on: \n left side - left arrow \n right side - right arrow \n',...
    'middle - up arrow \n\n Press SPACE to start the headphone check.'];
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

return