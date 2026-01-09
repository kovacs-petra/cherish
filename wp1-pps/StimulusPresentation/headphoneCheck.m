function headphoneCheck(init,fs)
% Usage: headphoneCheck(init,fs)
% Play a sound on the left, right, middle, and have the listener respond
% which side they heard it on (left arrow, right arrow, up arrow,
% respectively)
%
% Inputs:
% init - Psychtoolbox needs to be initiated? 1 - yes, 0 - no, the function that called it has done that already
% fs   - sampling frequency 

%% Generate the sounds (triangle waves)
addpath('..\StimulusGeneration\sig_triwave.m');
testDur = 2;
signal = sig_triwave(400,fs,testDur);
left = [signal; zeros(1,length(signal))];
right = [zeros(1,length(signal)); signal];
middle = [signal; signal];
testSounds = {left; right; middle}; % fixed order

if init

    %% Initialize PTB
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

    % screen params, screen selection
backGroundColor = [0 0 0];
textColor = [255 255 255];

    % open stimulus window
[win, rect] = Screen('OpenWindow', S.screen_number, backGroundColor);


    fs = 48000;
    pahandle = PsychPortAudio('Open', S.sound_device_idx, 1, 1, fs, []);
    PsychDefaultSetup(2);
    Screen('Preference', 'SkipSyncTests', 1);

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
end % if init

%% Define keys
    % Define the specific keys we use
    KbName('UnifyKeyNames');
    keys = struct;
    keys.abort = KbName('ESCAPE');
    keys.go = KbName('SPACE');
    keys.left = KbName('LeftArrow');
    keys.right = KbName('RightArrow');
    keys.up = KbName('UpArrow');
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
    ShowCursor(S.screen_number);
    return;
end

% user message
disp([char(10), 'Subject signalled she/he is ready, we start headphone check']);

% respInt = 4;
success = 0;
for check = 1:3
    Screen('CopyWindow', fixCrossWin, win);
    Screen('DrawingFinished', win);
    Screen('Flip', win);
    PsychPortAudio('FillBuffer', pahandle, testBuffer(check)); % bufferdata: 2xN
    PsychPortAudio('Start', pahandle, 1, [], 1);

    % Record a response and write it in a user msg
    while 1
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown % if subject responded
            if find(keyCode) == testResp(check) % and it was the correct button
                disp(strcat('Side ',num2str(check), ': Correct response.'));
                success = success+1;
                WaitSecs(0.5);
                break;
            else
                disp(strcat('Side ',num2str(check), ': Incorrect response.'));
                WaitSecs(0.5);
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
    % disp('Press ENTER to continue or ESC to abort.');
    WaitSecs(2);
    % while 1
    %     [keyIsDown, ~, keyCode] = KbCheck;
    %     if keyIsDown
    %         % if subject is ready to start
    %         if find(keyCode) == keys.go
    %             break;
    %         elseif find(keyCode) == keys.abort
    %             % if abort was requested
    %             abortFlag = 1;
    %             break;
    %         end
    %     end
    % end
    % if abortFlag
        ListenChar(0);
        Priority(0);
        RestrictKeysForKbCheck([]);
        PsychPortAudio('Close');
        Screen('CloseAll');
        ShowCursor(S.screen_number);
    %     return;
    % end

end % action after headphone check

return