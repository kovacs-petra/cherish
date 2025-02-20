function headphoneCheck
% Headphone check
% Play a sound on the left, right, middle, and have the listener respond
% which side they heard it on (left arrow, right arrow, up arrow,
% respectively)


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