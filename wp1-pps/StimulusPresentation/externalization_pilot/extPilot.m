function extPilot(subNum)
% Usage: extPilot(subNum)
%
% Externalization pilot for CherISH DC6 WP1
%
% Input:
%   - subNum:   subject number, integer
%
% - Distance judgement task with a GUI
% - Stimuli are presented on pages, every page contains 3 sounds that can be
% played by the listener. Each sound on a page has a different distance
% (out of 4 options: in the head, 20 cm, 60 cm, 1 m)
% - Source intensity changes every 4 pages and the listener is informed about
% this.
% - First half of the experiment has variable f0, second half constant f0, or
% the other way round.
% - Stimuli have to be arranged in an f0s X distances X azimuths cell
% array, with each cell containing the audio (2 channels) in position 1 and
% the specific distance in position 2
%
% Author: Petra Kovacs, 2025

%% Set trial number info
nTrials = 128; % how many trials in total
nTrialsPerF0 = nTrials/2; % how many trials in an f0 condition
nTrialsPerSI = nTrialsPerF0/2; % how many trials per source intensity condition
nTrialsPerPage = 4; % how many trials on a page

%% Initialize Psychtoolbox (screen and audio)
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1);
InitializePsychSound(1);

% Select audio device
device = [];
% tmpDevices = PsychPortAudio('GetDevices');
% for i = 1:numel(tmpDevices)
%     % if strcmp(tmpDevices(i).DeviceName, 'Headphones (Conexant HD Audio headphone)')
%     if strcmp(tmpDevices(i).DeviceName, 'Speakers/Headphones (Realtek(R) Audio)')
%         device = tmpDevices(i).DeviceIndex;
%     end
% end
% mode is simple playback
mode = 1;
% reqlatencyclass is set to low-latency
reqLatencyClass = 2;
% 2 channels output
nrChannels = 2;
fs = 48e3;

% user message
disp([char(10), 'Set audio parameters']);

% Define the specific keys we use
KbName('UnifyKeyNames');
keys = struct;
keys.abort = KbName('ESCAPE');
keys.go = KbName('SPACE');

% restrict keys to the ones we use
keysFields = fieldnames(keys);
keysVector = zeros(1, length(keysFields));
for f = 1:length(keysFields)
    keysVector(f) = keys.(keysFields{f});
end
RestrictKeysForKbCheck(keysVector);

% Force costly mex functions into memory to avoid latency later on
GetSecs; WaitSecs(0.1); KbCheck();

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

%% Open the stimuli structure
stimArrayName = 'stimStruct.mat';
stimArrayFileStruct = dir(stimArrayName);
stimArrayFile = [stimArrayFileStruct.folder, '/', stimArrayFileStruct.name];
load(stimArrayFile, 'stimStruct');

%% Initialize the graphical interface for Psychtoolbox
colBG=[200 200 200];
screens=Screen('Screens');
screenNumber=max(screens);

% ignore screen warning
% PsychDefaultSetup(2); % makes sure Screen is functional and unifies keyCodes across OS
[w,rect] = Screen('OpenWindow',screenNumber, colBG);
% debugRect = [10 10 400 1200];
%     [w,rect] = Screen('OpenWindow',screenNumber,colBG,debugRect);
%     ShowCursor;

[screenXpix,screenYpix] = Screen('WindowSize', w);
[x_center,y_center] = RectCenter(rect);

Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%% Listener instruction
instruction = [...
    'Abstände von Tönen einschätzen: \n\n\n',...
    'Pro Seite gibt es vier Töne, die Sie durch Mausklick auf die Lautsprechersymbole links im Bild so oft Sie wollen anhören können. \n',...
    'Rechts neben den Lautsprechersymbolen gibt es pro Ton eine Leiste, auf der Sie den Abstand des Tons per Mausklick stufenlos angeben können. \n',...
    'Die Töne können links oder rechts sein, oder sogar innerhalb Ihrem Kopf. \n\n',...
    'Wenn Sie die vier Töne ausreichend verglichen haben und mit Ihrer Einschätzung zufrieden sind, drücken Sie die Leertaste, um fortzufahren. \n',...
    'Manche Seiten habe leisere, manche lautere Töne. Außerdem kann die Tonquelle entweder gleich bleiben, oder sich ändern. \n', ...
    'Diese Informationen werden Sie bekommen, bitte beachten Sie darauf. \n',...
    'Das Experiment besteht insgesamt aus ' num2str(nTrials/nTrialsPerPage) ' Seiten. \n\n\n',...
    'Haben Sie noch Fragen? \n\n',...
    'Zum Starten die Leertaste drücken.'];

DrawFormattedText(w,instruction,.2*x_center,'center',[0 0 0],120,0,0,1.5);
Screen('Flip',w);

%% define colors, sizes and positions
% img of head
% path_head = "C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\head_trialversion.png";
path_head = "C:\Users\experimentator.KFS\Documents\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\head_trialversion.png";
[img, ~, alpha] = imread(path_head);
texture = Screen('MakeTexture', w, img);
widthimg=1500;
heightimg=100;

colG=[50,200,50];
colB=[50,50,200];
colR=[200,50,50];
colY=[200,200,50];

defaultTxtSize = Screen('TextSize',w);

leftbordersound=100;
leftborderscale=leftbordersound+150;
sizebutton=100;
widthscale=1500;
colscale=[0,0,0];
startheightbutton=100;
% startheightbutton=200;
% distancebuttons=320;
distancebuttons=250;
heightscale=10;
widthslider=5;
heightslider=50;
% sizeimg=150; % or if you want original size: % sizeimg=size(img); sizeimg(1);

heightB1=startheightbutton;
heightB2=startheightbutton+distancebuttons;
heightB3=startheightbutton+2*distancebuttons;
heightB4=startheightbutton+3*distancebuttons;

heightS1=startheightbutton+0.5*sizebutton-0.5*heightscale;
heightS2=startheightbutton+distancebuttons+0.5*sizebutton-0.5*heightscale;
heightS3=startheightbutton+2*distancebuttons+0.5*sizebutton-0.5*heightscale;
heightS4=startheightbutton+3*distancebuttons+0.5*sizebutton-0.5*heightscale;

posBut1=[leftbordersound, heightB1, leftbordersound+sizebutton, heightB1+sizebutton];
posBut2=[leftbordersound, heightB2, leftbordersound+sizebutton, heightB2+sizebutton];
posBut3=[leftbordersound, heightB3, leftbordersound+sizebutton, heightB3+sizebutton];
posBut4=[leftbordersound, heightB4, leftbordersound+sizebutton, heightB4+sizebutton];

posimg1=[leftborderscale, heightS1-heightimg, leftborderscale+widthimg, heightS1+heightimg];
posimg2=[leftborderscale, heightS2-heightimg, leftborderscale+widthimg, heightS2+heightimg];
posimg3=[leftborderscale, heightS3-heightimg, leftborderscale+widthimg, heightS3+heightimg];
posimg4=[leftborderscale, heightS4-heightimg, leftborderscale+widthimg, heightS4+heightimg];

% sound symbol img
% path_symb = "C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\soundsymbol.png";
path_symb = "C:\Users\experimentator.KFS\Documents\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\soundsymbol.png";
[imgSymb, ~, alphasymb] = imread(path_symb);
imgSymb(:, :, 4) = alphasymb;
textureS = Screen('MakeTexture', w, imgSymb);

%% experiment
% wait for key after instructions
while ~KbCheck; end
WaitSecs(0.6)

%% Frequencies
fmin = 310;
fmax = 2*fmin; % 1 octave
nF0 = 10; % unique f0 values
Freqs = round(linspace(fmin,fmax,nF0),-1);

% Problematic F0s: 380 and 410 Hz
goodIdxs = [1,2,5,6,7,8,9,10];
stimStruct.stim = stimStruct.stim(goodIdxs,:,:);
Freqs=Freqs(Freqs~=380); Freqs=Freqs(Freqs~=410);
nF0 = 8;
% Freqs = Freqs(randperm(numel(Freqs))); % randomize
constantf0 = Freqs(round(nF0/2));
% constantf0 = 410;

%% Azimuths
Azimuths = [90,-90];

%% Conditions: f0 and source intensity
constf0 = zeros(1,nTrialsPerF0/nTrialsPerPage); varf0 = ones(1,nTrialsPerF0/nTrialsPerPage);
f0conds = [constf0 varf0];
if randi([0 1],1) == 0
    f0conds = flip(f0conds); % randomize order
end

lowInt = zeros(1,nTrialsPerSI/nTrialsPerPage); highInt = ones(1,nTrialsPerSI/nTrialsPerPage);
sourceIntConds = repmat([lowInt highInt],1,length(unique(f0conds)));
if mod(subNum,2) == 0 
    sourceIntConds = flip(sourceIntConds); % randomize order
end

%% Counters
resp=[];
uu=1; % trial counter

for pp = 1:(nTrials/nTrialsPerPage) % Page

    % Handle f0 condition
    if f0conds(pp) == 0 % constant f0 throughtout the experiment
        f01=constantf0; f02=constantf0; f03=constantf0; f04=constantf0;
        f1=find(Freqs==f01); f2=find(Freqs==f02); f3=find(Freqs==f03);
        f4=find(Freqs==f04);

    elseif f0conds(pp) == 1 % varying f0 (even within page)
        % f0 chosen randomly for each trial in a page
        f1 = randi(nF0,1); f01 = Freqs(f1);
        f2 = randi(nF0,1); f02 = Freqs(f2);
        f3 = randi(nF0,1); f03 = Freqs(f3);
        f4 = randi(nF0,1); f04 = Freqs(f4);
    end

    % Handle source intensity condition
    if sourceIntConds(pp) == 0
        SPL = db(0.2);
    else
        SPL = db(1);
    end

    % azi chosen randomly for each trial in a page
    azi1 = randi(length(Azimuths),1);
    azi2 = randi(length(Azimuths),1);
    azi3 = randi(length(Azimuths),1);
    azi4 = randi(length(Azimuths),1);

    % distance counterbalanced within a page
    distFix = randperm(numel(1:3)); % 0, 0.2 and 1 m all have to be presented, in random order.
    % They are found in columns 1:3 of
    % the stim cell
    distRand = randi([4,7],1); % One instance of the distances 0.3-0.9 m has to be present

    % Load sound
    cellG = stimStruct.stim{f1,distFix(1),azi1};
    cellB = stimStruct.stim{f2,distFix(2),azi2};
    cellR = stimStruct.stim{f3,distFix(3),azi3};
    cellY = stimStruct.stim{f4,distRand,azi4};

    soundG = cellG{1}*10^(SPL/20); distG = cellG{2};
    soundB = cellB{1}*10^(SPL/20); distB = cellB{2};
    soundR = cellR{1}*10^(SPL/20); distR = cellR{2};
    soundY = cellY{1}*10^(SPL/20); distY = cellY{2};

    % phase to listen to sounds and choose response
    respG=[];
    respB=[];
    respR=[];
    respY=[];

    % clear screen
    Screen('FillRect', w ,colBG, [0 0 screenXpix screenYpix]);
    Screen('Flip', w, 0, 1);

    % write source intensity information if it has changed
    if pp == 1 || sourceIntConds(pp) ~= sourceIntConds(pp-1)
        if sourceIntConds(pp) == 0
            intensityInfo = 'Die nächsten Seiten haben LEISE Töne ';
            disp('Low intensity pages.')
        else
            intensityInfo = 'Die nächsten Seiten haben LAUTE Töne ';
            disp('High intensity pages.')
        end

        if pp == 1 || f0conds(pp) ~= f0conds(pp-1)
            if f0conds(pp) == 0
                f0Info = 'und die Tonquelle bleibt immer gleich...';
                disp('Constant F0 starting.')
            else
                f0Info = 'und die Tonquelle ändert sich immer...';
                disp('Variable F0 starting.')
            end
        end

    % display info
    Screen('TextSize',w,30);
    DrawFormattedText(w,[intensityInfo,'\n',f0Info],'center','center',[0 0 0],120,0,0,1.5);
    Screen('Flip', w, 0, 1);
    WaitSecs(5);
    Screen('TextSize',w,defaultTxtSize);
    end


    % % write f0 information if it has changed
    % if pp == 1 || f0conds(pp) ~= f0conds(pp-1)
    %     if f0conds(pp) == 0
    %         f0Info = 'Die Tonquelle bleibt jetzt immer gleich...';
    %         disp('Constant F0 starting.')
    %     else
    %         f0Info = 'Die Tonquelle ändert sich jetzt immer...';
    %         disp('Variable F0 starting.')
    %     end
    % end    

    Screen('FillRect', w ,colBG, [0 0 screenXpix screenYpix]);

    % draw rectangles as sound buttons
    Screen('FillRect', w ,colG, posBut1);
    Screen('FillRect', w ,colB, posBut2);
    Screen('FillRect', w ,colR, posBut3);
    Screen('FillRect', w ,colY, posBut4);

    % draw sound symbols in buttons
    Screen('DrawTexture', w, textureS, [], posBut1);
    Screen('DrawTexture', w, textureS, [], posBut2);
    Screen('DrawTexture', w, textureS, [], posBut3);
    Screen('DrawTexture', w, textureS, [], posBut4);

    % draw imgs of heads
    Screen('DrawTexture', w, texture, [], posimg1);
    Screen('DrawTexture', w, texture, [], posimg2);
    Screen('DrawTexture', w, texture, [], posimg3);
    Screen('DrawTexture', w, texture, [], posimg4);

    % draw rectangle as slider line
    Screen('FillRect', w ,colscale,[leftborderscale, heightS1, leftborderscale+widthscale, heightS1+heightscale]);
    Screen('FillRect', w ,colscale,[leftborderscale, heightS2, leftborderscale+widthscale, heightS2+heightscale]);
    Screen('FillRect', w ,colscale,[leftborderscale, heightS3, leftborderscale+widthscale, heightS3+heightscale]);
    Screen('FillRect', w ,colscale,[leftborderscale, heightS4, leftborderscale+widthscale, heightS4+heightscale]);

    % draw progress text
    progress = ['Seite ' num2str(fix(uu/4)+1) ' von ' num2str(nTrials/nTrialsPerPage),...
        '. Weiter per Leertaste.'];
    DrawFormattedText(w,progress,leftborderscale,heightB4+heightimg*2,[0 0 0],120,0,0,1.5);
    Screen('Flip', w, 0, 1);

    while not(KbCheck && isempty(respG)+isempty(respB)+isempty(respR)+isempty(respY) == 0) %check keyboard has not been pressed
        [x,y,button] = GetMouse(w);
        if any(button)
            if x > leftborderscale && x < leftborderscale+widthscale
                if y > heightB1 && y < heightB1+sizebutton
                    if isempty(respG)
                        Screen('FillRect', w ,colG,[x, heightS1-heightslider/2+heightscale/2, x+widthslider, heightS1+heightslider/2+heightscale/2]);
                        Screen('Flip', w, 0, 1);
                        respG=x-leftborderscale;
                    else
                        Screen('FillRect', w ,colBG, [leftborderscale,heightS1-heightimg,leftborderscale+widthscale,heightS1+heightimg]);
                        Screen('DrawTexture', w, texture, [], posimg1);
                        Screen('FillRect', w ,colscale,[leftborderscale, heightS1, leftborderscale+widthscale, heightS1+heightscale]);
                        Screen('FillRect', w ,colG,[x, heightS1-heightslider/2+heightscale/2, x+widthslider, heightS1+heightslider/2+heightscale/2]);
                        Screen('Flip', w, 0, 1);
                        respG=x-leftborderscale;
                    end
                elseif y > heightB2 && y < heightB2+sizebutton
                    if isempty(respB)
                        Screen('FillRect', w ,colB,[x, heightS2-heightslider/2+heightscale/2, x+widthslider, heightS2+heightslider/2+heightscale/2]);
                        Screen('Flip', w, 0, 1);
                        respB=x-leftborderscale;
                    else
                        Screen('FillRect', w ,colBG, [leftborderscale,heightS2-heightimg,leftborderscale+widthscale,heightS2+heightimg]);
                        Screen('DrawTexture', w, texture, [], posimg2);
                        Screen('FillRect', w ,colscale,[leftborderscale, heightS2, leftborderscale+widthscale, heightS2+heightscale]);
                        Screen('FillRect', w ,colB,[x, heightS2-heightslider/2+heightscale/2, x+widthslider, heightS2+heightslider/2+heightscale/2]);
                        Screen('Flip', w, 0, 1);
                        respB=x-leftborderscale;
                    end
                elseif y > heightB3 && y < heightB3+sizebutton
                    if isempty(respR)
                        Screen('FillRect', w ,colR,[x, heightS3-heightslider/2+heightscale/2, x+widthslider, heightS3+heightslider/2+heightscale/2]);
                        Screen('Flip', w, 0, 1);
                        respR=x-leftborderscale;
                    else
                        Screen('FillRect', w ,colBG, [leftborderscale,heightS3-heightimg,leftborderscale+widthscale,heightS3+heightimg]);
                        Screen('DrawTexture', w, texture, [], posimg3);
                        Screen('FillRect', w ,colscale,[leftborderscale, heightS3, leftborderscale+widthscale, heightS3+heightscale]);
                        Screen('FillRect', w ,colR,[x, heightS3-heightslider/2+heightscale/2, x+widthslider, heightS3+heightslider/2+heightscale/2]);
                        Screen('Flip', w, 0, 1);
                        respR=x-leftborderscale;
                    end
                elseif y > heightB4 && y < heightB4+sizebutton
                    if isempty(respY)
                        Screen('FillRect', w ,colY,[x, heightS4-heightslider/2+heightscale/2, x+widthslider, heightS4+heightslider/2+heightscale/2]);
                        Screen('Flip', w, 0, 1);
                        respY=x-leftborderscale;
                    else
                        Screen('FillRect', w ,colBG, [leftborderscale,heightS4-heightimg,leftborderscale+widthscale,heightS4+heightimg]);
                        Screen('DrawTexture', w, texture, [], posimg4);
                        Screen('FillRect', w ,colscale,[leftborderscale, heightS4, leftborderscale+widthscale, heightS4+heightscale]);
                        Screen('FillRect', w ,colY,[x, heightS4-heightslider/2+heightscale/2, x+widthslider, heightS4+heightslider/2+heightscale/2]);
                        Screen('Flip', w, 0, 1);
                        respY=x-leftborderscale;
                    end
                end
            elseif x > leftbordersound && x < leftbordersound+sizebutton
                if y > heightB1 && y < heightB1+sizebutton
                    PsychPortAudio('FillBuffer', pahandle, soundG');
                    PsychPortAudio('Start', pahandle, 1, 0, 1);
                    WaitSecs(0.5)
                elseif y > heightB2 && y < heightB2+sizebutton
                    PsychPortAudio('FillBuffer', pahandle, soundB');
                    PsychPortAudio('Start', pahandle, 1, 0, 1);
                    WaitSecs(0.5)
                elseif y > heightB3 && y < heightB3+sizebutton
                    PsychPortAudio('FillBuffer', pahandle, soundR');
                    PsychPortAudio('Start', pahandle, 1, 0, 1);
                    WaitSecs(0.5)
                elseif y > heightB4 && y < heightB4+sizebutton
                    PsychPortAudio('FillBuffer', pahandle, soundY');
                    PsychPortAudio('Start', pahandle, 1, 0, 1);
                    WaitSecs(0.5)
                end
            end
        end
    end

    resp(uu,1:6) = [sourceIntConds(pp) f0conds(pp) f01 Azimuths(azi1) distG respG];
    resp(uu+1,1:6)=[sourceIntConds(pp) f0conds(pp) f02 Azimuths(azi2) distB respB];
    resp(uu+2,1:6)=[sourceIntConds(pp) f0conds(pp) f03 Azimuths(azi3) distR respR];
    resp(uu+3,1:6)=[sourceIntConds(pp) f0conds(pp) f04 Azimuths(azi4) distY respY];

    save(['ext_resp_' num2str(subNum)],'resp');
    uu=uu+4; % trial counter
    pause(1)
end % page

results.participant=subNum;
results.cols={'sourceIntensity','f0condition','f0','azi','distance','externalization'};
results.res=resp;
save(['ext_resp_' num2str(subNum)],'results');

% End text for participant
Screen('FillRect', w ,colBG, [0 0 screenXpix screenYpix]);
DrawFormattedText(w,'Vielen Dank! Das Experiment ist abgeschlossen.',.2*x_center,'center',[0 0 0],120,0,0,1.5);
Screen('Flip',w);
WaitSecs(3);

PsychPortAudio('Close', pahandle);
sca;
end % function

%% Blocks loop

% start from the block specified by the parameters/settings parts
%     for block = startBlockNo:noBlocks
%
%         % If low source intensity block, scale down
%         if sourceInt(block) == 0
%             SPL = db(0.2);
%             disp('Low intensity block');
%         else
%             SPL = db(1);
%             disp('High intensity block');
%         end
%
%         % uniform background
%         Screen('FillRect', win, backGroundColor);
%         Screen('Flip', win);
%
%         % wait for releasing keys before going on
%         releaseStart = GetSecs;
%         KbReleaseWait([], releaseStart+2);
%
%         % fill a dynamic buffer with data for whole block
%         % get trial index list for current block
%         trialList = trialIdx(blockIdx==block);
%         buffer = [];
%         for trial = min(trialList):max(trialList)
%             audioData = stimStruct{trial, audioColumn};
%             % Intensity roving:
%             audioData = 10^(SPL/20)*audioData;
%             buffer(end+1) = PsychPortAudio('CreateBuffer', [], audioData');
%         end
%
%         % counter for trials in given block
%         trialCounterForBlock = 0;
%
%         % user message
%         disp([char(10), 'Buffered all stimuli for block ', num2str(block),...
%             ', showing block start message']);
%
%         % block starting text
%         if sourceInt(block) == 0
%             intensityInfo = 'This block will have relatively soft sounds.';
%         else
%             intensityInfo = 'This block will have relatively loud sounds.';
%         end
%
%         blockStartText = ['Beginning block ', num2str(block), '. \n\n',...
%             intensityInfo, '\n',...
%             'There will be ', num2str(length(trialList)), ' trials in the block.\n\n\n',...
%             'Ready? Press SPACE to begin.'];
%
%         % uniform background
%         Screen('FillRect', win, backGroundColor);
%         % draw block-starting text
%         DrawFormattedText(win, blockStartText, 'center', 'center', textColor);
%         Screen('Flip', win);
%         % wait for key press to start
%         while 1
%             [keyIsDown, ~, keyCode] = KbCheck;
%             if keyIsDown
%                 % if subject is ready to start
%                 if find(keyCode) == keys.go
%                     break;
%                 elseif find(keyCode) == keys.abort
%                     abortFlag = 1;
%                     break;
%                 end
%             end
%         end
%
%         if abortFlag
%             if flag.triggers
%                 % ppdev_mex('Close', 1);
%                 IOport('Close',TB);
%             end
%             ListenChar(0);
%             Priority(0);
%             RestrictKeysForKbCheck([]);
%             PsychPortAudio('Close');
%             Screen('CloseAll');
%             ShowCursor(screenNumber);
%             diary off
%             return;
%         end
%
%         %% Trials loop
%         % trial loop (over the trials for given block)
%         for trial = min(trialList):max(trialList)
%
%             % relative trial number (trial in given block)
%             trialCounterForBlock = trialCounterForBlock+1;
%
%             % background with fixation cross, get trial start timestamp
%             Screen('CopyWindow', fixCrossWin, win);
%             Screen('DrawingFinished', win);
%             trialStart = Screen('Flip', win);
%
%             % user message
%             disp([char(10), 'Starting trial ', num2str(trialCounterForBlock)]);
%
%             % fill audio buffer with next stimuli
%             PsychPortAudio('FillBuffer', pahandle, buffer(trialCounterForBlock));
%
%             % wait till we are 100 ms from the start of the playback
%             while GetSecs-trialStart <= iti-100
%                 WaitSecs(0.001);
%             end
%
%             % blocking playback start for precision
%             startTime = PsychPortAudio('Start', pahandle, 1, trialStart+iti, 1);
%
%             % user message
%             disp(['Audio started at ', num2str(startTime-trialStart), ' secs after trial start']);
%             disp(['(Target ITI was ', num2str(iti), ' secs)']);
%
%             % Wait for sound to end, then display 3AFC prompt
%             WaitSecs(totalDur(trial));
%             promptText = ['How far was the sound? \n\n',...
%                 '1 - in my head \n 2 - close to my ear \n 3 - farther from my ear'];
%             Screen('FillRect', win, backGroundColor);
%             DrawFormattedText(win, promptText, 'center', 'center', textColor);
%             Screen('Flip', win);
%
%             % Wait for response (target detection)
%             respFlag = 0;
%             while GetSecs <= startTime+totalDur(trial)+respInt % response interval hasn't ended
%                 [keyIsDown, ~, keyCode] = KbCheck;
%                 if keyIsDown % if subject responded
%                     if ismember(find(keyCode),respKeys) % and it was a response button
%                         respFlag = 1;
%                         response(trial) = str2double(KbName(find(keyCode)));
%                         if find(keyCode) == keys.intern % if response was in-the-head
%                             if cell2mat(stimStruct(trial,distanceColumn)) == 0 % and distance is 0
%                                 accuracy(trial) = 1;
%                             else
%                                 accuracy(trial) = 0;
%                             end
%                         elseif find(keyCode) == keys.near % if response was 20 cm
%                             if cell2mat(stimStruct(trial,distanceColumn)) == 0.2 % and distance is 20 cm
%                                 accuracy(trial) = 1;
%                             else
%                                 accuracy(trial) = 0;
%                             end
%                         elseif find(keyCode) == keys.far % if response is 1 m
%                             if cell2mat(stimStruct(trial,distanceColumn)) == 1 % and distance is 1 m
%                                 accuracy(trial) = 1;
%                             else
%                                 accuracy(trial) = 0;
%                             end
%                         end
%                         break;
%                     elseif find(keyCode) == keys.abort % if it was the escape button
%                         abortFlag = 1;
%                         break;
%                     end
%                 end
%             end
%             % if abort was requested, quit
%             if abortFlag
%                 ListenChar(0);
%                 Priority(0);
%                 RestrictKeysForKbCheck([]);
%                 PsychPortAudio('Close');
%                 Screen('CloseAll');
%                 ShowCursor(screenNumber);
%                 diary off
%                 return;
%             end
%
%             % user messages
%             if isnan(accuracy(trial))
%                 disp('Subject did not respond in time');
%             else
%                 disp(['Distance: ', num2str(cell2mat(stimStruct(trial,distanceColumn))),' m']);
%                 disp(['Response: ', num2str(response(trial))]);
%                 disp(['Accuracy: ', num2str(accuracy(trial))]);
%             end
%
%             % Cumulative accuracy and RT in block
%             blockAcc = sum(accuracy(trial-trialCounterForBlock+1:trial), 'omitnan')/trialCounterForBlock*100;
%             disp(['Overall accuracy in block so far is ', num2str(blockAcc), '%']);
%
%             % accumulating results in logging / results variable
%             logVar(trial+1,strcmp(logHeader, 'sourceInt')) = {sourceInt(block)};
%             logVar(trial+1,strcmp(logHeader, 'accuracy')) = {accuracy(trial)};
%             logVar(trial+1,strcmp(logHeader, 'response')) = {response(trial)};
%
%             % save logging/results variable
%             if f0cond == 1
%                 save(subLogF1, 'logVar');
%             elseif f0cond == 2
%                 save(subLogF2, 'logVar');
%             end
%
%         end  % trial for loop
%
%         % user messages
%         disp([char(10), char(10), 'Block no. ', num2str(block), ' has ended,'...
%             'showing block-ending text to participant']);
%         disp([char(10), 'Overall accuracy in block was ', num2str(blockAcc),'%']);
%
%         %% Feedback to subject at the end of block
%         % if not last block
%         if block ~= noBlocks
%             % block ending text
%             blockEndText = ['End of block ', num2str(block), '! \n\n\n',...
%                 'You had ',num2str(blockAcc),'% correct responses in this block. \n',...
%                 'Press SPACE to begin next block.'];
%             % uniform background
%             Screen('FillRect', win, backGroundColor);
%             % draw block-starting text
%             DrawFormattedText(win, blockEndText, 'center', 'center', textColor);
%             Screen('Flip', win);
%             % wait for key press to start
%             while 1
%                 [keyIsDown, ~, keyCode] = KbCheck;
%                 if keyIsDown
%                     % if subject is ready to start
%                     if find(keyCode) == keys.go
%                         break;
%                     elseif find(keyCode) == keys.abort
%                         abortFlag = 1;
%                         break;
%                     end
%                 end
%             end
%             if abortFlag
%                 ListenChar(0);
%                 Priority(0);
%                 RestrictKeysForKbCheck([]);
%                 PsychPortAudio('Close');
%                 Screen('CloseAll');
%                 ShowCursor(screenNumber);
%                 diary off
%                 return;
%             end
%         elseif block == noBlocks % if this was the last block
%             % user message
%             disp([char(10), 'The task has ended.']);
%             blockEndText = ['End of last block. \n',...
%                 'You had ',num2str(blockAcc),'% correct responses in this block. \n',...
%                 'We can take a break...'];
%             % uniform background
%             Screen('FillRect', win, backGroundColor);
%             % draw block-starting text
%             DrawFormattedText(win, blockEndText, 'center', 'center', textColor);
%             Screen('Flip', win);
%             WaitSecs(5);
%
%         end  % if last block
%
%     end  % block for loop
% end % round for loop
% diary off
%
% %% Ending, cleaning up
% ListenChar(0);
% Priority(0);
% RestrictKeysForKbCheck([]);
% PsychPortAudio('Close');
% Screen('CloseAll');
% ShowCursor;
% return
%
% end % function