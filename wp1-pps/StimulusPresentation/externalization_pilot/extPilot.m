function extPilot(subNum)
% Usage: extPilot(subNum)
%
% Externalization pilot for CherISH DC6 WP1
%
% Input:
%   - subNum:   subject number, integer
%
% - Distance judgement task with a GUI
% - Stimuli are presented on pages, every page contains 4 sounds that can be
% played by the listener. Each sound on a page has a different distance
% (options: 20 cm, 40/90/130/180 cm, 2 m)
% - Source intensity changes at halftime from low to high, with a break
% inbetween
% - Two stimulus types are tested separately: triangle wave, white noise
% - Stimuli have to be arranged in an f0s X distances X azimuths X stimulus types cell
% array, with each cell containing the audio (2 channels) in position 1 and
% the specific distance in position 2
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
% Author: Petra Kovacs, 2025

%% Set trial number info
nTrials = 192; % how many trials in total
nTrialsPerSI = nTrials/2; % how many trials per source intensity condition
nTrialsPerType = nTrialsPerSI/2; % how many trials per stimulus type
nTrialsPerPage = 4; % how many trials on a page

%% Initialize PsychPortAudio and Screen Number
%%% Thanks to David Meijer for providing code

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

if strcmp(getenv('computername'),'EXPGRUEN')  % Manually overwrite the default sampling rate on the lab computer (default is 44100, but that leads to auditory delays relative to visual)
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

fs = 48000;
pahandle = PsychPortAudio('Open', S.sound_device_idx, 1, 1, fs, 2);
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1);
% InitializePsychSound(1);

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

% disp(pahandleStatus);

% initial start & stop of audio device to avoid potential initial latencies
tmpSound = zeros(2, fs/10);  % silence
tmpBuffer = PsychPortAudio('CreateBuffer', [], tmpSound);  % create buffer
PsychPortAudio('FillBuffer', pahandle, tmpBuffer);  % fill the buffer of audio device with silence
PsychPortAudio('Start', pahandle, 1);  % start immediately
PsychPortAudio('Stop', pahandle, 1);  % stop when playback is over

%% Open the stimuli structure
stimArrayName = 'stimStructEq.mat';
stimArrayFileStruct = dir(stimArrayName);
stimArrayFile = [stimArrayFileStruct.folder, '/', stimArrayFileStruct.name];
load(stimArrayFile, 'stimStructEq');

%% Initialize the graphical interface for Psychtoolbox
colBG=[200 200 200];

% ignore screen warning
[w,rect] = Screen('OpenWindow',S.screen_number, colBG);
% debugRect = [10 10 400 1200];
%     [w,rect] = Screen('OpenWindow',S.screen_number,colBG,debugRect);
%     ShowCursor;

[screenXpix,screenYpix] = Screen('WindowSize', w);
[x_center,y_center] = RectCenter(rect);

Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%% Listener instruction
instruction = [...
    'Wie weit ist die Fliege entfernt? \n\n\n',...
    'Pro Seite gibt es vier Töne, die Sie durch Mausklick auf die Lautsprechersymbole links im Bild so oft Sie wollen anhören können. \n',...
    'Stellen Sie sich vor, diese Töne stammen von verschiedenen Fliegen – allerdings sind es außerirdische Fliegen mit seltsamen Stimmen. \n',...
    'Rechts neben den Lautsprechersymbolen gibt es pro Ton eine Leiste, auf der Sie den Abstand der Fliege per Mausklick stufenlos angeben können. \n',...
    'Die Fliegen können sich innerhalb oder außerhalb Ihrer Reichweite befinden, entweder links oder rechts. Ihre Stimme – also die Tonhöhe – kann sich bei jedem Mausklick verändern. \n\n',...
    'Die Skala der Leiste reicht bis etwa zwei Meter Entfernung. \n',...
    'Wenn Sie die Abstände der vier Fliegen ausreichend verglichen haben und mit Ihrer Einschätzung zufrieden sind, drücken Sie die Leertaste, um fortzufahren. \n',...
    'Die Hälfte der Seiten enthält leise Fliegen, die andere Hälfte laute. Dazwischen machen wir eine Pause. \n', ...
    'Die relative Lautstärke der Fliegen auf derselben Seite ist dabei informativ. \n',...
    'Das Experiment besteht insgesamt aus ' num2str(nTrials/nTrialsPerPage) ' Seiten. \n\n\n',...
    'Haben Sie noch Fragen? \n\n',...
    'Zum Starten die Leertaste drücken.'];

DrawFormattedText(w,instruction,.2*x_center,'center',[0 0 0],120,0,0,1.5);
Screen('Flip',w);

%% define colors, sizes and positions
% img of head
path_head = "C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\head.png";
% path_head = "C:\Users\experimentator.KFS\Documents\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\head_trialversion.png";
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
path_symb = "C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\soundsymbol.png";
% path_symb = "C:\Users\experimentator.KFS\Documents\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\soundsymbol.png";
[imgSymb, ~, alphasymb] = imread(path_symb);
imgSymb(:, :, 4) = alphasymb;
textureS = Screen('MakeTexture', w, imgSymb);

% fly img
path_fly = "C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\fly.png";
[imgFly, ~, alphaFly] = imread(path_fly);
imgFly(:, :, 4) = alphaFly;
textureF = Screen('MakeTexture', w, imgFly);
heightFly = 40;
widthFly = 30;
gapFly = 5;


%% experiment
% wait for key after instructions
while ~KbCheck; end
WaitSecs(0.6);

%% Azimuths
Azimuths = [90,-90];

%% Set order of conditions
% Source intensity conditions
lowInt = zeros(1,nTrialsPerSI/nTrialsPerPage); highInt = ones(1,nTrialsPerSI/nTrialsPerPage);
sourceIntConds = [lowInt highInt];
% if mod(subNum,2) == 0 
%     sourceIntConds = flip(sourceIntConds); % counterbalance the order
% end

% Stimulus type conditions
TW = zeros(1,nTrialsPerType/nTrialsPerPage);
WN = ones(1,nTrialsPerType/nTrialsPerPage);
typeConds = [TW,WN,TW,WN];
if mod(subNum,2) == 0 
    typeConds = flip(typeConds); % counterbalance the order
end

%% Counters
resp=[];
uu=1; % trial counter

for pp = 1:(nTrials/nTrialsPerPage) % Page

    % Handle source intensity condition
    if sourceIntConds(pp) == 0 % low intensity condition
        SPL = db(0.2);
    else % high intensity condition
        SPL = db(2);
    end

    % Handle stimulus type condition
    if typeConds(pp) == 0 % triangle wave condition
        tt = 1; % stim type index in the stimStruct
    else % white noise condition
        tt = 2;
    end

    % azi chosen randomly for each trial in a page
    azi1 = randi(length(Azimuths),1);
    azi2 = randi(length(Azimuths),1);
    azi3 = randi(length(Azimuths),1);
    azi4 = randi(length(Azimuths),1);

    % distance counterbalanced within a page
    distFix = 1:2; % 0.2 and 2 m all have to be presented;
    % they are found in columns 1:2 of
    % the stim cell
    distRand = randperm(4,2)+2; % One instance of the distances 0.3-1.9 m has to be present
    distances = [distRand,distFix];
    distances = distances(randperm(length(distances)));

    % Load sound
    cellG = stimStructEq.stim(:,distances(1),azi1,tt);
    cellB = stimStructEq.stim(:,distances(2),azi2,tt);
    cellR = stimStructEq.stim(:,distances(3),azi3,tt);
    cellY = stimStructEq.stim(:,distances(4),azi4,tt);

    distG = cellG{1,1}{1,2};
    distB = cellB{1,1}{1,2};
    distR = cellR{1,1}{1,2};
    distY = cellY{1,1}{1,2};

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
            intensityInfo = 'Die nächsten Seiten haben LEISE Töne...';
            disp('Low intensity pages.')
        else
            intensityInfo = 'Die nächsten Seiten haben LAUTE Töne...';
            disp('High intensity pages.')
        end

    % display info
    Screen('TextSize',w,30);
    DrawFormattedText(w,intensityInfo,'center','center',[0 0 0],120,0,0,1.5);
    Screen('Flip', w, 0, 1);
    WaitSecs(3);
    Screen('TextSize',w,defaultTxtSize);
    end  

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
                        sliderTopY = heightS1 - heightslider/2 + heightscale/2;
                        flyLeft = x - widthFly/2;
                        flyTop = sliderTopY - gapFly - heightFly;
                        flyRight = x + widthFly/2;
                        flyBottom = sliderTopY - gapFly;
                        
                        % Draw the fly texture
                        Screen('DrawTexture', w, textureF, [], [flyLeft, flyTop, flyRight, flyBottom]);
                        Screen('Flip', w, 0, 1);
                        respG=x-leftborderscale;
                    else
                        Screen('FillRect', w ,colBG, [leftborderscale,heightS1-heightimg,leftborderscale+widthscale,heightS1+heightimg]);
                        Screen('DrawTexture', w, texture, [], posimg1);
                        Screen('FillRect', w ,colscale,[leftborderscale, heightS1, leftborderscale+widthscale, heightS1+heightscale]);
                        Screen('FillRect', w ,colG,[x, heightS1-heightslider/2+heightscale/2, x+widthslider, heightS1+heightslider/2+heightscale/2]);
                        sliderTopY = heightS1 - heightslider/2 + heightscale/2;
                        flyLeft = x - widthFly/2;
                        flyTop = sliderTopY - gapFly - heightFly;
                        flyRight = x + widthFly/2;
                        flyBottom = sliderTopY - gapFly;
                        Screen('DrawTexture',w,textureF,[],[flyLeft, flyTop, flyRight, flyBottom]);
                        Screen('Flip', w, 0, 1);
                        respG=x-leftborderscale;
                    end
                elseif y > heightB2 && y < heightB2+sizebutton
                    if isempty(respB)
                        Screen('FillRect', w ,colB,[x, heightS2-heightslider/2+heightscale/2, x+widthslider, heightS2+heightslider/2+heightscale/2]);
                        sliderTopY = heightS2 - heightslider/2 + heightscale/2;
                        flyLeft = x - widthFly/2;
                        flyTop = sliderTopY - gapFly - heightFly;
                        flyRight = x + widthFly/2;
                        flyBottom = sliderTopY - gapFly;                        
                        % Draw the fly texture
                        Screen('DrawTexture', w, textureF, [], [flyLeft, flyTop, flyRight, flyBottom]);
                        Screen('Flip', w, 0, 1);
                        respB=x-leftborderscale;
                    else
                        Screen('FillRect', w ,colBG, [leftborderscale,heightS2-heightimg,leftborderscale+widthscale,heightS2+heightimg]);
                        Screen('DrawTexture', w, texture, [], posimg2);
                        Screen('FillRect', w ,colscale,[leftborderscale, heightS2, leftborderscale+widthscale, heightS2+heightscale]);
                        Screen('FillRect', w ,colB,[x, heightS2-heightslider/2+heightscale/2, x+widthslider, heightS2+heightslider/2+heightscale/2]);
                        sliderTopY = heightS2 - heightslider/2 + heightscale/2;
                        flyLeft = x - widthFly/2;
                        flyTop = sliderTopY - gapFly - heightFly;
                        flyRight = x + widthFly/2;
                        flyBottom = sliderTopY - gapFly;                        
                        % Draw the fly texture
                        Screen('DrawTexture', w, textureF, [], [flyLeft, flyTop, flyRight, flyBottom]);
                        Screen('Flip', w, 0, 1);
                        respB=x-leftborderscale;
                    end
                elseif y > heightB3 && y < heightB3+sizebutton
                    if isempty(respR)
                        Screen('FillRect', w ,colR,[x, heightS3-heightslider/2+heightscale/2, x+widthslider, heightS3+heightslider/2+heightscale/2]);
                        sliderTopY = heightS3 - heightslider/2 + heightscale/2;
                        flyLeft = x - widthFly/2;
                        flyTop = sliderTopY - gapFly - heightFly;
                        flyRight = x + widthFly/2;
                        flyBottom = sliderTopY - gapFly;                        
                        % Draw the fly texture
                        Screen('DrawTexture', w, textureF, [], [flyLeft, flyTop, flyRight, flyBottom]);
                        Screen('Flip', w, 0, 1);
                        respR=x-leftborderscale;
                    else
                        Screen('FillRect', w ,colBG, [leftborderscale,heightS3-heightimg,leftborderscale+widthscale,heightS3+heightimg]);
                        Screen('DrawTexture', w, texture, [], posimg3);
                        Screen('FillRect', w ,colscale,[leftborderscale, heightS3, leftborderscale+widthscale, heightS3+heightscale]);
                        Screen('FillRect', w ,colR,[x, heightS3-heightslider/2+heightscale/2, x+widthslider, heightS3+heightslider/2+heightscale/2]);
                        sliderTopY = heightS3 - heightslider/2 + heightscale/2;
                        flyLeft = x - widthFly/2;
                        flyTop = sliderTopY - gapFly - heightFly;
                        flyRight = x + widthFly/2;
                        flyBottom = sliderTopY - gapFly;                        
                        % Draw the fly texture
                        Screen('DrawTexture', w, textureF, [], [flyLeft, flyTop, flyRight, flyBottom]);
                        Screen('Flip', w, 0, 1);
                        respR=x-leftborderscale;
                    end
                elseif y > heightB4 && y < heightB4+sizebutton
                    if isempty(respY)
                        Screen('FillRect', w ,colY,[x, heightS4-heightslider/2+heightscale/2, x+widthslider, heightS4+heightslider/2+heightscale/2]);
                        sliderTopY = heightS4 - heightslider/2 + heightscale/2;
                        flyLeft = x - widthFly/2;
                        flyTop = sliderTopY - gapFly - heightFly;
                        flyRight = x + widthFly/2;
                        flyBottom = sliderTopY - gapFly;                        
                        % Draw the fly texture
                        Screen('DrawTexture', w, textureF, [], [flyLeft, flyTop, flyRight, flyBottom]);
                        Screen('Flip', w, 0, 1);
                        respY=x-leftborderscale;
                    else
                        Screen('FillRect', w ,colBG, [leftborderscale,heightS4-heightimg,leftborderscale+widthscale,heightS4+heightimg]);
                        Screen('DrawTexture', w, texture, [], posimg4);
                        Screen('FillRect', w ,colscale,[leftborderscale, heightS4, leftborderscale+widthscale, heightS4+heightscale]);
                        Screen('FillRect', w ,colY,[x, heightS4-heightslider/2+heightscale/2, x+widthslider, heightS4+heightslider/2+heightscale/2]);
                        sliderTopY = heightS4 - heightslider/2 + heightscale/2;
                        flyLeft = x - widthFly/2;
                        flyTop = sliderTopY - gapFly - heightFly;
                        flyRight = x + widthFly/2;
                        flyBottom = sliderTopY - gapFly;                        
                        % Draw the fly texture
                        Screen('DrawTexture', w, textureF, [], [flyLeft, flyTop, flyRight, flyBottom]);
                        Screen('Flip', w, 0, 1);
                        respY=x-leftborderscale;
                    end
                end
            elseif x > leftbordersound && x < leftbordersound+sizebutton
                if y > heightB1 && y < heightB1+sizebutton
                    soundG = cellG{randi([1,length(cellG)],1),1};
                    soundG = soundG{1,1}*10^(SPL/20);
                    PsychPortAudio('FillBuffer', pahandle, soundG');
                    PsychPortAudio('Start', pahandle, 1, 0, 1);
                    WaitSecs(0.5);
                elseif y > heightB2 && y < heightB2+sizebutton
                    soundB = cellB{randi([1,length(cellB)],1),1};
                    soundB = soundB{1,1}*10^(SPL/20);
                    PsychPortAudio('FillBuffer', pahandle, soundB');
                    PsychPortAudio('Start', pahandle, 1, 0, 1);
                    WaitSecs(0.5);
                elseif y > heightB3 && y < heightB3+sizebutton
                    soundR = cellR{randi([1,length(cellR)],1),1};
                    soundR = soundR{1,1}*10^(SPL/20);
                    PsychPortAudio('FillBuffer', pahandle, soundR');
                    PsychPortAudio('Start', pahandle, 1, 0, 1);
                    WaitSecs(0.5);
                elseif y > heightB4 && y < heightB4+sizebutton
                    soundY = cellY{randi([1,length(cellY)],1),1};
                    soundY = soundY{1,1}*10^(SPL/20);
                    PsychPortAudio('FillBuffer', pahandle, soundY');
                    PsychPortAudio('Start', pahandle, 1, 0, 1);
                    WaitSecs(0.5);
                end
            end
        end
    end

    resp(uu,1:5) = [sourceIntConds(pp) tt Azimuths(azi1) distG respG];
    resp(uu+1,1:5)=[sourceIntConds(pp) tt Azimuths(azi2) distB respB];
    resp(uu+2,1:5)=[sourceIntConds(pp) tt Azimuths(azi3) distR respR];
    resp(uu+3,1:5)=[sourceIntConds(pp) tt Azimuths(azi4) distY respY];

    save(['ext_resp_' num2str(subNum)],'resp');
    uu=uu+4; % trial counter
    WaitSecs(1);

    if pp == nTrialsPerSI/nTrialsPerPage % intensity change
        Screen('FillRect', w ,colBG, [0 0 screenXpix screenYpix]);
        Screen('Flip', w, 0, 1);
        % Break message and weiter mit leertaste
        Screen('TextSize',w,30);
        DrawFormattedText(w,'Wir machen jetzt ein Paar Minuten Pause.','center','center',[0 0 0],120,0,0,1.5);
        Screen('Flip', w, 0, 1);
        Screen('TextSize',w,defaultTxtSize);
        % wait for key after msg
        while ~KbCheck; end
        WaitSecs(0.6);
    end % intensity change

end % page

results.participant=subNum;
results.cols={'sourceIntensity','stimType','azi','distance','externalization'};
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
