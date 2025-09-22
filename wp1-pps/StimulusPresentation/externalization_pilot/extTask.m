function extTask(subNum,session)

% Externalization task for CherISH WP1
% usage: extTask(subNum,session)
% subNum - subject number
% session - 1 or 2 

%% Set volume
system('\\KFS\Fileserver\Users\pkovacs\SetVol40.bat.lnk');

%% Open the stimuli structure
stimStructName = 'stimStructExt.mat';
stimFileStruct = dir(stimStructName);
stimStructFile = [stimFileStruct.folder, '/', stimFileStruct.name];
load(stimStructFile, 'stimStructExt');

nTrials = 80; 
nTrialsPerPage = 4;

if mod(subNum,2) == 0 % even subject number
    if session == 1 % 1st session high intensity
        SPL = db(1);
    elseif session == 2 % 2nd session low intensity
        SPL = db(0.2/2);
    end
else % odd subject number
    if session == 1 % 1st session low intensity
        SPL = db(0.2/2);
    else % 2nd session high intensity
        SPL = db(1);
    end
end

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
% % % % % pahandle = PsychPortAudio('Open', S.sound_device_idx, 1, 1, fs, 2);
pahandle = PsychPortAudio('Open', S.sound_device_idx, 1, 1, [], 2);
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1);
% InitializePsychSound(1);

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
    'Stellen Sie sich vor, diese Töne stammen von verschiedenen Fliegen. \n',...
    'Rechts neben den Lautsprechersymbolen gibt es pro Ton eine Leiste, auf der Sie den Abstand der Fliege per Mausklick stufenlos angeben können. \n',...
    'Die Fliegen können sich innerhalb oder außerhalb Ihrer Reichweite befinden, entweder links oder rechts. Ihre Stimme – also die Tonhöhe – kann sich bei jedem Mausklick verändern. \n\n',...
    'Die Skala der Leiste reicht bis etwa zwei Meter Entfernung. \n',...
    'Wenn Sie die Abstände der vier Fliegen ausreichend verglichen haben und mit Ihrer Einschätzung zufrieden sind, drücken Sie die Leertaste, um fortzufahren. \n',...
    'Das Experiment besteht insgesamt aus ' num2str(nTrials/nTrialsPerPage) ' Seiten. \n\n\n',...
    'Haben Sie noch Fragen? \n\n',...
    'Zum Starten die Leertaste drücken.'];

DrawFormattedText(w,instruction,.2*x_center,'center',[0 0 0],120,0,0,1.5);
Screen('Flip',w);

%% define colors, sizes and positions
% img of head
% path_head = "C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\head.png";
path_head = "\\kfs\fileserver\ProjektDaten\CherISH\code\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\head.png";
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
% path_symb = "C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\soundsymbol.png";
path_symb = "\\kfs\fileserver\ProjektDaten\CherISH\code\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\soundsymbol.png";
[imgSymb, ~, alphasymb] = imread(path_symb);
imgSymb(:, :, 4) = alphasymb;
textureS = Screen('MakeTexture', w, imgSymb);

% fly img
% path_fly = "C:\Users\pkovacs\Documents\GitHub\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\fly.png";
path_fly = "\\kfs\fileserver\ProjektDaten\CherISH\code\cherish\wp1-pps\StimulusPresentation\externalization_pilot\img\fly.png";
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

%% Counters
resp=[];
uu=1; % trial counter
respType = zeros(1,nTrials);

for pp = 1:nTrials/nTrialsPerPage
    
    % azi chosen randomly for each trial in a page
    azi1 = randi(length(Azimuths),1);
    azi2 = randi(length(Azimuths),1);
    azi3 = randi(length(Azimuths),1);
    azi4 = randi(length(Azimuths),1);


    % Load sound
    if pp <= 4 % 1st 4 pages have stationary sounds, choose based on distance

        % dummy trajectories variable
        trajectories = ones(1,4);

        % distance counterbalanced within a page
        distFix = 1:2; % 0.2 and 2 m all have to be presented;
        % they are found in columns 1:2 of
        % the stim cell
        distRand = randperm(2,2)+2; % One instance of the distances 0.3-1.9 m has to be present
        distances = [distRand,distFix];
        distances = distances(randperm(length(distances)));

        cellG = stimStructExt.stim(:,distances(1),azi1,1);
        cellB = stimStructExt.stim(:,distances(2),azi2,1);
        cellR = stimStructExt.stim(:,distances(3),azi3,1);
        cellY = stimStructExt.stim(:,distances(4),azi4,1);

        startDistG = cellG{1,1}{1,2}; endDistG = startDistG; 
        startDistB = cellB{1,1}{1,2}; endDistB = startDistB;
        startDistR = cellR{1,1}{1,2}; endDistR = startDistR;
        startDistY = cellY{1,1}{1,2}; endDistY = startDistY;

        respType(uu:uu+3) = 1; % stationary response

    else % rest of the pages have moving sounds, choose based on trajectory

        % trajectory counterbalanced within a page
        trajectories = randperm(4)+1;

        cellG = stimStructExt.stim(:,:,azi1,trajectories(1));
        cellB = stimStructExt.stim(:,:,azi2,trajectories(2));
        cellR = stimStructExt.stim(:,:,azi3,trajectories(3));
        cellY = stimStructExt.stim(:,:,azi4,trajectories(4));

        startDistG = cellG{1,1}{1,2}; 
        switch trajectories(1) 
            case {1,3}; endDistG = 0.2; % looming or PPS
            case {2,4}; endDistG = 2; % receding or EPS            
        end

        startDistB = cellB{1,1}{1,2};
        switch trajectories(2)
            case {1,3}; endDistB = 0.2; % looming or PPS
            case {2,4}; endDistB = 2; % receding or EPS            
        end

        startDistR = cellR{1,1}{1,2};
        switch trajectories(3)
            case {1,3}; endDistR = 0.2; % looming or PPS
            case {2,4}; endDistR = 2; % receding or EPS            
        end

        startDistY = cellY{1,1}{1,2};
        switch trajectories(4)
            case {1,3}; endDistY = 0.2; % looming or PPS
            case {2,4}; endDistY = 2; % receding or EPS            
        end

        if pp <= 12
            respType(uu:uu+3) = 2; % onset response

            % instruction 
            if respType(uu-1) ~= respType(uu)
                instrOnset = [...
                    'Auf den nächsten Seiten können sich die Fliegen bewegen! Bitte schätzen Sie ein, wo jede Fliege BEGINNT sich zu bewegen. \n\n',...
                    'Bereit? Zum Fortfahren die Leertaste drücken.'];

                % clear screen, show instruction
                Screen('FillRect', w ,colBG, [0 0 screenXpix screenYpix]);
                Screen('Flip', w, 0, 1);
                DrawFormattedText(w,instrOnset,.2*x_center,'center',[0 0 0],120,0,0,1.5);
                Screen('Flip',w);

                % wait for key after instructions
                while ~KbCheck; end
                WaitSecs(0.6);
            end

        else
            respType(uu:uu+3) = 3; % offset response

            % instruction
            if respType(uu-1) ~= respType(uu)
                instrOnset = [...
                    'Bitte schätzen Sie auf den nächsten Seiten ein, wo jede Fliege AUFHÖRT sich zu bewegen. \n\n',...
                    'Zum Fortfahren die Leertaste drücken.'];

                % clear screen, show instruction
                Screen('FillRect', w ,colBG, [0 0 screenXpix screenYpix]);
                Screen('Flip', w, 0, 1);
                DrawFormattedText(w,instrOnset,.2*x_center,'center',[0 0 0],120,0,0,1.5);
                Screen('Flip',w);

                % wait for key after instructions
                while ~KbCheck; end
                WaitSecs(0.6);
            end
        end
    end

    % phase to listen to sounds and choose response
    respG=[];
    respB=[];
    respR=[];
    respY=[];

    % clear screen
    Screen('FillRect', w ,colBG, [0 0 screenXpix screenYpix]);
    Screen('Flip', w, 0, 1);

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

    resp(uu,1:7) = [SPL azi1 trajectories(1)-1 startDistG endDistG respG respType(uu)];
    resp(uu+1,1:7)=[SPL azi2 trajectories(2)-1 startDistB endDistB respB respType(uu+1)];
    resp(uu+2,1:7)=[SPL azi3 trajectories(3)-1 startDistR endDistR respR respType(uu+2)];
    resp(uu+3,1:7)=[SPL azi4 trajectories(4)-1 startDistY endDistY respY respType(uu+3)];

    save(['ext_resp_' num2str(subNum),'_',num2str(session)],'resp');
    uu=uu+4; % trial counter
    WaitSecs(1);

end

results.participant=subNum;
results.cols={'SPL','azi','trajectory','startDistance','endDistance','externalization',' respType'};
results.res=resp;
save(['ext_resp_' num2str(subNum),'_',num2str(session)],'results');

PsychPortAudio('Close');
sca;
end % function
