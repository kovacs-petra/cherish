% Example script to replicate BRT bug

%% Generate a sound input for BRT
fs = 48e3;

% Distances used, m
PPS = 0.2;
EPS = 1;

% F0
fmin = 310;
fmax = 2*fmin; % 1 octave
f0_options = round(linspace(fmin,fmax,10),-1); % 10 unique f0 values

% Duration (s)
dur = 5;
updateRate = 2/1000; % update position every 2 ms during recording

% Spatial position (trajectory)
azi = linspace(-90,90,round(dur/updateRate)); % right to left
ele = linspace(0,0,round(dur/updateRate));
r = linspace(1,1,round(dur/updateRate));

% Set connection between Matlab and BRT through OSC; set some settings
u = pnet('udpsocket',10017); % Listen port in BRT
oscsend(u, '/control/connect', 'si', 'localhost',10019);
oscsend(u, '/listener/enableSpatialization', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableInterpolation', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableNearFieldEffect', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableITD', 'sB', 'DefaultListener',1);
oscsend(u, '/environment/enableModel', 'sB', 'FreeField',1);
oscsend(u, '/environment/enableDirectPath', 'sB', 'FreeField',1);
oscsend(u, '/environment/enableReverbPath', 'sB', 'FreeField',1);
oscsend(u, '/listener/enableParallaxCorrection', 'sB', 'DefaultListener',1);
oscsend(u, '/listener/enableModel', 'sB', 'DirectPath', 1);
oscsend(u, '/listener/enableModel', 'sB', 'ReverbPath', 1);
oscsend(u, '/environment/enablePropagationDelay', 'sB', 'FreeField', 1);
oscsend(u, '/environment/enableDistanceAttenuation', 'sB', 'FreeField', 1);

% Generate 10 stimuli
for stimulus = 1:10
    % Progress update
    clc; disp(['Processing stimulus ',num2str(stimulus),'/10...']);

    % Vary the f0
    f0 = f0_options(stimulus);

    % Make triangle wave
    t = linspace(0,dur,dur*fs);
    sig = sawtooth(2 * pi * f0 * tMoving, 0.5);

    % Save as wav for BRT input
    wavname = fullfile(pwd,strcat('sample',num2str(stimulus),'.wav'));
    audiowrite(wavname, sig, fs); 
    matname = fullfile(pwd,strcat('sample',num2str(stimulus),'.mat'));

    % Spatialize with BRT
    [sigSpat, sigParams] = BRTspat(wavname,dur,azi,ele,r,matname,updateRate,u);

    % Add on-offset ramp to the spatialized stimulus
    win = tukeywin(size(sigSpat,2),(0.1*fs)/size(sigSpat,2)); % 0.1 s 
    sigSpatWin = sigSpat.*repmat(win',2,1);

    % Scale sound (PPS sound would be scaled to 1, EPS is scaled to be
    % lower in intensity)
    scaleto = PPS/EPS;
    scaledb = 20 * log10(scaleto / max(sigSpatWin,[],"all"));
    scalesig1 = scaletodbspl(sigSpatWin(1,:)',dbspl(sigSpatWin(1,:)')+scaledb);
    scalesig2 = scaletodbspl(sigSpatWin(2,:)',dbspl(sigSpatWin(2,:)')+scaledb);
    scalesig = [scalesig1 scalesig2];




end


