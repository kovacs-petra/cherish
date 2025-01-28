function rampAndFilter(trajectory,offsetAzimuth)
% Inputs:
% trajectory    - 1: looming, 2: receding, 3: rotating near, 4: rotating
%                 far
% offsetAzimuth - integer, e.g 90 or -90

%% Load HRTF dataset for spatialization
if not(exist("SOFAdbPath.m","file"))
    sofaPath = '\\kfs\fileserver\ProjektDaten\CherISH\code\SOFAtoolbox\SOFAtoolbox';
    addpath(sofaPath);
    SOFAstart;
end
database = 'scut';
HRTFfilename = 'SCUT_KEMAR_radius_all_interp2.sofa';
fullfn = fullfile(SOFAdbPath, 'database', database, HRTFfilename);
Obj = SOFAload(fullfn);
fs = Obj.Data.SamplingRate;

savename = ['sample',num2str(trajectory),'_',num2str(offsetAzimuth),'_scut.mat'];

PPS = 0.2;
EPS = 2.0;
v = 1;

% frequency = (randi(4,1)+5)*100; % 600-900 Hz with round 100 values
frequency = 700;

    % Set radii
    switch trajectory
        case 1 % Looming
            rOnset = [EPS EPS];  % stationary portion at the onset
            rOffset = [PPS PPS]; % stationary portion at the offset
            rMoving = [EPS PPS]; % moving portion

        case 2 % Receding
            rOnset = [PPS PPS];  % stationary portion at the onset
            rOffset = [EPS EPS]; % stationary portion at the offset
            rMoving = [PPS EPS]; % moving portion

        case 3 % Rotating near
            rOnset = [PPS PPS];  % stationary portion at the onset
            rOffset = rOnset;    % stationary portion at the offset
            rMoving = rOnset;    % moving portion

        case 4 % Rotating far
            rOnset = [EPS EPS];  % stationary portion at the onset
            rOffset = rOnset;    % stationary portion at the offset
            rMoving = rOnset;    % moving portion
    end

    %% Set durations
    % Stationary portions: 600-800 ms with round 100 values, in s
    durStatOnset = ((randi(4,1)+4)/10); 
    durStatOffset = ((randi(4,1)+4)/10);

    % Moving portion: duration calculated from distance and velocity
    durMoving = abs(EPS-PPS)/v;

    % Calculate total duration WARNING should be done later prolly
    totalDur = durStatOnset+durMoving+durStatOffset; % in s

    ele = linspace(0,0,round(totalDur*fs));

    % Set trajectory vector based on the durations
    r = [linspace(rOnset(1), rOnset(2), round(durStatOnset*fs)), ...
            linspace(rMoving(1), rMoving(2), round(durMoving*fs)), ...
            linspace(rOffset(1), rOffset(2), round(durStatOffset*fs))];

    %% Set azimuth
    if trajectory <= 2 % looming or receding
        % aziMain = [offsetAzimuth offsetAzimuth];
        azi = [linspace(offsetAzimuth, offsetAzimuth, round(durStatOnset*fs)), ...
            linspace(offsetAzimuth, offsetAzimuth, round(durMoving*fs)), ...
            linspace(offsetAzimuth, offsetAzimuth, round(durStatOffset*fs))];
    elseif trajectory >= 3 % rotating near or far
        % aziMain = [0 offsetAzimuth];
        azi = [linspace(0, 0, round(durStatOnset*fs)), ...
            linspace(0, offsetAzimuth, round(durMoving*fs)), ...
            linspace(offsetAzimuth, offsetAzimuth, round(durStatOffset*fs))];        
    end

    %% Generate triangle waves for each portion
    % Generate square waves for stationary portions
    tStatOnset = linspace(0,durStatOnset,durStatOnset*fs); % Time vector from 0 to dur with sampling interval 1/fs
    statOnset = sawtooth(2 * pi * frequency * tStatOnset, 0.5); % 0.5 for a symmetric triangular wave
    tStatOffset = linspace(0,durStatOffset,durStatOffset*fs);
    statOffset  = sawtooth(2 * pi * frequency * tStatOffset, 0.5);

    % Generate square waves for moving portion
    tMoving = linspace(0,durMoving,durMoving*fs);
    moving = sawtooth(2 * pi * frequency * tMoving, 0.5);

    %% Generate targets
    targetITI = 0.1; % 100 ms
    target = sin(2*pi*2400*(0:(1/fs):targetITI-1/fs)); % 100 ms, 2400 Hz
    gap = zeros(targetITI*fs,1);

    % Ramp the target
    rampSamples = round(fs*0.01);
    rampOnset = sin(linspace(0,pi/2,rampSamples));
    rampOffset = cos(linspace(0,pi/2,rampSamples));
    target(1:rampSamples) = target(1:rampSamples)'.*rampOnset';
    target((end-rampSamples+1):end) = target((end-rampSamples+1):end)'.*rampOffset';

    %% Concatenate the sound portions and spatialize the stimuli
    C = [statOnset moving statOffset]; % cue
    win = tukeywin(size(C,2),(0.1*fs)/size(C,2)); % 0.1 s on-offset ramp
    C = C.*win';

    % if targetTrials(stimNo) == 1
        T = [gap' target; gap' target]; % target
    % else
    %     T = [gap' gap'; gap' gap'];
    % end

    [cue, aziA, eleA, rA, ~] = local_SOFAspat(C',Obj,azi,ele,r);
    if trajectory < 4
        dbChange = 20 * log10(1 / max(cue,[],"all")); % scale to [-1,1]
 
    else
        dbChange = 20 * log10(0.1 / max(cue,[],"all")); % scale to [-0.1,0.1] because it's far
    end
    cue = db2mag(dbChange)*cue;

    out = [cue' T];
    save(savename, "out", "fs");

end

function [out, aziActual, eleActual, rActual, idx] = local_SOFAspat(signal,Obj,azi,ele,r)
if length(r)~=length(signal)
    errorText = ['Signal (length: ', num2str(length(signal)), ') and r (length: ',...
        num2str(length(r)),') need to have the same length!'];
    error(errorText)
end
[out, aziActual, eleActual, rActual, idx] = SOFAspat(signal./r(:),Obj,azi,ele,r);
end