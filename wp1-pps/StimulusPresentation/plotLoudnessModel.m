function plotLoudnessModel(saved_variables)
% Calculates moore2016 loudness model for each stimulus in the CherISH WP1
% externalization task and main task.
% Input: saved_variables - path for saved variables to open

if exist(saved_variables,"file")
    load(saved_variables);
end 

% Add AMT to the path
if ~exist("amt_start.m","file")
    amtPath = '\\kfs.oeaw.ac.at\Fileserver\ProjektDaten\CherISH\code\amtoolbox-full-1.4.0';
    addpath(amtPath);
    amt_start;
end

%% Main task stimuli

% Load main task stimuli
if ~exist("stimArray","var")
    main_stim_path = "\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\code\cherish\wp1-pps\StimulusPresentation\stimArray.mat";
    tic
    load(main_stim_path,"stimArray");
    disp('Loaded stimulus array.');
    toc
end

% Sort according to f0, azi and trajectory
stimArray = sortrows(stimArray,[2,9,8]);

% Each stimulus is duplicated (to be played with both intensities) - take
% only one copy of each here
stimArray = stimArray(1:2:end,:);

% Define the SPL value that was used to differentiate low and high
% intensity conditions
SPL_lo = db(0.2/2);
SPL_hi = db(1);

% Define the variables in the stimulus array (for indexing)
stimCol = 15;
stimCondCol = 8;
fsCol = 14;
% f0Col = 2;
% sideCol = 9;
% targCol = 10;
% congCol = 11;

% Define fs
fs = stimArray{1,fsCol};

% Define minimum length w/o 200 ms target/silence at the end
minL = min(cell2mat(stimArray(:,3))-0.2);

if ~exist("stl_loom_hi","var")
stl_loom_hi = []; stl_loom_lo = [];
stl_rec_hi = [];  stl_rec_lo = [];
stl_pps_hi = [];  stl_pps_lo = [];
stl_eps_hi = [];  stl_eps_lo = [];
end

startRow = rr+1;
for rr = startRow:size(stimArray,1)
    try
        stimulus_hi = 10^(SPL_hi/20)*stimArray{rr,stimCol}';
        stimulus_lo = 10^(SPL_lo/20)*stimArray{rr,stimCol}';

        % Take the stimulus w/o the last 200 ms (target or silence)
        stimulus_hi = stimulus_hi(1:end-(fs*0.2),:);
        stimulus_lo = stimulus_lo(1:end-(fs*0.2),:);

        % Take the last 2600 ms of the stimulus
        stimulus_hi = stimulus_hi(end-(fs*minL):end,:);
        stimulus_lo = stimulus_lo(end-(fs*minL):end,:);

        % Calculate loudness model for the loud stimulus
        disp(['Calculating STL for stimulus ',num2str(rr),'/400 (high int)...']);
        tic
        [stl,~,~] = moore2016(resample(stimulus_hi,2,3));
        toc

        % Put the short-term loudness result into the corresponding vector
        switch stimArray{rr,stimCondCol}
            case 1; stl_loom_hi(:,end+1) = stl;
            case 2; stl_rec_hi(:,end+1) = stl;
            case 3; stl_pps_hi(:,end+1) = stl;
            case 4; stl_eps_hi(:,end+1) = stl;
        end

        % Calculate loudness model for the low-intensity version of the same
        % stimulus
        disp(['Calculating stl for stimulus ',num2str(rr),'/400 (low int)...']);
        tic
        [stl,~,~] = moore2016(resample(stimulus_lo,2,3));
        toc

        % Put the stimulus into the corresponding vector
        switch stimArray{rr,stimCondCol}
            case 1; stl_loom_lo(:,end+1) = stl;
            case 2; stl_rec_lo(:,end+1) = stl;
            case 3; stl_pps_lo(:,end+1) = stl;
            case 4; stl_eps_lo(:,end+1) = stl;
        end

    catch
        warning(['There was an error for stimulus',num2str(rr)]);
        continue
    end
end

    % Average the short-term loudness estimates for each condition
    avg_stl_loom_hi = mean(stl_loom_hi,2);
    avg_stl_loom_lo = mean(stl_loom_lo,2);
    avg_stl_rec_hi = mean(stl_rec_hi,2);
    avg_stl_rec_lo = mean(stl_rec_lo,2);
    avg_stl_pps_hi = mean(stl_pps_hi,2);
    avg_stl_pps_lo = mean(stl_pps_lo,2);
    avg_stl_eps_hi = mean(stl_eps_hi,2);
    avg_stl_eps_lo = mean(stl_eps_lo,2);

    figure; hold on;
    plot(avg_stl_loom_hi,'r-','LineWidth',2);
    plot(avg_stl_loom_lo,'r--','LineWidth',2);
    plot(avg_stl_rec_hi,'b-','LineWidth',2);
    plot(avg_stl_rec_lo,'b--','LineWidth',2);
    plot(avg_stl_pps_hi,'y-','LineWidth',2);
    plot(avg_stl_pps_lo,'y--','LineWidth',2);
    plot(avg_stl_eps_hi,'g-','LineWidth',2);
    plot(avg_stl_eps_lo,'g--','LineWidth',2);
    legend(["loom hi","loom lo","rec hi","rec lo",...
        "pps hi","pps lo","eps hi","eps lo"],'Location','northeastoutside');

    save('loudnessModelResults\saved_variables.mat');
end