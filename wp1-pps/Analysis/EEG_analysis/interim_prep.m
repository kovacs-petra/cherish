function interim_prep
% Prepare CherISH WP-1 data for interim analysis in Python

%% Set up paths
path_eeglab = 'C:\Users\pkovacs\Documents\MATLAB\eeglab2024.2';
path_epoched_data = '\\kfs\fileserver\Projektdaten\CherISH\data\wp-1\EEG\08_epoched\PPSlo-EPShi\';
path_ft = 'C:\Users\pkovacs\Documents\MATLAB\fieldtrip-20250523';
dir_epoched_data = dir(path_epoched_data);
addpath(path_eeglab,path_epoched_data,path_ft);
eeglab;

savedir = '\\kfs\fileserver\Projektdaten\CherISH\data\wp-1\EEG\08_epoched_ft\PPSlo-EPShi\';

%% Initialize variables
% PPS_hi = {};
% EPS_hi = {};
% PPS_lo = {};
% EPS_lo = {};

for dd = 3:length(dir_epoched_data) % for each epoched file
    filename = dir_epoched_data(dd).name;

    % Work only on artefact rejected .set files
    if contains(filename,'.set') && contains(filename, 'artrej')
        ss = str2double(filename(2:3)); % Subject number (sXX)

        % Load EEGlab .set file
        EEG = pop_loadset('filename',filename,'filepath', path_epoched_data);
        
        % Convert EEGlab to FieldTrip
        dataFT = eeglab2fieldtrip(EEG,'preprocessing');

        % Baseline correction
        cfg = [];
        cfg.demean = 'yes';
        cfg.baselinewindow = [-0.2 0];
        [dataFT] = ft_preprocessing(cfg, dataFT);

        % Create configuration variable
        cfg = [];
        % cfg.keeptrials = 'no';
        % cfg.trials = find(ismember(string(dataFT.trialinfo{:,8}), {'163', '164', '165', '166'}));

        % Calculate avg ERPs for each condition
        % Choose smaller epoch: e.g. -100 to 500 ms?
        if contains(filename,'pps_hi')
            PPS_hi = ft_timelockanalysis(cfg,dataFT);
            savename = strcat(filename(1:end-11), '.mat');
            save(strcat(savedir,savename),'PPS_hi');

        elseif contains(filename,'eps_hi')
            EPS_hi = ft_timelockanalysis(cfg,dataFT);
            savename = strcat(filename(1:end-11), '.mat');
            save(strcat(savedir,savename),'EPS_hi');

        elseif contains(filename,'pps_lo')
            PPS_lo = ft_timelockanalysis(cfg,dataFT);
            savename = strcat(filename(1:end-11), '.mat');
            save(strcat(savedir,savename),'PPS_lo');

        elseif contains(filename,'eps_lo')
            EPS_lo = ft_timelockanalysis(cfg,dataFT);
            savename = strcat(filename(1:end-11), '.mat');
            save(strcat(savedir,savename),'EPS_lo');

        end

        % Save subject erp's as .fiff, which is the input to the Python analysis
        fiff_name = strcat(filename(1:end-11),'-ave.fif'); % has to end in -ave.fif
        fiff_path = '\\kfs\fileserver\Projektdaten\CherISH\data\wp-1\EEG\09_fifData\PPSlo-EPShi\';
        if contains(filename,'pps_hi')
            fieldtrip2fiff(strcat(fiff_path, fiff_name), PPS_hi);
        elseif contains(filename,'eps_hi')
            fieldtrip2fiff(strcat(fiff_path, fiff_name), EPS_hi);
        elseif contains(filename,'pps_lo')
            fieldtrip2fiff(strcat(fiff_path, fiff_name), PPS_lo);
        elseif contains(filename,'eps_lo')
            fieldtrip2fiff(strcat(fiff_path, fiff_name), EPS_lo);
        end
  
    end
end

%% Plotting group data
% Set up 
dir_ft_data = dir(savedir);

for gg = 3:length(dir_ft_data)
    filename = dir_ft_data(gg).name;
    ss = str2double(filename(2:3));
    folder = dir_ft_data(gg).folder;

    if contains(filename,'pps_hi')
            load(strcat(folder,'\',filename),'PPS_hi');
            PPS_hi_all(ss) = {PPS_hi};

    elseif contains(filename,'eps_hi')
            load(strcat(folder,'\',filename),'EPS_hi');
            EPS_hi_all(ss) = {EPS_hi};

    elseif contains(filename,'pps_lo')
            load(strcat(folder,'\',filename),'PPS_lo');
            PPS_lo_all(ss) = {PPS_lo};

    elseif contains(filename,'eps_lo')
            load(strcat(folder,'\',filename),'EPS_lo');
            EPS_lo_all(ss) = {EPS_lo};

    end

end

PPS_hi_all = PPS_hi_all(3:end);
EPS_hi_all = EPS_hi_all(3:end);
PPS_lo_all = PPS_lo_all(3:end);
EPS_lo_all = EPS_lo_all(3:end);

% For plotting: calculate grand avg ERPs for each of the four conditions
% Warning: this only works if cfg.keeptrials was set to 'no' in the subject
% level data
cfg = [];
cfg.channel = 'all';
cfg.parameter = 'avg';
gravg_PPS_hi = ft_timelockgrandaverage(cfg, PPS_hi_all{:});
gravg_EPS_hi = ft_timelockgrandaverage(cfg, EPS_hi_all{:});
gravg_PPS_lo = ft_timelockgrandaverage(cfg, PPS_lo_all{:});
gravg_EPS_lo = ft_timelockgrandaverage(cfg, EPS_lo_all{:});

gravg_PPS = ft_timelockgrandaverage(cfg, PPS_hi_all{:}, PPS_lo_all{:});
gravg_EPS = ft_timelockgrandaverage(cfg, EPS_hi_all{:}, EPS_lo_all{:});

cfg = [];
cfg.channel = [24];
cfg.xlim = [-0.2 1];
cfg.ylim = [-2.5 3];
cfg.showlegend = 'yes';

% figure;
% ft_singleplotER(cfg, gravg_PPS_lo, gravg_EPS_lo);
% figure;
% ft_singleplotER(cfg, gravg_PPS_hi, gravg_EPS_hi);

figure;
ft_singleplotER(cfg, gravg_PPS, gravg_EPS);


    

