function ft2fif
% Convert FieldTrip ERP files to MNE Python .fiff files

%% Set up paths
path_epoched_set = '\\KFS\Fileserver\ProjektDaten\CherISH\data\wp-1\EEG\08_epoched';
path_epoched_ft = '\\kfs\fileserver\Projektdaten\CherISH\data\wp-1\EEG\08_epoched_ft\';
path_ft = 'C:\Users\pkovacs\Documents\MATLAB\fieldtrip-20250523';
path_fiff = '\\kfs\fileserver\Projektdaten\CherISH\data\wp-1\EEG\09_fifData\';
path_eeglab = 'C:\Users\pkovacs\Documents\MATLAB\eeglab2024.2';
dir_epoched_set = dir(path_epoched_set);
addpath(path_eeglab,path_epoched_ft,path_ft,path_fiff);
eeglab;

for dd = 3:length(dir_epoched_set) % for each epoched file
    filename = dir_epoched_set(dd).name;

    if contains(filename, 'artrej.set')
        ss = str2double(filename(2:3)); % Subject number (sXX)
        if ss > 12

            % Load EEGlab .set file
            EEG = pop_loadset('filename',filename,'filepath', path_epoched_set);

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
            if contains(filename,'loom_hi')
                loom_hi = ft_timelockanalysis(cfg,dataFT);
                savename = strcat(filename(1:end-11), '.mat');
                save(strcat(path_epoched_ft,savename),'loom_hi');

            elseif contains(filename,'rec_hi')
                rec_hi = ft_timelockanalysis(cfg,dataFT);
                savename = strcat(filename(1:end-11), '.mat');
                save(strcat(path_epoched_ft,savename),'rec_hi');

            elseif contains(filename,'loom_lo')
                loom_lo = ft_timelockanalysis(cfg,dataFT);
                savename = strcat(filename(1:end-11), '.mat');
                save(strcat(path_epoched_ft,savename),'loom_lo');

            elseif contains(filename,'rec_lo')
                rec_lo = ft_timelockanalysis(cfg,dataFT);
                savename = strcat(filename(1:end-11), '.mat');
                save(strcat(path_epoched_ft,savename),'rec_lo');

            end

            % Save subject erp's as .fiff, which is the input to the Python analysis
            fiff_name = strcat(filename(1:end-4),'-ave.fif');
            if contains(filename,'loom_hi')
                fieldtrip2fiff(strcat(path_fiff, fiff_name), loom_hi);
            elseif contains(filename,'rec_hi')
                fieldtrip2fiff(strcat(path_fiff, fiff_name), rec_hi);
            elseif contains(filename,'loom_lo')
                fieldtrip2fiff(strcat(path_fiff, fiff_name), loom_lo);
            elseif contains(filename,'rec_lo')
                fieldtrip2fiff(strcat(path_fiff, fiff_name), rec_lo);
            end
        end
    end
end
