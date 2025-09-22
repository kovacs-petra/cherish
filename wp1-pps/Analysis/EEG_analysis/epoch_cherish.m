function epoch_cherish

% Epoching script for CherISH WP1
% Author: Petra Kovacs, Austrian Academy of Sciences, 2025

% % % Triggers: % % %
%
% EEG   GSR     Trajectory  Target  Offset-azi
% --------------------------------------------
% 151   171     loom        no      -90
% 152   172     loom        no       90
% 153           loom        yes     -90
% 154           loom        yes      90
% 155   175     rec         no      -90
% 156   176     rec         no       90
% 157           rec         yes     -90
% 158           rec         yes      90
% 159   179     PPS         no      -90
% 160   180     PPS         no       90
% 161           PPS         yes     -90
% 162           PPS         yes      90
% 163   183     EPS         no      -90
% 164   184     EPS         no       90
% 165           EPS         yes     -90
% 166           EPS         yes      90

diary epoching_diaries.txt; 

%% Flags
flags.eeg = 1;
flags.gsr = 0;

%% Epoch the EEG data
if flags.eeg
    dir_name = uigetdir; % 07_interpolated > high_int or low_int
    dir_file = dir(dir_name);
    num_file = size(dir_file);
    savepath = '\\KFS\Fileserver\ProjektDaten\CherISH\data\wp-1\EEG\08_epoched';

    if contains(dir_name,'high_int')
        intCond = 'hi';
    else
        intCond = 'lo';
    end

    addpath(dir_name)

    for i=3:num_file(1)
        filename = dir_file(i).name;
        set_e = strfind(filename, '.set');

        if (set_e > 0)
            subNum = filename(7:8);

            %% Load .set file
            global EEG;
            [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
            EEG = pop_loadset( 'filename', filename, 'filepath', dir_name);
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = eeg_checkset( EEG );

            %% Looming epochs
            EEG2 = pop_epoch( EEG, {'151'  '152'  '153'  '154'}, [-0.2 1], 'newname', 'loom', 'epochinfo', 'yes');
            [ALLEEG EEG2 CURRENTSET] = pop_newset(ALLEEG, EEG2, 1,'overwrite','on','gui','off');
            EEG2 = eeg_checkset( EEG2 );
            % EEG2 = pop_rmbase( EEG2, [-200 0]); % remove baseline
            EEG2 = pop_resample( EEG2, 100 ); % downsample to 100 Hz
            EEG2 = pop_saveset( EEG2,  'filename', strcat('s',subNum,'_loom_',intCond,'.set'), 'filepath', savepath); % save dataset

            % Reject artifacts (+- 100 uV)
            EEG2 = pop_eegthresh(EEG2,1,[1:4,6:9,11:20,22:26,28:63],-100,100,0,1,0,0); % except channels 27 21 10 5
            if ismember(0,EEG2.reject.rejthresh) % if there is anything NOT to reject
                EEG2 = pop_rejepoch(EEG2,[EEG2.reject.rejthresh],0);
                EEG2 = pop_saveset( EEG2,  'filename', strcat('s',subNum,'_loom_',intCond,'_artrej.set'), 'filepath', savepath);
            else
                warning(['Subject ',subNum,' has no epochs left after artefact rejection.']);
            end

            %% Receding epochs
            EEG3 = pop_epoch( EEG, {'155'  '156'  '157'  '158'}, [-0.2 1], 'newname', 'rec', 'epochinfo', 'yes');
            [ALLEEG EEG3 CURRENTSET] = pop_newset(ALLEEG, EEG3, 1,'overwrite','on','gui','off');
            EEG3 = eeg_checkset( EEG3 );
            % EEG3 = pop_rmbase( EEG3, [-200 0]); % remove baseline
            EEG3 = pop_resample( EEG3, 100 ); % downsample to 100 Hz
            EEG3 = pop_saveset( EEG3,  'filename', strcat('s',subNum,'_rec_',intCond,'.set'), 'filepath', savepath);% save dataset

            % Reject artifacts (+- 100 uV)
            EEG3 = pop_eegthresh(EEG3,1,[1:4,6:9,11:20,22:26,28:63],-100,100,0,1,0,0); %except channels 27 21 10 5
            if ismember(0,EEG3.reject.rejthresh) % if there is anything NOT to reject
                EEG3 = pop_rejepoch(EEG3,[EEG3.reject.rejthresh],0);
                EEG3 = pop_saveset( EEG3,  'filename', strcat('s',subNum,'_rec_',intCond,'_artrej.set'), 'filepath', savepath);
            else
                warning(['Subject ',subNum,' has no epochs left after artefact rejection.']);
            end

            %% Rotating in PPS epochs
            EEG4 = pop_epoch( EEG, {'159'  '160'  '161'  '162'}, [-0.2 1], 'newname', 'PPS', 'epochinfo', 'yes');
            [ALLEEG EEG4 CURRENTSET] = pop_newset(ALLEEG, EEG4, 1,'overwrite','on','gui','off');
            EEG4 = eeg_checkset( EEG4 );
            % EEG4 = pop_rmbase( EEG4, [-200 0]); % remove baseline
            EEG4 = pop_resample( EEG4, 100 ); % downsample to 100 Hz
            EEG4 = pop_saveset( EEG4,  'filename', strcat('s',subNum,'_pps_',intCond,'.set'), 'filepath', savepath);% save dataset

            % Reject artifacts (+- 100 uV)
            EEG4 = pop_eegthresh(EEG4,1,[1:4,6:9,11:20,22:26,28:63],-100,100,0,1,0,0); %except channels 27 21 10 5
            if ismember(0,EEG4.reject.rejthresh) % if there is anything NOT to reject
                EEG4 = pop_rejepoch(EEG4,[EEG4.reject.rejthresh],0);
                EEG4 = pop_saveset( EEG4,  'filename', strcat('s',subNum,'_pps_',intCond,'_artrej.set'), 'filepath', savepath);
            else
                warning(['Subject ',subNum,' has no epochs left after artefact rejection.']);
            end

            %% Rotating in EPS epochs
            EEG5 = pop_epoch( EEG, {'163'  '164'  '165'  '166'}, [-0.2 1], 'newname', 'EPS', 'epochinfo', 'yes');
            [ALLEEG EEG5 CURRENTSET] = pop_newset(ALLEEG, EEG5, 1,'overwrite','on','gui','off');
            EEG5 = eeg_checkset( EEG5 );
            % EEG5 = pop_rmbase( EEG5, [-200 0]); % remove baseline
            EEG5 = pop_resample( EEG5, 100 ); % downsample to 100 Hz
            EEG5 = pop_saveset( EEG5,  'filename', strcat('s',subNum,'_eps_',intCond,'.set'), 'filepath', savepath);% save dataset

            % Reject artifacts (+- 100 uV)
            EEG5 = pop_eegthresh(EEG5,1,[1:4,6:9,11:20,22:26,28:63],-100,100,0,1,0,0); %except channels 27 21 10 5
            if ismember(0,EEG5.reject.rejthresh) % if there is anything NOT to reject
                EEG5 = pop_rejepoch(EEG5,[EEG5.reject.rejthresh],0);
                EEG5 = pop_saveset( EEG5,  'filename', strcat('s',subNum,'_eps_',intCond,'_artrej.set'), 'filepath', savepath);
            else
                warning(['Subject ',subNum,' has no epochs left after artefact rejection.']);
            end

        end
    end
end

%% Epoch the GSR data
if flags.gsr
    dir_name = uigetdir; % 03_gsr > high_int or low_int      
    dir_file = dir(dir_name);
    num_file = size(dir_file);
    gsrpath = '\\KFS\Fileserver\ProjektDaten\CherISH\data\wp-1\EEG\08_epochedGSR';

    if contains(dir_name,'high_int')
        intCond = 'hi';
    else
        intCond = 'lo';
    end

    addpath(dir_name)

    for i=3:num_file(1)
        filename = dir_file(i).name;
        set_e = strfind(filename, '.set');

        if (set_e > 0)
            subNum = filename(7:8);

            %% Load .set file
            global EEG;
            [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
            EEG = pop_loadset( 'filename', filename, 'filepath', dir_name);
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = eeg_checkset( EEG );

            %% looming epochs
            EEG_gsr2 = pop_epoch( EEG, {'171'  '172'  '173'  '174'}, [-1 5], 'newname', 'loom', 'epochinfo', 'yes');
            [ALLEEG EEG_gsr2 CURRENTSET] = pop_newset(ALLEEG, EEG_gsr2, 1,'overwrite','on','gui','off');
            EEG_gsr2 = eeg_checkset( EEG_gsr2 );
            EEG_gsr2 = pop_rmbase( EEG_gsr2, [-1000 0]); % remove baseline
            EEG_gsr2 = pop_resample( EEG_gsr2, 100 ); % downsample to 100 Hz
            EEG_gsr2 = pop_saveset( EEG_gsr2,  'filename', strcat('s',subNum,'_loom_',intCond,'.set'), 'filepath', gsrpath); % save dataset

            %% Receding epochs
            EEG_gsr3 = pop_epoch( EEG, {'175'  '176'  '177'  '178'}, [-1 5], 'newname', 'rec', 'epochinfo', 'yes');
            [ALLEEG EEG_gsr3 CURRENTSET] = pop_newset(ALLEEG, EEG_gsr3, 1,'overwrite','on','gui','off');
            EEG_gsr3 = eeg_checkset( EEG_gsr3 );
            EEG_gsr3 = pop_rmbase( EEG_gsr3, [-1000 0]); % remove baseline
            EEG_gsr3 = pop_resample( EEG_gsr3, 100 ); % downsample to 100 Hz
            EEG_gsr3 = pop_saveset( EEG_gsr3,  'filename', strcat('s',subNum,'_rec_',intCond,'.set'), 'filepath', gsrpath);% save dataset

            %% Rotating in PPS epochs
            EEG_gsr4 = pop_epoch( EEG, {'179'  '180'  '181'  '182'}, [-1 5], 'newname', 'PPS', 'epochinfo', 'yes');
            [ALLEEG EEG4 CURRENTSET] = pop_newset(ALLEEG, EEG4, 1,'overwrite','on','gui','off');
            EEG_gsr4 = eeg_checkset( EEG_gsr4 );
            EEG_gsr4 = pop_rmbase( EEG_gsr4, [-1000 0]); % remove baseline
            EEG_gsr4 = pop_resample( EEG_gsr4, 100 ); % downsample to 100 Hz
            EEG_gsr4 = pop_saveset( EEG_gsr4,  'filename', strcat('s',subNum,'_pps_',intCond,'.set'), 'filepath', gsrpath);% save dataset

            %% Rotating in EPS epochs
            EEG_gsr5 = pop_epoch( EEG, {'183'  '184'  '185'  '186'}, [-1 5], 'newname', 'EPS', 'epochinfo', 'yes');
            [ALLEEG EEG_gsr5 CURRENTSET] = pop_newset(ALLEEG, EEG_gsr5, 1,'overwrite','on','gui','off');
            EEG_gsr5 = eeg_checkset( EEG_gsr5 );
            EEG_gsr5 = pop_rmbase( EEG_gsr5, [-1000 0]); % remove baseline
            EEG_gsr5 = pop_resample( EEG_gsr5, 100 ); % downsample to 100 Hz
            EEG_gsr5 = pop_saveset( EEG_gsr5,  'filename', strcat('s',subNum,'_eps_',intCond,'.set'), 'filepath', gsrpath);% save dataset
        end
    end
end

diary off;
