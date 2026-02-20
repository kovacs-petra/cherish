% Collect trial-by-trial amplitude maxima and store them together with the
% behavioral (promptness) results. CherISH WP1.

% #Author: Petra Kovács, Austrian Academy of Sciences, 2026

%% Initialize required files and tools
% Open EEGlab
path_eeglab = '\\kfs.oeaw.ac.at\Fileserver\Projektdaten\CherISH\tools\eeglab2025.0.0';
addpath(path_eeglab);
eeglab;

% Load file containing behavioral (promptness) results
if exist("promptness_ERP_artrej.csv","file")
    B = readtable("promptness_ERP_artrej.csv");
else
    path_behav = '\\KFS.oeaw.ac.at\Fileserver\Projektdaten\CherISH\data\wp-1\Behav\main_task\mainTaskData_dbFS.csv';
    B = readtable(path_behav);

    % create columns for urevent and mean ERP
    B.urevent = zeros(height(B),1);
    B.meanERP = nan(height(B),1);
end

% Define included subjects
allSub = 3:24;

% Define motion onset EEG triggers for target trials
trigs = [153,154,157,158,161,162,165,166];

%% Tag each trial in the behavioral data with the urevent number from the
%% original (latency-corrected) EEG recording
% For each EEG file
lat_folder = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\data\wp-1\EEG\02_latencyCorrected\';
lat_file = dir(lat_folder);
int_folder = '\\KFS\Fileserver\ProjektDaten\CherISH\data\wp-1\EEG\07_interpolated';
int_file = dir(int_folder);
addpath(lat_folder,int_folder);

% Use the channel and time clusters that are significant in the looming-receding
% comparison
ch = [ 1,  5, 23, 24, 27, 29, 30, 31, 52, 53, 54, 55, 56] + 1; % copied from python, but that starts with 0
sig_onset = 1440; % ms
stim_offset = 2000; % significance continues, but we only go until end of stimulus

for ff = 3:size(lat_file,1)
    filename_latcorr = lat_file(ff).name;
    search = strfind(filename_latcorr, '.set');

    if search > 0
        global EEG;
        current_sub = str2double(filename_latcorr(7:8));
        if ismember(current_sub,allSub)
            % Load latency-corrected EEG recording
            EEG = pop_loadset('filename',filename_latcorr,'filepath', lat_folder);
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname','','gui','off');
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = eeg_checkset( EEG );

            block1_start_trigs = find([EEG.event.type] == 101); % trigger for 1st block
            if length(block1_start_trigs) > 1
                block1_start = block1_start_trigs(end); % some Ss restarted the 1st block
            else
                block1_start = 1;
            end

            % Define rows in B which belong to target trials
            rows = find(B.subNum == current_sub & B.targetTrial == 1);
            rr = 1; % row of logVar
            % Look for the next motion onset trigger
            for ee = block1_start:length(EEG.event)
                if ismember(EEG.event(ee).type, trigs)
                    % Verify that the condition is the same as in the next promptness datapoint
                    switch EEG.event(ee).type
                        case {trigs(1),trigs(2)}
                            current_cond = 1; % loom
                        case {trigs(3),trigs(4)}
                            current_cond = 2; % rec
                        case {trigs(5),trigs(6)}
                            current_cond = 3; % pps
                        case {trigs(7),trigs(8)}
                            current_cond = 4; % eps
                    end

                    if isequal(current_cond,B.trajectory(rows(rr)))
                        % Put event number in the last column of B
                        B.urevent(rows(rr)) = EEG.event(ee).urevent;
                        rr = rr+1;
                    else
                        error(['EEG and behav conditions don''t match.',...
                            'Subject: ',num2str(current_sub),'; '...
                            'Urevent: ',num2str(EEG.event(ee).urevent)]);
                    end
                end
            end

            % Clear latency-corrected data from memory
            ALLEEG = []; EEG = []; CURRENTSET = [];

            % Get interpolated file for the same subject
            subNumString = num2str(current_sub,'%02d');
            filename_interp = strcat('\WP1_00',subNumString,'_interp.set');

            %% Extract amplitude maxima from EEG data
            % Load interpolated data
            EEG = pop_loadset('filename',filename_interp,'filepath', int_folder);
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname','','gui','off');
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = eeg_checkset( EEG );

            block1_start_trigs = find(contains([EEG.event.type],'101'));
            if length(block1_start_trigs) > 1
                block1_start = block1_start_trigs(end); % some Ss restarted the 1st block
            else
                block1_start = 1;
            end

            for rr = 1:length(rows)
                % In the EEG data, look for the urevent which corresponds to the one found
                % in the promptness data
                for ee = block1_start:length(EEG.event)
                    if isequal(EEG.event(ee).urevent,B.urevent(rows(rr)))
                        start_event = EEG.event(ee).latency;
                        break
                    end
                end

                % Baseline correct the EEG epoch
                EEG.data(:,round(start_event:start_event+stim_offset)) = ...
                    pop_rmbase(EEG, [-200 0]);

                % Calculate mean amplitude and log it in B
                t = round(start_event+sig_onset:start_event+stim_offset);

                if length(EEG.data) < t(end)
                    break
                else
                    a = mean(EEG.data(ch, t), [1,2]);
                end

                if B.urevent(rows(rr)) == 0
                    a = nan; % handle missing blocks
                end

                % Reject any artifacts > +-100 uV
                noEyeChannels = [1:4,6:9,11:20,22:26,28:63];
                if ~isempty(find(abs(EEG.data(noEyeChannels,t))>100))
                    a = nan;
                end

                B.meanERP(rows(rr)) = a;
            end
        end % current subject
    end
%% Save .csv
writetable(B,"promptness_ERP_artrej.csv");
end
