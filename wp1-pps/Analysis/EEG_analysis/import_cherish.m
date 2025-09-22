function import_cherish

%% Flags
flags.imp = 1; % import raw data into .set file?
flags.latCorr = 1; % apply latency correction?
flags.filt = 1; % apply bandpass filter on the EEG data?
flags.gsr = 1; % save and filter GSR data?
flags.trigCorr = [ones(1,10) zeros(1,16)]; % triggers need to be corrected? by subject
% flags.trigCorr = 0;

% diary;
dir_name = uigetdir; % select folder with raw data
dir_file = dir(dir_name);
num_file = size(dir_file);
addpath(dir_name)
saved_set = 0;

% Save paths have to be 'char', not "string"!
imp_folder = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\data\wp-1\EEG\01_imported\';
corr_folder = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\data\wp-1\EEG\01_triggerCorrected\';
lat_folder = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\data\wp-1\EEG\02_latencyCorrected\';
filt_folder = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\data\wp-1\EEG\03_filtered\';
gsr_folder = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\data\wp-1\EEG\03_gsr\';

% EEGLab path
% addpath 'C:\Users\pkovacs\Documents\MATLAB\eeglab2024.2'

for i=3:num_file(1)
    filename = dir_file(i).name;
    search = strfind(filename, '.vhdr'); % if newly importing raw data
    % search = strfind(filename, '.set'); % if data already imported before
    if (search > 0)
        global EEG;
        current_subj = str2double(filename(end-6:end-5));
        [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

        %% Importing
        if flags.imp
        % Load BrainVision file; 1st 64 electrodes plus StimTrak, GSR
        EEG = pop_loadbv(dir_name, filename, []);
        % EEG = pop_loadset('filename',filename,'filepath', dir_name);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname','','gui','off');
        % EEG = pop_loadset( 'filename', filename, 'filepath', dir_name);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = eeg_checkset( EEG );

        % load("chanlocs66.mat", "chanlocs");
        % pop_chanedit(chanlocs);
        % EEG=pop_chanedit(EEG, 'changefield',{64,'theta',''},'changefield',{64,'radius',''},'changefield',{64,'X',''},'changefield',{64,'Y',''},'changefield',{64,'Z',''},'changefield',{64,'sph_theta',''},'changefield',{64,'sph_phi',''},'changefield',{64,'sph_radius',''},'changefield',{64,'X','0.37049'},'changefield',{64,'Y','0'},'changefield',{64,'Z','0.85717'});
        % EEG = pop_chanedit(EEG, 'changefield',{64,'labels','FCz'},'changefield',{64,'X','0.37049'},'changefield',{64,'Y','0'},'changefield',{64,'Z','0.85717'});

        % save imported
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG,  'filename', strcat(filename(1:end-5),'.set'), 'filepath', imp_folder);
        end

        %% Correcting trigger error
        EEG = pop_loadset('filename',strcat(filename(1:end-5),'.set'),'filepath', imp_folder);
        if flags.trigCorr(current_subj) 
            disp(['Subject ',num2str(current_subj),' had too many triggers sent out, we correct the mistake.']);
            
            EEG = correctTriggerMistake(EEG,current_subj);

            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG,  'filename', strcat(filename(1:end-5),'_trCorr.set'), 'filepath', corr_folder);
        else % remove "T" entry from trigger type
            tmp = {EEG.event.type}; % Initial char array
            tmp = cellfun(@(x) x(2:end),tmp,'UniformOutput',false); % Reject first entry ("T")
            tmp = str2double(tmp);

            tmp(isnan(tmp)) = '0';
            tmpEventList = tmp;

            % Keep the adjusted stimuli in the EEG list
            for ii = 1:length(tmpEventList)
            	EEG.event(ii).type = tmpEventList(:,ii);
            end
            EEG = pop_saveset( EEG,  'filename', strcat(filename(1:end-5),'_trCorr.set'), 'filepath', corr_folder);
        end

        %% Latency correction based on StimTrak
        if flags.latCorr
            
            % Adjust trigger latencies based on the StimTrak channel
            disp('Latency correction requested.')
            EEG = pop_loadset('filename',strcat(filename(1:end-5),'_trCorr.set'),'filepath', corr_folder);
            EEG = latencyCorrection(EEG);
            
            % save latency corrected
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG,  'filename', strcat(filename(1:end-5),'_latCorr','.set'), 'filepath', lat_folder);
        end

        %% Bandpass filtering
        if flags.filt
            
            % Separate EEG and GSR datasets
            disp('Filtering requested.')
            EEG = pop_loadset('filename',strcat(filename(1:end-5),'_latCorr','.set'),'filepath', lat_folder);
            
            EEG_GSR = pop_select(EEG,'channel',{'GSR'}); % Keep only this channel
            [ALLEEG, EEG_GSR, CURRENTSET] = eeg_store( ALLEEG, EEG_GSR, 0 );
            EEG_GSR = eeg_checkset( EEG_GSR );

            EEG = pop_select(EEG,'rmchannel',{'StimTrak','GSR'}); % Given EEG dataset, return it excluding channel StimTrak and GSR
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname','','gui','off');
            EEG = eeg_checkset( EEG );
            eeglab redraw;

            %% EEG channels
            % Bandpass filter EEG signal
            EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'hicutoff',30,'plotfreqz',0);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG = eeg_checkset( EEG );

            % Calculate average reference excluding eye channels; add a channel of zeros and the ref channel first
            num_chans = EEG.nbchan;

            if current_subj >= 12
                % loc_file = 'C:\Users\pkovacs\Documents\MATLAB\eeglab2024.2\plugins\dipfit\standard_BEM\elec\standard_1020.elc';
                loc_file = 'C:\Users\pkovacs\Documents\MATLAB\eeglab2025.0.0\plugins\dipfit\standard_BEM\elec\standard_1020.elc';
                EEG = pop_chanedit(EEG, 'append', num_chans, 'changefield', {num_chans+1, 'labels', 'FCz'}, 'lookup', loc_file, 'setref', {1:num_chans, 'FCz'});
                EEG.data(end+1,:) = zeros(1,length(EEG.data));
                EEG = pop_reref( EEG, [],'exclude',[5 10 21 27],'keepref','on','refloc', struct('labels', {'FCz'}, 'type', {''}, 'theta', {0.7867}, 'radius', {0.095376}, 'X', {27.39}, 'Y', {-0.3761}, 'Z', {88.668}, 'sph_theta', {-0.7867}, 'sph_phi', {72.8323}, 'sph_radius', {92.8028}, 'urchan', {num_chans+1}, 'ref', {''}, 'datachan', {0}));
                [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
                EEG = eeg_checkset( EEG );

            else
                load("chanlocs.mat", "chanlocs");
                pop_chanedit(chanlocs);
                EEG = pop_chanedit(EEG, 'changefield',{64,'labels','FCz'},'changefield',{64,'X','0.37049'},'changefield',{64,'Y','0'},'changefield',{64,'Z','0.85717'});
                EEG = pop_select(EEG, 'channel', 1:64);
                EEG.data(end+1,:) = zeros(1,length(EEG.data));
                EEG = pop_reref( EEG, [],'exclude',[5 10 21 27],'keepref','on');
            end


            % Save filtered EEG data
            EEG = pop_saveset( EEG,  'filename', strcat(filename(1:end-5),'_filt','.set'), 'filepath', filt_folder);
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        end

        if flags.gsr
            
            %% GSR channel
            % Bandpass filter signal (0.1-5 Hz after Bach, 2009)
            EEG_GSR = pop_eegfiltnew(EEG_GSR, 'locutoff',0.1,'hicutoff',5,'plotfreqz',1);
            [ALLEEG EEG_GSR CURRENTSET] = pop_newset(ALLEEG, EEG_GSR, 1,'overwrite','on','gui','off');
            [ALLEEG EEG_GSR] = eeg_store(ALLEEG, EEG_GSR, CURRENTSET);
            EEG_GSR = eeg_checkset( EEG_GSR );

            % Save filtered GSR data
            % [ALLEEG EEG_GSR] = eeg_store(ALLEEG, EEG_GSR, CURRENTSET);
            % EEG_GSR = eeg_checkset( EEG_GSR );
            EEG_GSR = pop_saveset( EEG_GSR,  'filename', strcat(filename(1:end-5),'_gsr','.set'), 'filepath', gsr_folder);
            [ALLEEG EEG_GSR] = eeg_store(ALLEEG, EEG_GSR, CURRENTSET);
        end
    end

end

if (saved_set == 0)

    disp ('DONE');
    diary off;

end