function ica(paths)

% Perform ICA decpmposition in EEGlab. 
%
% Input:
% paths - struct with field "ica_folder"

dir_name = uigetdir; 
dir_file = dir(dir_name);
num_file = size(dir_file);
addpath(dir_name)
saved_set = 0;

for i=3:num_file(1)
    filename = dir_file(i).name;
    search = strfind(filename, '.set'); 
    if (search > 0)
        global EEG;
        % current_subj = str2double(filename(end-6:end-5));
        [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

        %% Load set
        EEG = pop_loadset('filename',filename,'filepath',dir_name);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname','','gui','off');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = eeg_checkset( EEG );

        %% ICA decomposition
        % EEG = pop_runica(EEG, 'icatype','runica');
        % EEG = pop_runica(EEG, 'pca',32,'interupt','on');
        EEG = pop_runica(EEG, 'extended',1,'interupt','on');
        EEG = pop_saveset( EEG,  'filename', strcat(filename(1:8),'_ica','.set'), 'filepath', paths.ica_folder);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    end
end
