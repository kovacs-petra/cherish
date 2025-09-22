dir_name = uigetdir; 
dir_file = dir(dir_name);
num_file = size(dir_file);
addpath(dir_name)
saved_set = 0;

% br_folder = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\data\wp-1\EEG\04_breaksRejected\';
ica_folder = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\data\wp-1\EEG\05_ica\';

% EEGLab path
% addpath 'C:\Users\pkovacs\Documents\MATLAB\eeglab2024.2'

for i=3:num_file(1)
    filename = dir_file(i).name;
    search = strfind(filename, '.set'); % if data already imported before
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
        EEG = pop_saveset( EEG,  'filename', strcat(filename(1:8),'_ica','.set'), 'filepath', ica_folder);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    end
end
