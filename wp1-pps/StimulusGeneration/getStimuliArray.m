function stimArray = getStimuliArray
% Function to summarize CherISH pilot stimuli into one big stimuli mat file
%
% USAGE: stimArray = getStimuliArray(folders)
%
% Loads all audio data and corresponding parameters from the supplied list
% of folders and concatenates them into one cell matrix. 
%
% Output:
% stimArray    - Cell array of cells, each with size "no. of stimuli" X "no. of parameters". 
%               The last column contains raw audio data matrices
%
% Notes         (1) We expect the param csv files to have the columns
%               specified in paramFields; check if that's really what you
%               need
%

folders = {'30-May-2025-114','30-May-2025-1126','30-May-2025-1149',...
        '30-May-2025-1211'}; 

% check for existence of folders
for i = 1:length(folders)
    if ~exist(folders{i}, 'dir')
        error(['Could not find target folder ', folders{i}, '!']);
    end
end

% user message
disp([char(10), 'Called function getStimuliArray with inputs: ', char(10),...
    'folders: ']);
for i = 1:length(folders)
    disp(folders{i});
end

%% Loop through folders, load the parameters file and audio data from each
% fields we expect in params files
paramFields = {
    'filename', ...
    'f0', ...                % Cue f0 in Hz
    'totalDur', ...          % Cue duration in s
    'durStatOnset', ...      % Duration of stationary onset in the cue
    'onsetDistance', ...     % Cue onset distance in m
    'offsetDistance', ...    % Cue offset distance in m
    'direction',...          % 1 - radial, 2 - angular
    'trajectory', ...        % 1 - loom, 2 - rec, 3 - rotate near, 4 - rotate far
    'offsetAzimuth', ...     % Side of the cue (at offset): 90 - left, -90 - right
    'targetTrial', ...       % 1 - target trial, 0 - nontarget trial
    'congruence', ...        % 1 - congruent target, 0 - incongruent target
    'targetAzimuth', ...     % 90 or -90 (deg)
    'sourceInt', ...         % 1 - high source intensity, 0 - low source intensity (placeholder)
    'fs', ...                % sampling rate
};

% output variable collecting stimuli sets
stimArray = cell(length(folders), 1);
% repeat = 4; % repeat each stimulus in the stimArray

for f = 1:length(folders)    
    disp([char(10), 'Loading params and audio from: ', folders{f}]);
    
    % find corresponding csv file
    csvFile = dir([folders{f},'/*.csv']);
    paramFile = [folders{f}, '/', csvFile.name];
    % load params from csv file
    T = readtable(paramFile, 'Delimiter', ',');
    % transform into cell
    myCell = table2cell(T);
    
    % check params fields (and implicitly, size)
    if ~isequal(T.Properties.VariableNames, paramFields)
        error(['CSV params file ', csvFile, ' had unexpected variables/columns']);
    end

    % load first audio to get audio size
    % audioF = [folders{f}, '/', myCell{1, 1}, '.mat'];
    audioF = [myCell{1,1}, '.mat']; 
    % [data, ~] = audioread(audioF);
    load(audioF, 'out');
    
    % preallocate all memory we need to include the audio as well (add
    % extra column for raw audio data)
    myCell = [myCell, repmat({zeros(size(out))}, size(myCell, 1), 1)];    
    % myCell = [myCell; repmat(myCell,repeat,1)];
    
    %% Loop through mat files containing audio data   
    for audio = 1:size(myCell, 1)
        % load audio into cell array
        audioF = [myCell{audio, 1}, '.mat'];
        load(audioF, 'out');
        myCell{audio, length(paramFields)+1} = out;

    end  % audio files for loop

    % accumulate all data into final 
    stimArray{f, 1} = myCell;    
    disp(['Done with folder ', folders{f}]);    
    
end  % folders for loop

stimArray = vertcat(stimArray{:});

%% Ending
stimArrayName = 'stimArray.mat';
save(stimArrayName,"stimArray");
disp([char(10), 'Finished, returning']);

return
    
    
    
    
    
    
    
    
    
    
    
    
    