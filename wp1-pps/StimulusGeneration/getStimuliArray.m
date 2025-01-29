function stimArray = getStimuliArray(conditions)
% Function to summarize CherISH pilot stimuli into one big stimuli mat file
%
% USAGE: stimArray = getStimuliArray(folders)
%
% Loads all audio data and corresponding parameters from the supplied list
% of folders and concatenates them into one cell matrix. 
%
% Input:
% folders       - Cell array of folder names/paths. Folders are expected to
%               be created with generateStimuli script,
%               that is, to contain wav files (one per stimulus) and a csv
%               file with metadata. If there is only one folder, simple
%               char array is sufficient as well.
% conditions    - 1: looming and receding blocks, 2: rotating blocks
%
% Output:
% stimArray    - Cell array of cells, each with size "no. of stimuli" X "no. of parameters". 
%               The last column contains raw audio data matrices
%
% Notes         (1) We expect the param csv files to have the columns
%               specified in paramFields; check if that's really what you
%               need
%

%% Input check
if conditions == 1
    folders = {'23-Jan-2025-1411','23-Jan-2025-1412'}; % task 1, loom and rec
else
    folders = {'23-Jan-2025-1336','23-Jan-2025-1346'}; % task 2, near and far
end

% check for existence of folders
if iscell(folders)
    for i = 1:length(folders)
        if ~exist(folders{i}, 'dir')
            error(['Could not find target folder ', folders{i}, '!']);
        end
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
    'frequency', ...         % Cue frequency in Hz
    'totalDur', ...          % Cue duration in s
    'durStatOnset', ...      % Duration of stationary onset
    'durStatOffset', ...     % Duration of stationary offset
    'onsetDistance', ...     % Cue onset distance in m
    'offsetDistance', ...    % Cue offset distance in m
    'trajectory', ...        % 1 - loom, 2 - rec, 3 - rotate near, 4 - rotate far
    'offsetAzimuth', ...     % Side of the cue (at offset): 90 - left, -90 - right
    'target', ...            % 1 - target trial, 0 - nontarget trial
};

% output variable collecting stimuli sets
stimArray = cell(length(folders), 1);
repeat = 5; % repeat each stimulus in the stimArray

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
    audioF = [myCell{1, 1}, '-full.mat'];
    % [data, ~] = audioread(audioF);
    load(audioF, 'scalestimT');
    
    % preallocate all memory we need to include the audio as well (add
    % extra column for raw audio data)
    myCell = [myCell, repmat({zeros(size(scalestimT))}, size(myCell, 1), 1)];    
    myCell = [myCell; repmat(myCell,repeat,1)];
    
    %% Loop through mat files containing audio data   
    for audio = 1:size(myCell, 1)
        % load audio into cell array
        audioF = [myCell{audio, 1}, '-full.mat'];
        load(audioF, 'scalestimT');
        myCell{audio, length(paramFields)+1} = scalestimT;

    end  % audio files for loop

    % accumulate all data into final 
    stimArray{f, 1} = myCell;    
    disp(['Done with folder ', folders{f}]);    
    
end  % folders for loop

stimArray = vertcat(stimArray{:});

%% Ending
stimArrayName = ['stimArray', num2str(conditions),'.mat'];
save(stimArrayName,"stimArray");
disp([char(10), 'Finished, returning']);

return
    
    
    
    
    
    
    
    
    
    
    
    
    