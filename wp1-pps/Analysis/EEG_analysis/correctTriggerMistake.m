function EEG = correctTriggerMistake(EEG,subNum)
% Initial setup
logIndex = 2; % Skip header row of logVar
deleteIndices = [];        % To collect indices of events to delete
changefieldList = {};      % For changing field values
insertList = {};           % For inserting new events

% Load subject's logVar
logfile = strcat('\\kfs\fileserver\Projektdaten\CherISH\data\wp-1\Behav\subject',...
    num2str(subNum),'\sub',num2str(subNum),'Log.mat');
load(logfile,"logVar");

% Adjust event list: Convert to binary and take only last 6 bits
tmp = {EEG.event.type}; % Initial char array representing the binary
tmp = cellfun(@(x) x(2:end),tmp,'UniformOutput',false); % Reject first entry ("T")
tmp = str2double(tmp);

tmp(isnan(tmp)) = '0';
% tmpEventList = dec2bin(tmp);
tmpEventList = tmp;

% Keep the binary stimuli in the EEG list
for ii = 1:length(tmpEventList)
	EEG.event(ii).type = tmpEventList(:,ii);
end

i = 1;
while i <= length(EEG.event)
    % Convert type to double if needed
    if ischar(EEG.event(i).type)
        thisType = str2double(EEG.event(i).type);
    else
        thisType = EEG.event(i).type;
    end

    if thisType == 210
        % Next event exists and logVar still has values
        if i + 1 <= length(EEG.event) && logIndex <= size(logVar, 1)
            % Replace next event's type with value from logVar
            newType = logVar{logIndex, 20};
            newLat = EEG.event(i+1).latency / EEG.srate;

            % changefieldList{end+1} = { i+1, 'type', newType };
            EEG = pop_editeventvals(EEG,'changefield',{i+1,'type',newType}, ...
                'changefield',{i+1,'latency',newLat});

            % If there is no target in the trial,
            % insert an additional event of type +20 (i.e., GSR trigger)
            if logVar{logIndex, 12} == 0
                % Create a new event copied from i+1
                newEvent = EEG.event(i+1);
                newEvent.type = newType + 20;

                % Make sure latency is defined (in seconds!)
                if isfield(newEvent, 'latency')
                    latencyInSec = newEvent.latency / EEG.srate;
                else
                    latencyInSec = (i+1) / EEG.srate; % fallback
                end

                % insertList{end+1} = { i+1, 'latency',latencyInSec,'duration',0.001,'type',newType+20 };
                EEG = pop_editeventvals(EEG,'insert',{i+1,[],[],[],[],[],[],[],[],[]}, ...
                    'changefield',{i+1,'latency',latencyInSec}, ...
                    'changefield',{i+1,'duration',0.001}, ...
                    'changefield',{i+1,'type',newType+20}, ...
                    'changefield',{i+1,'code','T'}, ...
                    'changefield',{i+1,'channel',0}, ...
                    'changefield',{i+1,'bvmknum',i+1});
                N = 3; 
            else
                N = 2;
            end

            logIndex = logIndex + 1;
        end

        i = i + N;

        % Delete events until next type >= 200
        while i <= length(EEG.event)
            if ischar(EEG.event(i).type)
                t = str2double(EEG.event(i).type);
            else
                t = EEG.event(i).type;
            end

            if t < 140 || t >= 199
                break;  % Don't delete this one
            else
                deleteIndices(end+1) = i;
                i = i + 1;
            end
        end

    elseif thisType == 200
        % If next exists and is not 210, delete it
        if i + 1 <= length(EEG.event)
            if ischar(EEG.event(i+1).type)
                nextType = str2double(EEG.event(i+1).type);
            else
                nextType = EEG.event(i+1).type;
            end
            if nextType ~= 210
                deleteIndices(end+1) = i + 1;
            end
        end
        i = i + 1;

    else
        i = i + 1;
    end
end

% Apply changes via pop_editeventvals
% if ~isempty(changefieldList)
%     for k = 1:length(changefieldList)
%         EEG = pop_editeventvals(EEG, 'changefield', changefieldList{k});
%         disp(['Changing field ',num2str(k),'/',num2str(length(changefieldList)),'...']);
%     end
% end

% if ~isempty(insertList)
%     for k = 1:length(insertList)
%         % insertList: {index, type, latency (in seconds)}
%         EEG = pop_editeventvals(EEG, 'insert', insertList{k});
%         disp(['Inserting trigger ',num2str(k),'/',num2str(length(insertList)),'...']);
% 
%     end
% end

if ~isempty(deleteIndices)
    EEG = pop_editeventvals(EEG, 'delete', deleteIndices);
    disp(['Deleting ',num2str(length(deleteIndices)),' triggers...']);
end

% Re-check consistency
EEG = eeg_checkset(EEG, 'eventconsistency');
eeglab redraw;