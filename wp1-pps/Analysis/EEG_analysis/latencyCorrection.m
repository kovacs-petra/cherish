function [EEG] = latencyCorrection(EEGin)
% Previous usage: [EEG] = latencyCorrection(EEG,eventList,subjID)
% No latencies added at
% - start/end of block
% - The datafiles with the reconstructed triggers

% StimTrak used for which subjects?
StimTrak_logical = [0 0 0 0 0 0 0 0 0 0 1 0 1 1 1 0 1 1 1 1 1 1 0 1 1 1]; 
current_subject = str2double(EEGin.filename(7:8));
% current_subject = 3;

% % Adjust event list: Convert to binary and take only last 6 bits
% tmp = {EEGin.event.type}; % Initial char array representing the binary
% tmp = cellfun(@(x) x(2:end),tmp,'UniformOutput',false); % Reject first entry ("T")
% tmp = str2double(tmp);
% 
% tmp(isnan(tmp)) = '0';
% % tmpEventList = dec2bin(tmp);
% tmpEventList = tmp;
% 
% % Keep the binary stimuli in the EEG list
% for ii = 1:length(tmpEventList)
% 	EEGin.event(ii).type = tmpEventList(:,ii);
% end
% 
% remove events without event/trigger value
[~,triggerStart]=ismember({'T'},{EEGin.event.code}); % find where events start to exclude first boundary segment
idx = find(cell2mat({EEGin.event.type})==0); % find indices of nonEvents
noEvent = idx(triggerStart:end);

EEG = pop_editeventvals(EEGin, 'delete', noEvent); % delete nonEvents from the list (saved in EEG.event, not in EEG.urevent

for ii = 1:length(EEG.event)
	eventList{ii,:} = EEG.event(ii).type;
end

% responseCodes={'T198','T196','T202','T200'}; \\ {'looming hit','looming
% miss','receding hit','receding miss'}
% blockCodes = {'T100'}; % Code denoting start of block
% blockCodes = {100:110};
onsetCodes = [210]; % trigger code where StimTrak signal was also sent
changeCodes = [150:190,199]; % other sound codes to latency correct, but ST was not sent out

% Create logical arrays of onset and change indicators
onsetList=zeros(length(eventList),1);
changeList=zeros(length(eventList),1);
% idleIndicator=zeros(length(eventList),1);

[~,eventStart]=ismember("T",{EEG.event.code}); % Find where actual data starts

for ii=eventStart:length(eventList) %Start after the encoding of stimulus infos
    if ismember(cell2mat(eventList(ii,1)),onsetCodes)
    % if eventList{2,1}{1,1}(1,end-6)==0 % All stim start by "10": motion onset
        % if str2double(eventList(ii,end-5))==0 % Onset
            onsetList(ii)=1;
%         elseif str2double(eventList(ii,end-5))==1 % Change
%             changeList(ii)=1;
        % end
%     elseif str2double(eventList(ii,end-6))==1 && ~ismember(EEG.urevent(ii).type,responseCodes) % If no stim and no response
%         idleIndicator(ii)=1;
    elseif ismember(cell2mat(eventList(ii,1)),changeCodes)
        changeList(ii)=1;
    end
end

onsetList = logical(onsetList);
changeList = logical(changeList);

if StimTrak_logical(current_subject)
    [~,StimTrakChan]=ismember("StimTrak",{EEG.chanlocs.labels});
    simGradient = gradient(EEG.data(StimTrakChan,:)); % Calculate instantaneous gradient on StimTrak
end
eventOnsLats = [EEG.event(onsetList).latency]'; % Trigger onset latencies
eventCngLats = [EEG.event(changeList).latency]'; % Trigger change latencies

maxGrads = zeros(length(eventOnsLats),1); % Store samplePoint of max simGradient in trial , (nr. of onsets = nr of trials)

% logName=fullfile('D:\WP1c\data_processed',subjID,[subjID '_ICA_log']); % for ICA log
% % logName=fullfile('E:\WP1c\data_processed',subjID,[subjID '_ICA_log']); % for ICA log
% diary (logName)

for ii=1:length(eventOnsLats)
    % %     if ii~=length(eventOnsLats) 
    % %         try
    % %             time_idx = eventOnsLats(ii):(eventOnsLats(ii+1));
    % %
    % %             if length(time_idx) > 1
    % %                 time_idx = time_idx(1:end-1);
    % %             end
    % %
    % %             [~,minGrad]= min(simGradient(time_idx)); % -1 for last sample before next trial begins
    % %
    % %             assert(numel(maxGrads(ii)) == numel(max(simGradient(eventOnsLats(ii):(eventOnsLats(ii)+minGrad)))), "something bad")
    % %
    % %             [~,maxGrads(ii)]= max(simGradient(eventOnsLats(ii):(eventOnsLats(ii)+minGrad))); % Find max before middle (in case offset pulse gives bigger gradient)
    % %         catch
    % %             warning('problem')
    % %             maxGrads(ii) = 0;
    % %         end
    % %         % maxGrads is in samples from trigger latency
    % %     else
    % %         [~,minGrad]= min(simGradient(eventOnsLats(ii):end)); % -1 for last sample before next trial begins
    % %         [~,maxGrads(ii)]= max(simGradient(eventOnsLats(ii):(eventOnsLats(ii)+minGrad)));
    % %     end
%     try
if StimTrak_logical(current_subject)
        if ii~=length(eventOnsLats)
            [~,minGrad]= min(simGradient(eventOnsLats(ii):(eventOnsLats(ii+1)-1))); % -1 for last sample before next trial begins
            [~,maxGrads(ii)]= max(simGradient(eventOnsLats(ii):(eventOnsLats(ii)+minGrad))); % Find max before middle (in case offset pulse gives bigger gradient)
            % maxGrads is in samples from trigger latency
            % disp(ii)
        else
            [~,minGrad]= min(simGradient(eventOnsLats(ii):end)); % -1 for last sample before next trial begins
            [~,maxGrads(ii)]= max(simGradient(eventOnsLats(ii):(eventOnsLats(ii)+minGrad)));
        end
%     catch
%         ID = ['Trial ',num2str(ii)];
%         warning(ID)
%     end
end
end

% diary off

% If urevent has numeric values, convert them to char
cc = 1;
while cc <= length(EEG.urevent)
    if isnumeric(EEG.urevent(cc).type)
        EEG.urevent(cc).type = num2str(EEG.urevent(cc).type);
    end
    cc = cc+1;
end
% blockIdx = find(ismember({EEG.urevent.type}',blockCodes{1,1}));
% blockIdx = find(isequal({EEG.urevent.type}',blockCodes));
% blockIdx = [blockIdx;length(EEG.event)];
onsetIdx = find(onsetList);
changeIdx = find(changeList);

% Augments all instances in a trial with the trial-specific latency, except
% for triggers denoting start/end of a block and resting state.
cnt=1; % Counter for trial-latency correction
for ii = 1:length(onsetList)
    if ismember(ii,onsetIdx) 
        if StimTrak_logical(current_subject)
            EEG.event(ii).latency=EEG.event(ii).latency+maxGrads(cnt); % FIX maxGrads
            if ismember(ii+1,changeIdx)
                EEG.event(ii+1).latency=EEG.event(ii+1).latency+maxGrads(cnt); % FIX maxGrads
                if ismember(ii+2,changeIdx)
                    EEG.event(ii+2).latency=EEG.event(ii+2).latency+maxGrads(cnt); % FIX maxGrads
                end
            end
        else % No StimTrak was used for current subject
            lat = 59.76; % hard-coded latency value based on subject 11 and previous measurements
            EEG.event(ii).latency=EEG.event(ii).latency+lat; % FIX maxGrads
            if ismember(ii+1,changeIdx)
                EEG.event(ii+1).latency=EEG.event(ii+1).latency+lat; % FIX maxGrads
                if ismember(ii+2,changeIdx)
                    EEG.event(ii+2).latency=EEG.event(ii+2).latency+lat; % FIX maxGrads
                end
            end
        end
        cnt=cnt+1;
    end
end

% Clean up EEG.event boundary events
for c = 1:length(EEG.event)
    EEG.event(c).duration = 1;
end

% Save latency corrected events for later importing
outFile = strcat('02_latencyCorrected\cleaned_events_s',num2str(current_subject),'.txt');
eventTable = struct2table(EEG.event);
writetable(eventTable, outFile, 'Delimiter', '\t');

end


% Initial code, only augmenting onset and change latencies
% eventOnsLats = eventOnsLats + (maxGrads-1)';
% % eventCngLats = eventCngLats + maxGrads-1;
% 
% onsetIdx = find(onsetList);
% % changeIdx = find(changeList);
% 
% for ii=1:length(onsetIdx) % Substitute the corrected latencies in EEG.event structure
%     EEG.event(onsetIdx(ii)).latency = eventOnsLats(ii);
% %     EEG.event(changeIdx(ii)).latency = eventCngLats(ii);
% end
