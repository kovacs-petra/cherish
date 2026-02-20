function exportCleanedEventsToTxt(EEG, subNum)

outFile = strcat('01_triggerCorrected\cleaned_events_s',num2str(subNum),'.txt');

% Load subject's logVar
logfile = strcat('\\kfs\fileserver\Projektdaten\CherISH\data\wp-1\Behav\subject', ...
    num2str(subNum), '\sub', num2str(subNum), 'Log.mat');
load(logfile, "logVar");
logIndex = 2;

% Adjust event list
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

% Create struct, fill first row with boundary event
cleanedEvents = struct;
cleanedEvents.latency = EEG.event(1).latency;
cleanedEvents.duration = EEG.event(1).duration;
cleanedEvents.channel = EEG.event(1).channel;
cleanedEvents.bvtime = EEG.event(1).bvtime;
cleanedEvents.bvmknum = EEG.event(1).bvmknum;
cleanedEvents.visible = EEG.event(1).visible;
cleanedEvents.type = EEG.event(1).type;
cleanedEvents.code = EEG.event(1).code;
cleanedEvents.urevent = EEG.event(1).urevent;

i = 2;

while i <= length(EEG.event)
    
    thisType = EEG.event(i).type;
    if isnan(thisType)
        i = i + 1;
        break;
    end

    if thisType == 210 && logIndex <= size(logVar,1)
        % Get latency of the current event
        % latency210 = EEG.event(i).latency;
        cleanedEvents(end+1).latency = EEG.event(i).latency;
        cleanedEvents(end).duration = EEG.event(i).duration;
        cleanedEvents(end).channel = EEG.event(i).channel;
        cleanedEvents(end).bvtime = EEG.event(i).bvtime;
        cleanedEvents(end).bvmknum = EEG.event(i).bvmknum;
        cleanedEvents(end).visible = EEG.event(i).visible;
        cleanedEvents(end).type = EEG.event(i).type;
        cleanedEvents(end).code = EEG.event(i).code;
        cleanedEvents(end).urevent = EEG.event(i).urevent;

        % Next event (usually the motion trigger) is corrected from logVar
        if i + 1 <= length(EEG.event)
            correctType = logVar{logIndex, 20};
            latencyMotion = EEG.event(i+1).latency;

            cleanedEvents(end+1).latency = EEG.event(i+1).latency;
            cleanedEvents(end).duration = EEG.event(i+1).duration;
            cleanedEvents(end).channel = EEG.event(i+1).channel;
            cleanedEvents(end).bvtime = EEG.event(i+1).bvtime;
            cleanedEvents(end).bvmknum = EEG.event(i+1).bvmknum;
            cleanedEvents(end).visible = EEG.event(i+1).visible;
            cleanedEvents(end).type = correctType;
            cleanedEvents(end).code = EEG.event(i+1).code;
            cleanedEvents(end).urevent = EEG.event(i+1).urevent;

            % If no target, add GSR trigger (type + 20)
            if logVar{logIndex, 12} == 0
                cleanedEvents(end+1).latency = latencyMotion + 1;  % Add 1 sample offset
                cleanedEvents(end).duration = EEG.event(i+1).duration;
                cleanedEvents(end).channel = EEG.event(i+1).channel;
                cleanedEvents(end).bvtime = EEG.event(i+1).bvtime;
                cleanedEvents(end).bvmknum = EEG.event(i+1).bvmknum;
                cleanedEvents(end).visible = EEG.event(i+1).visible;
                cleanedEvents(end).type = correctType + 20;
                cleanedEvents(end).code = EEG.event(i+1).code;
                cleanedEvents(end).urevent = EEG.event(i+1).urevent;
            end

            logIndex = logIndex + 1;
        end

        % Skip bad motion triggers
        j = i + 2;
        while j <= length(EEG.event)
            t = EEG.event(j).type;
            if isnan(t)
                j = j + 1;
                continue;
            end

            if t >= 151 && t <= 186
                j = j + 1;
            elseif t < 140 || t >= 200
                break;
            else
                break;
            end
        end
        i = j;

    elseif thisType == 200
        cleanedEvents(end+1).latency = EEG.event(i).latency;
        cleanedEvents(end).type = 200;
        cleanedEvents(end).duration = EEG.event(i).duration;
        cleanedEvents(end).channel = EEG.event(i).channel;
        cleanedEvents(end).bvtime = EEG.event(i).bvtime;
        cleanedEvents(end).bvmknum = EEG.event(i).bvmknum;
        cleanedEvents(end).visible = EEG.event(i).visible;
        cleanedEvents(end).code = EEG.event(i).code;
        cleanedEvents(end).urevent = EEG.event(i).urevent;

        % Keep the 210 only if it immediately follows
        if i + 1 <= length(EEG.event)
            nextType = EEG.event(i+1).type;
            if nextType == 210
                % Let next loop handle it
            else
                i = i + 1;  % Skip the wrong event
            end
        end
        i = i + 1;

    elseif ismember(thisType, [48, 100:110, 190, 220])
        % Keep other important events
        cleanedEvents(end+1).latency = EEG.event(i).latency;
        cleanedEvents(end).type = thisType;
        cleanedEvents(end).duration = EEG.event(i).duration;
        cleanedEvents(end).channel = EEG.event(i).channel;
        cleanedEvents(end).bvtime = EEG.event(i).bvtime;
        cleanedEvents(end).bvmknum = EEG.event(i).bvmknum;
        cleanedEvents(end).visible = EEG.event(i).visible;
        cleanedEvents(end).code = EEG.event(i).code;
        cleanedEvents(end).urevent = EEG.event(i).urevent;
        i = i + 1;

    else
        i = i + 1;
    end
end

% Convert to table
% latencies = [cleanedEvents.latency]';
% types = [cleanedEvents.type]';
% eventTable = table(latencies, types, 'VariableNames', {'latency', 'type'});
eventTable = struct2table(cleanedEvents);

% Write to .txt file
writetable(eventTable, outFile, 'Delimiter', '\t');

fprintf('Exported %d cleaned events to %s\n', height(eventTable), outFile);
type = [cleanedEvents.type].';
loom = sum(histc(type,[151:154]));
rec = sum(histc(type,[155:158]));
pps = sum(histc(type,[159:162]));
eps = sum(histc(type,[163:166]));
fprintf('Looming: %d\nReceding: %d\nPPS: %d\nEPS: %d\n',loom,rec,pps,eps);

% EEG = pop_importevent( EEG, 'append','no','event', outFile, ...
%     'fields',{'latency','duration','channel','bvtime','bvmknum','visible','type','code', ...
%     'urevent'},'skipline',2,'timeunit',0.001,'align',1);

end
