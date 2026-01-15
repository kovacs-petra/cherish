function [loom,rec,pps,eps...
            ...,resp_loom,resp_rec,resp_pps,resp_eps...
            ] = ...
    prepStim(EEG,NTrialsPerCond,loom,rec,pps,eps,dDim,aDim...
            ...,resp_loom,resp_rec,resp_pps,resp_eps
            )

NVar = length(loom.names);
% dDim = 1; % col in which distance data is stored
% aDim = 2; % col in which azi data is stored

% Init stimulus position coordinates
% aziData = repmat(1e10,size(EEG.data,2),1);
% dData = aziData;

% Set event triggers
onsetCode = 210;
tr_loom = 151:154;
tr_rec = 155:158;
tr_pps = 159:162;
tr_eps = 163:166;

% Counter for trials in specific conditions
loomTrial = 1; loomTestTrial = randi(NTrialsPerCond/2,1);
recTrial = 1; recTestTrial = randi(NTrialsPerCond/2,1);
ppsTrial = 1; ppsTestTrial = randi(NTrialsPerCond/2,1);
epsTrial = 1; epsTestTrial = randi(NTrialsPerCond/2,1);

% Exclude eye channels
noEyeChans = [1:4,6:9,11:20,22:26,28:size(EEG.data,1)];
NChan = length(noEyeChans);
noEyeLocs = EEG.chanlocs(1,noEyeChans);
loom.chanlocs = noEyeLocs;
rec.chanlocs  = noEyeLocs;
eps.chanlocs  = noEyeLocs;
pps.chanlocs  = noEyeLocs;

% Init training and test sets
nfold = NTrialsPerCond; % how many partitions to make
ntest = 1; % how many of these partitions to hold out for testing

loom.strain = cell(nfold-ntest,1); loom.stest  = zeros(1,NVar);
rec.strain  = cell(nfold-ntest,1); rec.stest   = zeros(1,NVar);
pps.strain  = cell(nfold-ntest,1); pps.stest   = zeros(1,NVar);
eps.strain  = cell(nfold-ntest,1); eps.stest   = zeros(1,NVar);

loom.rtrain = cell(nfold-ntest,1); loom.rtest  = zeros(1,NChan);
rec.rtrain  = cell(nfold-ntest,1); rec.rtest   = zeros(1,NChan);
pps.rtrain  = cell(nfold-ntest,1); pps.rtest   = zeros(1,NChan);
eps.rtrain  = cell(nfold-ntest,1); eps.rtest   = zeros(1,NChan);

% Set necessary values
fs = EEG.srate;
durMov = 2*fs; % in samples
PPS = 0.2; % in m
EPS = 2; % in m

i = 1;
while i <= length(EEG.event) && i+1 <= length(EEG.event)
    if ischar(EEG.event(i).type)
        trigger = str2double(EEG.event(i).type);
    else
        trigger = EEG.event(i).type;
    end

    if ischar(EEG.event(i+1).type)
        nextType = str2double(EEG.event(i+1).type);
    else
        nextType = EEG.event(i+1).type;
    end

    if trigger == onsetCode
        durStatOnset = round(EEG.event(i+1).latency-EEG.event(i).latency,-2); % in ms (samples)
        eventStartSamp = round(EEG.event(i).latency);
        
        %%%%%%%%
        % LOOM %
        %%%%%%%%
        if ismember(nextType,tr_loom)            
            % Decode offset azi from the trigger
            if not(mod(nextType,2)) % odd trigger
                azi = -90;
            else                    % even trigger
                azi = 90;
            end

            % Record stimulus and response data 
            if loomTrial  == loomTestTrial % Random trial goes into the test set
                % % Stationary part % %
                % Radius:
                loom.stest(1:durStatOnset,dDim) = linspace(EPS,EPS,durStatOnset);
                % Azi:
                loom.stest(1:durStatOnset,aDim) = linspace(azi,azi,durStatOnset);
                
                % % Moving part % %
                % Radius:
                loom.stest(durStatOnset+1:durStatOnset+durMov,dDim) = ...
                    linspace(EPS,PPS,durMov);
                % Azi:
                loom.stest(durStatOnset+1:durStatOnset+durMov,aDim) = ...
                    linspace(azi,azi,durMov);

                % % EEG data for this trial % % 
                loom.rtest(1:durStatOnset+durMov,:) = EEG.data(noEyeChans,...
                eventStartSamp:eventStartSamp+durStatOnset+durMov-1)';

            else % Remaining trials go into the training set
                % % Stationary part % %
                % Radius:
                loom.strain{loomTrial,1}(1:durStatOnset,dDim) = linspace(EPS,EPS,durStatOnset);
                % Azi:
                loom.strain{loomTrial,1}(1:durStatOnset,aDim) = linspace(azi,azi,durStatOnset);

                % % Moving part % %
                % Radius:
                loom.strain{loomTrial,1}(durStatOnset+1:durStatOnset+durMov,dDim) = ...
                    linspace(EPS,PPS,durMov);
                % Azi:
                loom.strain{loomTrial,1}(durStatOnset+1:durStatOnset+durMov,aDim) = ...
                    linspace(azi,azi,durMov);

                % % EEG data for this trial % % 
                loom.rtrain{loomTrial,1}(1:durStatOnset+durMov,:) = EEG.data(noEyeChans,...
                eventStartSamp:eventStartSamp+durStatOnset+durMov-1)';
            end

            loomTrial = loomTrial+1;

        %%%%%%%
        % REC %
        %%%%%%%
        elseif ismember(nextType,tr_rec)
            % Decode offset azi from the trigger
            if not(mod(nextType,2)) % odd trigger
                azi = -90;
            else                    % even trigger
                azi = 90;
            end

            % Record stimulus and response data 
            if recTrial  == recTestTrial % Random trial goes into the test set
                % % Stationary part % %
                % Radius:
                rec.stest(1:durStatOnset,dDim) = linspace(PPS,PPS,durStatOnset);
                % Azi:
                rec.stest(1:durStatOnset,aDim) = linspace(azi,azi,durStatOnset);
                
                % % Moving part % %
                % Radius:
                rec.stest(durStatOnset+1:durStatOnset+durMov,dDim) = ...
                    linspace(PPS,EPS,durMov);
                % Azi:
                rec.stest(durStatOnset+1:durStatOnset+durMov,aDim) = ...
                    linspace(azi,azi,durMov);

                % % EEG data for this trial % % 
                rec.rtest(1:durStatOnset+durMov,:) = EEG.data(noEyeChans,...
                eventStartSamp:eventStartSamp+durStatOnset+durMov-1)';

            else % Remaining trials go into the training set
                % % Stationary part % %
                % Radius:
                rec.strain{recTrial,1}(1:durStatOnset,dDim) = linspace(PPS,PPS,durStatOnset);
                % Azi:
                rec.strain{recTrial,1}(1:durStatOnset,aDim) = linspace(azi,azi,durStatOnset);

                % % Moving part % %
                % Radius:
                rec.strain{recTrial,1}(durStatOnset+1:durStatOnset+durMov,dDim) = ...
                    linspace(PPS,EPS,durMov);
                % Azi:
                rec.strain{recTrial,1}(durStatOnset+1:durStatOnset+durMov,aDim) = ...
                    linspace(azi,azi,durMov);

                % % EEG data for this trial % % 
                rec.rtrain{recTrial,1}(1:durStatOnset+durMov,:) = EEG.data(noEyeChans,...
                eventStartSamp:eventStartSamp+durStatOnset+durMov-1)';
            end

            recTrial = recTrial+1;

        %%%%%%%
        % PPS %
        %%%%%%%
        elseif ismember(nextType,tr_pps)
            % Decode offset azi from the trigger
            if not(mod(nextType,2)) % odd trigger
                azi = -90;
            else                    % even trigger
                azi = 90;
            end

            % Record stimulus and response data 
            if ppsTrial == ppsTestTrial % Random trial goes into the test set
                % % Stationary part % %
                % Radius:
                pps.stest(1:durStatOnset,dDim) = linspace(PPS,PPS,durStatOnset);
                % Azi:
                pps.stest(1:durStatOnset,aDim) = linspace(-azi,-azi,durStatOnset);
                
                % % Moving part % %
                % Radius:
                pps.stest(durStatOnset+1:durStatOnset+durMov,dDim) = ...
                    linspace(PPS,PPS,durMov);
                % Azi:
                pps.stest(durStatOnset+1:durStatOnset+durMov,aDim) = ...
                    linspace(-azi,azi,durMov);

                % % EEG data for this trial % % 
                pps.rtest(1:durStatOnset+durMov,:) = EEG.data(noEyeChans,...
                eventStartSamp:eventStartSamp+durStatOnset+durMov-1)';

            else % Remaining trials go into the training set
                % % Stationary part % %
                % Radius:
                pps.strain{ppsTrial,1}(1:durStatOnset,dDim) = linspace(PPS,PPS,durStatOnset);
                % Azi:
                pps.strain{ppsTrial,1}(1:durStatOnset,aDim) = linspace(-azi,-azi,durStatOnset);

                % % Moving part % %
                % Radius:
                pps.strain{ppsTrial,1}(durStatOnset+1:durStatOnset+durMov,dDim) = ...
                    linspace(PPS,PPS,durMov);
                % Azi:
                pps.strain{ppsTrial,1}(durStatOnset+1:durStatOnset+durMov,aDim) = ...
                    linspace(-azi,azi,durMov);

                % % EEG data for this trial % % 
                pps.rtrain{ppsTrial,1}(1:durStatOnset+durMov,:) = EEG.data(noEyeChans,...
                eventStartSamp:eventStartSamp+durStatOnset+durMov-1)';
            end

            ppsTrial = ppsTrial+1;

        %%%%%%%
        % EPS %
        %%%%%%%
        elseif ismember(nextType,tr_eps)
            % Decode offset azi from the trigger
            if not(mod(nextType,2)) % odd trigger
                azi = -90;
            else                    % even trigger
                azi = 90;
            end

            % Record stimulus and response data 
            if epsTrial  == epsTestTrial % Random trial goes into the test set
                % % Stationary part % %
                % Radius:
                eps.stest(1:durStatOnset,dDim) = linspace(EPS,EPS,durStatOnset);
                % Azi:
                eps.stest(1:durStatOnset,aDim) = linspace(-azi,-azi,durStatOnset);
                
                % % Moving part % %
                % Radius:
                eps.stest(durStatOnset+1:durStatOnset+durMov,dDim) = ...
                    linspace(EPS,EPS,durMov);
                % Azi:
                eps.stest(durStatOnset+1:durStatOnset+durMov,aDim) = ...
                    linspace(-azi,azi,durMov);

                % % EEG data for this trial % % 
                eps.rtest(1:durStatOnset+durMov,:) = EEG.data(noEyeChans,...
                eventStartSamp:eventStartSamp+durStatOnset+durMov-1)';

            else % Remaining trials go into the training set
                % % Stationary part % %
                % Radius:
                eps.strain{epsTrial,1}(1:durStatOnset,dDim) = linspace(EPS,EPS,durStatOnset);
                % Azi:
                eps.strain{epsTrial,1}(1:durStatOnset,aDim) = linspace(-azi,-azi,durStatOnset);

                % % Moving part % %
                % Radius:
                eps.strain{epsTrial,1}(durStatOnset+1:durStatOnset+durMov,dDim) = ...
                    linspace(EPS,EPS,durMov);
                % Azi:
                eps.strain{epsTrial,1}(durStatOnset+1:durStatOnset+durMov,aDim) = ...
                    linspace(-azi,azi,durMov);

                % % EEG data for this trial % % 
                eps.rtrain{epsTrial,1}(1:durStatOnset+durMov,:) = EEG.data(noEyeChans,...
                eventStartSamp:eventStartSamp+durStatOnset+durMov-1)';
            end

            epsTrial = epsTrial+1;

        end

    end

    i = i + 1;
end

% If there's empty elements (missing trials) in strain (and rtrain), get
% rid of these elements
for idx = 1:nfold-ntest
    if idx <= length(loom.strain) && isempty(loom.strain{idx,1})
        loom.strain = loom.strain([1:idx-1,idx+1:end],1);
        loom.rtrain = loom.rtrain([1:idx-1,idx+1:end],1);
    end
    if idx <= length(rec.strain) && isempty(rec.strain{idx,1})
        rec.strain = rec.strain([1:idx-1,idx+1:end,[]],1);
        rec.rtrain = rec.rtrain([1:idx-1,idx+1:end,[]],1);
    end
    if idx <= length(pps.strain) && isempty(pps.strain{idx,1})
        pps.strain = pps.strain([1:idx-1,idx+1:end,[]],1);
        pps.rtrain = pps.rtrain([1:idx-1,idx+1:end,[]],1);
    end
    if idx <= length(eps.strain) && isempty(eps.strain{idx,1})
        eps.strain = eps.strain([1:idx-1,idx+1:end,[]],1);
        eps.rtrain = eps.rtrain([1:idx-1,idx+1:end,[]],1);
    end
end

% Normalize the EEG data
loom.rtest = loom.rtest/std(loom.rtest(:));
rec.rtest = rec.rtest/std(rec.rtest(:));
pps.rtest = pps.rtest/std(pps.rtest(:));
eps.rtest = eps.rtest/std(eps.rtest(:));

for i = 1:length(loom.rtrain)
    loom.rtrain{i,1} = loom.rtrain{i,1}/std(loom.rtrain{i,1}(:));
end
for i = 1:length(rec.rtrain)
    rec.rtrain{i,1} = rec.rtrain{i,1}/std(rec.rtrain{i,1}(:));
end
for i = 1:length(pps.rtrain)
    pps.rtrain{i,1} = pps.rtrain{i,1}/std(pps.rtrain{i,1}(:));
end
for i = 1:length(eps.rtrain)
    eps.rtrain{i,1} = eps.rtrain{i,1}/std(eps.rtrain{i,1}(:));
end

% Resample from 1000 Hz to 100 Hz to reduce computation time
% aziData_res = resample(aziData,1,10);
% dData_res = resample(dData,1,10);
% stim.data{1,:} = dData_res;
% stim.data{2,:} = aziData_res;
% stim.fs = 0.1*stim.fs;

% Turn nans back to inf if needed
% for i = 1:length(stim.data{1})
%     if isnan(stim.data{1}(i,1)) 
%         stim.data{1}(i,1) = Inf;
%     end
% 
%     if isnan(stim.data{1}(i,2)) 
%         stim.data{1}(i,2) = Inf;
%     end
% end
