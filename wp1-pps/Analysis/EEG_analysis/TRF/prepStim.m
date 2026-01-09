function [stim_loom,stim_rec,stim_pps,stim_eps,...
            resp_loom,resp_rec,resp_pps,resp_eps] = ...
    prepStim(EEG,sub,stim_loom,stim_rec,stim_pps,stim_eps,...
            resp_loom,resp_rec,resp_pps,resp_eps)

aDim             = 1; % dimension in which azi data is stored
dDim             = 2; % dimension in which distance data is stored

% Fill up stim.data with the stimulus position coordinates
% aziData = repmat(1e10,size(EEG.data,2),1);
% dData = aziData;

% Set event triggers
onsetCode = 210;
loom = 151:154;
rec = 155:158;
pps = 159:162;
eps = 163:166;

% Counter for trials in specific conditions
loomTrial = 1;
recTrial = 1;
ppsTrial = 1;
epsTrial = 1;

% Set necessary values
fs = EEG.srate;
durMov = 2*fs; % in samples
PPS = 0.2; % in m
EPS = 2; % in m

% Exclude eye channels
noEyeChans = [1:4,6:9,11:20,22:26,28:size(EEG.data,1)];
noEyeLocs = EEG.chanlocs(1,noEyeChans);

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
        if ismember(nextType,loom)
            %%% Radius: 2 m to 0.2 m %%%
            %   Stationary part: 2 m to 2 m
            % dData(eventStartSamp:(eventStartSamp+durStatOnset)-1,1) = ...
            %     linspace(EPS,EPS,durStatOnset);
            stim_loom.data{sub,loomTrial}(1:durStatOnset,dDim) = ...
                linspace(EPS,EPS,durStatOnset);

            %   Moving part: 2 m to 0.2 m
            % dData(eventStartSamp+durStatOnset:...
            %     (eventStartSamp+durStatOnset+durMov)-1,1) = ...
            %     linspace(EPS,PPS,durMov);
            stim_loom.data{sub,loomTrial}(durStatOnset+1:durStatOnset+durMov,dDim) = ...
                linspace(EPS,PPS,durMov);

            %%% Azimuth: radial %%%
            if not(mod(nextType,2)) % odd trigger
                azi = -90;
            else                    % even trigger
                azi = 90;
            end
        
            % aziData(eventStartSamp:(eventStartSamp+durStatOnset+...
            %     durMov)-1,1) = linspace(azi,azi,durStatOnset+durMov);
            stim_loom.data{sub,loomTrial}(1:durStatOnset+durMov,aDim) = ...
                linspace(azi,azi,durStatOnset+durMov);

            % Collect the EEG data for this trial in the resp structure
            resp_loom.data{sub,loomTrial} = EEG.data(noEyeChans,...
                eventStartSamp:eventStartSamp+durStatOnset+durMov-1);

            loomTrial = loomTrial+1;

        %%%%%%%
        % REC %
        %%%%%%%
        elseif ismember(nextType,rec)
            %%% Radius: 0.2 m to 2 m %%%
            %   Stationary part: 0.2 m to 0.2 m
            % dData(eventStartSamp:(eventStartSamp+durStatOnset)-1,1) = ...
            %     linspace(PPS,PPS,durStatOnset);
            stim_rec.data{sub,recTrial}(1:durStatOnset,dDim) = ...
                linspace(PPS,PPS,durStatOnset);

            %   Moving part: 0.2 m to 2 m
            % dData(eventStartSamp+durStatOnset:...
            %     (eventStartSamp+durStatOnset+durMov)-1,1) = ...
            %     linspace(PPS,EPS,durMov);
            stim_rec.data{sub,recTrial}(durStatOnset+1:durStatOnset+durMov,dDim) = ...
                linspace(PPS,EPS,durMov);

            %%% Azimuth: radial %%%
            if not(mod(nextType,2)) % odd trigger
                azi = -90;
            else                    % even trigger
                azi = 90;
            end
        
            % aziData(eventStartSamp:(eventStartSamp+durStatOnset+...
            %     durMov)-1,1) = linspace(azi,azi,durStatOnset+durMov);
            stim_rec.data{sub,recTrial}(1:durStatOnset+durMov,aDim) = ...
                linspace(azi,azi,durStatOnset+durMov);

            % Collect the EEG data for this trial in the resp structure
            resp_rec.data{sub,recTrial} = EEG.data(noEyeChans,...
                eventStartSamp:eventStartSamp+durStatOnset+durMov-1);

            recTrial = recTrial+1;

        %%%%%%%
        % PPS %
        %%%%%%%
        elseif ismember(nextType,pps)
            %%% Radius: 0.2 m to 0.2 m %%%
            %   Stationary part: 0.2 m to 0.2 m
            % dData(eventStartSamp:(eventStartSamp+durStatOnset)-1,1) = ...
            %     linspace(PPS,PPS,durStatOnset);
            stim_pps.data{sub,ppsTrial}(1:durStatOnset,dDim) = ...
                linspace(PPS,PPS,durStatOnset);

            %   Moving part: 0.2 m to 0.2 m
            % dData(eventStartSamp+durStatOnset:...
            %     (eventStartSamp+durStatOnset+durMov)-1,1) = ...
            %     linspace(PPS,PPS,durMov);
            stim_pps.data{sub,ppsTrial}(durStatOnset+1:durStatOnset+durMov,dDim) = ...
                linspace(PPS,PPS,durMov);

            %%% Azimuth: angular %%%
            if not(mod(nextType,2)) % odd trigger
                azi = -90;
            else                    % even trigger
                azi = 90;
            end
        
            % aziData(eventStartSamp:(eventStartSamp+durStatOnset)-1,1) = ...
            %     linspace(azi,azi,durStatOnset);
            stim_pps.data{sub,ppsTrial}(1:durStatOnset,aDim) = ...
                linspace(azi,azi,durStatOnset);
            % aziData(eventStartSamp+durStatOnset:(eventStartSamp+durStatOnset+...
            %     durMov)-1,1) = linspace(-azi,azi,durMov);
            stim_pps.data{sub,ppsTrial}(durStatOnset+1:durStatOnset+durMov,aDim) = ...
                linspace(azi,-azi,durMov);

            % Collect the EEG data for this trial in the resp structure
            resp_pps.data{sub,ppsTrial} = EEG.data(noEyeChans,...
                eventStartSamp:eventStartSamp+durStatOnset+durMov-1);

            ppsTrial = ppsTrial+1;

        %%%%%%%
        % EPS %
        %%%%%%%
        elseif ismember(nextType,eps)
            %%% Radius: 2 m to 2 m %%%
                % Stationary part: 2 m to 2 m
            % dData(eventStartSamp:(eventStartSamp+durStatOnset)-1,1) = ...
            %     linspace(EPS,EPS,durStatOnset);
            stim_eps.data{sub,epsTrial}(1:durStatOnset,dDim) = ...
                linspace(EPS,EPS,durStatOnset);

                % Moving part: 0.2 m to 0.2 m
            % dData(eventStartSamp+durStatOnset:...
            %     (eventStartSamp+durStatOnset+durMov)-1,1) = ...
            %     linspace(EPS,EPS,durMov);
            stim_eps.data{sub,epsTrial}(durStatOnset+1:durStatOnset+durMov,dDim) = ...
                linspace(EPS,EPS,durMov);

            %%% Azimuth: angular %%%
            if not(mod(nextType,2)) % odd trigger
                azi = -90;
            else                    % even trigger
                azi = 90;
            end
        
            % aziData(eventStartSamp:(eventStartSamp+durStatOnset)-1,1) = ...
            %     linspace(azi,azi,durStatOnset);
            stim_eps.data{sub,epsTrial}(1:durStatOnset,aDim) = ...
                linspace(azi,azi,durStatOnset);
            % aziData(eventStartSamp+durStatOnset:(eventStartSamp+...
            %     +durStatOnset+durMov)-1,1) = linspace(-azi,azi,durMov);
            stim_eps.data{sub,epsTrial}(durStatOnset+1:durStatOnset+durMov,aDim) = ...
                linspace(azi,-azi,durMov);

            % Collect the EEG data for this trial in the resp structure
            resp_eps.data{sub,epsTrial} = EEG.data(noEyeChans,...
                eventStartSamp:eventStartSamp+durStatOnset+durMov-1);

            epsTrial = epsTrial+1;

        end

    end

    i = i + 1;
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
