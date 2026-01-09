function [stim_loom,stim_rec,stim_pps,stim_eps,...
         resp_loom,resp_rec,resp_pps,resp_eps] = prepTRF
% Usage: [stim_loom,stim_rec,stim_pps,stim_eps,...
%         resp_loom,resp_rec,resp_pps,resp_eps] = prepTRF
%
% Prepares stimuli and EEG files for analysis in the mTRF toolbox.
%
% Outputs:
%   stim        struct with fields of metadata and a "data" cell.
%               stim.data has one row for each stimulus feature (e.g.
%               distance, loudness) and one column for each subject. Each
%               cell is a column vector with NSamples rows.
%   resp        struct with fields of metadata and a "data" cell.
%               resp.data has one row and NSub columns. Each cell has
%               size NSamples X NChannels

% #Author: Petra Kovacs, 2025-2026, Austrian Academy of Sciences

% Set number of observations
NSub = 22;
allSub = 3:24;
NTrialsPerCond = 200;
% NConds = 4;

% Init model variable
% model_azi = cell(1,NSub);
% model_d = model_azi;

% Manage paths
dir_name = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\data\wp-1\EEG\03_filtered\';
path_eeglab = '\\kfs.oeaw.ac.at\Fileserver\Projektdaten\CherISH\tools\eeglab2025.0.0';
path_trf = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\tools\mTRF-Toolbox\mtrf';
dir_file = dir(dir_name);
addpath(dir_name,path_eeglab,path_trf);
[ALLEEG, EEG, CURRENTSET] = eeglab;

% Init stimulus structures
stim_loom        = struct;
stim_loom.names  = {'distance','azimuth'};
stim_loom.data   = cell(NSub,NTrialsPerCond);

stim_rec         = struct;
stim_rec.names   = {'distance','azimuth'};
stim_rec.data    = cell(NSub,NTrialsPerCond);

stim_pps         = struct;
stim_pps.names   = {'distance','azimuth'};
stim_pps.data    = cell(NSub,NTrialsPerCond);

stim_eps         = struct;
stim_eps.names   = {'distance','azimuth'};
stim_eps.data    = cell(NSub,NTrialsPerCond);

% Init response (EEG) structures
resp_loom                     = struct;
resp_loom.dataType            = 'EEG';
resp_loom.deviceName          = 'BrainVision';
resp_loom.data                = cell(NSub,NTrialsPerCond);

resp_rec                     = struct;
resp_rec.dataType            = 'EEG';
resp_rec.deviceName          = 'BrainVision';
resp_rec.data                = cell(NSub,NTrialsPerCond);
resp_rec.chanlocs            = resp_rec.data;

resp_pps                     = struct;
resp_pps.dataType            = 'EEG';
resp_pps.deviceName          = 'BrainVision';
resp_pps.data                = cell(NSub,NTrialsPerCond);
resp_pps.chanlocs            = resp_pps.data;

resp_eps                     = struct;
resp_eps.dataType            = 'EEG';
resp_eps.deviceName          = 'BrainVision';
resp_eps.data                = cell(NSub,NTrialsPerCond);
resp_eps.chanlocs            = resp_eps.data;

for i = 3:size(dir_file,1)
    filename = dir_file(i).name;
    sub = str2double(filename(7:8));

    if ismember(sub,allSub) && contains(filename,'.set')

        sub = sub-2;

        % Load EEG data
        EEG = pop_loadset('filename',filename,'filepath', dir_name);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = eeg_checkset( EEG );

        % Resample EEG
        EEG = pop_resample(EEG,100);
        EEG.data = double(EEG.data);

        % Create stimuli and response structures
        [stim_loom,stim_rec,stim_pps,stim_eps,...
         resp_loom,resp_rec,resp_pps,resp_eps] = ...
            prepStim(EEG,sub,stim_loom,stim_rec,stim_pps,stim_eps,...
                    resp_loom,resp_rec,resp_pps,resp_eps);

        % Save sampling rate and channel locations
        stim_loom.fs    = EEG.srate;
        stim_rec.fs     = EEG.srate;
        stim_pps.fs     = EEG.srate;
        stim_eps.fs     = EEG.srate;

        % Clear memory
        ALLEEG = []; EEG = []; CURRENTSET = [];
        
    end
end
save('stim_resp_data_20260108.mat','stim_loom','stim_rec','stim_pps','stim_eps',...
    'resp_loom','resp_rec','resp_pps','resp_eps');

% % Run the analysis
%         model_d{sub-2} = mTRFtrain(stim.data{1,:},resp.data,stim.fs,1,-100,1000,0.1);
%         model_azi{sub-2} = mTRFtrain(stim.data{2,:},resp.data,stim.fs,1,-100,1000,0.1);
% 
% % Average over subjects
% % model_azi = cell(1,NSub);
% % model_d = model_a;
% % for ss = 1:NSub
% %     if isstruct(model_azi{1,ss})
% %         dist_models{1,ss} = model_azi{1,ss}.w(1,:,:);
% %         azi_models{1,ss} = model_azi{1,ss}.w(2,:,:);
% %     end
% % end
% 
% mean_dist = mean([model_d{1,1}.w;model_d{1,2}.w;model_d{1,3}.w;model_d{1,4}.w;...
%     model_d{1,5}.w;model_d{1,6}.w;model_d{1,7}.w;model_d{1,8}.w;...
%     model_d{1,9}.w;model_d{1,10}.w;model_d{1,11}.w;model_d{1,12}.w;...
%     model_d{1,13}.w;model_d{1,14}.w;model_d{1,15}.w;model_d{1,16}.w;...
%     model_d{1,17}.w;model_d{1,18}.w;model_d{1,19}.w;model_d{1,20}.w;...
%     model_d{1,21}.w;model_d{1,22}.w]);
% 
% mean_azi = mean([model_azi{1,1}.w;model_azi{1,2}.w;model_azi{1,3}.w;model_azi{1,4}.w;...
%     model_azi{1,5}.w;model_azi{1,6}.w;model_azi{1,7}.w;model_azi{1,8}.w;...
%     model_azi{1,9}.w;model_azi{1,10}.w;model_azi{1,11}.w;model_azi{1,12}.w;...
%     model_azi{1,13}.w;model_azi{1,14}.w;model_azi{1,15}.w;model_azi{1,16}.w;...
%     model_azi{1,17}.w;model_azi{1,18}.w;model_azi{1,19}.w;model_azi{1,20}.w;...
%     model_azi{1,21}.w;model_azi{1,22}.w]);
% 
% m = struct;
% m.w = [mean_azi;mean_dist];
% m.b = model_azi{1,1}.b;
% m.t = model_azi{1,1}.t;
% m.fs = model_azi{1,1}.fs;
% m.Dir = model_azi{1,1}.Dir;
% m.type = model_azi{1,1}.type;
% 
% % Separate models for azi and distance
% m_a = struct;
% m_a.w = mean_azi;
% m_a.b = model_azi{1,1}.b;
% m_a.t = model_azi{1,1}.t;
% m_a.fs = model_azi{1,1}.fs;
% m_a.Dir = model_azi{1,1}.Dir;
% m_a.type = model_azi{1,1}.type;
% 
% m_d = struct;
% m_d.w = mean_dist;
% m_d.b = model_azi{1,1}.b;
% m_d.t = model_azi{1,1}.t;
% m_d.fs = model_azi{1,1}.fs;
% m_d.Dir = model_azi{1,1}.Dir;
% m_d.type = model_azi{1,1}.type;
% 
% % Plot azi model
% figure;
% subplot(1,2,1)
% mTRFplot(m_a,'trf');
% title('TRF weights, azimuth model')
% grid on
% 
% % Plot GFP
% subplot(1,2,2)
% mTRFplot(m_a,'gfp',[],'all');
% title('TRF Global Field Power, azimuth model')
% grid on
% 
% % Plot distance model
% figure;
% subplot(1,2,1)
% mTRFplot(m_d,'trf');
% title('TRF weights, distance model')
% grid on
% 
% % Plot GFP
% subplot(1,2,2)
% mTRFplot(m_d,'gfp',[],'all');
% title('TRF Global Field Power, distance model')
% grid on
% 
% % Plot multivariate model
% figure;
% subplot(1,2,1)
% mTRFplot(m,'mtrf',[],'all');
% title('TRF weights')
% grid on
% 
% % Plot GFP
% subplot(1,2,2)
% mTRFplot(m,'mgfp',[],'all');
% title('TRF Global Field Power')
% grid on

