function [loom,rec,pps,eps] = runTRF

% addpath '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\tools\mTRF-Toolbox\mtrf';

% Set number of observations
NSub = 22;
allSub = 3:24;
NTrialsPerCond = 200;

% Manage paths
dir_name = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\data\wp-1\EEG\07_interpolated\';
path_eeglab = '\\kfs.oeaw.ac.at\Fileserver\Projektdaten\CherISH\tools\eeglab2025.0.0';
path_trf = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\tools\mTRF-Toolbox\mtrf';
dir_file = dir(dir_name);
addpath(dir_name,path_eeglab,path_trf);
[ALLEEG, EEG, CURRENTSET] = eeglab;

% Init stimulus structures
loom        = struct;
loom.names  = {'distance','azimuth'};
% loom.data   = cell(NSub,1);

rec         = struct;
rec.names   = {'distance','azimuth'};
% rec.data    = cell(NSub,1);

pps         = struct;
pps.names   = {'distance','azimuth'};
% pps.data    = cell(NSub,1);

eps         = struct;
eps.names   = {'distance','azimuth'};
% eps.data    = cell(NSub,1);

% Init response (EEG) structures
% resp_loom                     = struct;
% resp_loom.dataType            = 'EEG';
% resp_loom.deviceName          = 'BrainVision';
% resp_loom.data                = cell(NSub,NTrialsPerCond);
% 
% resp_rec                     = struct;
% resp_rec.dataType            = 'EEG';
% resp_rec.deviceName          = 'BrainVision';
% resp_rec.data                = cell(NSub,1);
% resp_rec.chanlocs            = resp_rec.data;
% 
% resp_pps                     = struct;
% resp_pps.dataType            = 'EEG';
% resp_pps.deviceName          = 'BrainVision';
% resp_pps.data                = cell(NSub,1);
% resp_pps.chanlocs            = resp_pps.data;
% 
% resp_eps                     = struct;
% resp_eps.dataType            = 'EEG';
% resp_eps.deviceName          = 'BrainVision';
% resp_eps.data                = cell(NSub,1);
% resp_eps.chanlocs            = resp_eps.data;

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
        [loom,rec,pps,eps...
            ...,resp_loom,resp_rec,resp_pps,resp_eps...
            ] = ...
        prepStim(EEG,NTrialsPerCond,loom,rec,pps,eps...
            ...,resp_loom,resp_rec,resp_pps,resp_eps...
            );

        % Save sampling rate and channel locations
        loom.fs    = EEG.srate;
        rec.fs     = EEG.srate;
        pps.fs     = EEG.srate;
        eps.fs     = EEG.srate;

        % Clear memory
        ALLEEG = []; EEG = []; CURRENTSET = [];

        % Prepare stimulus and eeg structures if not already inputted
        % if ~nargin
        % [stim_loom,stim_rec,stim_pps,stim_eps,...
        %     resp_loom,resp_rec,resp_pps,resp_eps] = prepTRF;
        % end

        % Optimize the encoder's ability to predict EEG features from new stimulus
        % data: tune the regularization parameter using an efficient leave-one-out
        % cross-validation (CV) procedure
        fs = loom.fs;
        Dir = -1;
        tmin = 0;
        tmax = 500;
        lambda = 10.^(-6:2:6);
        loom.cv = mTRFcrossval(loom.strain,loom.rtrain,fs,Dir,tmin,tmax,lambda,'zeropad',0,'fast',1);
        rec.cv = mTRFcrossval(rec.strain,rec.rtrain,fs,Dir,tmin,tmax,lambda,'zeropad',0,'fast',1);
        pps.cv = mTRFcrossval(pps.strain,pps.rtrain,fs,Dir,tmin,tmax,lambda,'zeropad',0,'fast',1);
        eps.cv = mTRFcrossval(eps.strain,eps.rtrain,fs,Dir,tmin,tmax,lambda,'zeropad',0,'fast',1);

        % Get the optimal regularization value from the cv structure
        [loom.rmax,loom.idx] = max(mean(loom.cv.r,[1,3]));
        [rec.rmax,rec.idx]  = max(mean(rec.cv.r,[1,3]));
        [pps.rmax,pps.idx]  = max(mean(pps.cv.r,[1,3]));
        [eps.rmax,eps.idx]  = max(mean(eps.cv.r,[1,3]));

        % Train the model using this lambda value
        loom.model = mTRFtrain(loom.strain,loom.rtrain,fs,Dir,tmin,tmax,lambda(loom.idx),'zeropad',0);
        rec.model  = mTRFtrain(rec.strain,rec.rtrain,fs,Dir,tmin,tmax,lambda(rec.idx),'zeropad',0);
        pps.model  = mTRFtrain(pps.strain,pps.rtrain,fs,Dir,tmin,tmax,lambda(pps.idx),'zeropad',0);
        eps.model  = mTRFtrain(eps.strain,eps.rtrain,fs,Dir,tmin,tmax,lambda(eps.idx),'zeropad',0);

        % Test the model
        [loom.pred,loom.test] = mTRFpredict(loom.stest,loom.rtest,loom.model,'zeropad',0);
        [rec.pred,rec.test]   = mTRFpredict(rec.stest,rec.rtest,rec.model,'zeropad',0);
        [pps.pred,pps.test]   = mTRFpredict(pps.stest,pps.rtest,pps.model,'zeropad',0);
        [eps.pred,eps.test]   = mTRFpredict(eps.stest,eps.rtest,eps.model,'zeropad',0);

        % Save whatever needs to be saved (model, pred, test, anything
        % else?)
    end
end

% % Plotting the data % %
nfold = 200;

% Looming
figure
subplot(2,2,1), errorbar(1:numel(lambda),mean(loom.cv.r,[1,3]),std(loom.cv.r,0,[1,3])/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',-6:2:6), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')

% Plot CV error
subplot(2,2,2), errorbar(1:numel(lambda),mean(loom.cv.err,[1,3]),std(loom.cv.err,0,[1,3])/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',-6:2:6), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')

% Plot reconstruction
subplot(2,2,3), plot((1:length(loom.stest))/fs,loom.stest(:,2),'linewidth',2), hold on
plot((1:length(loom.pred(:,2)))/fs,mean(loom.pred(:,1),2),'linewidth',2), hold off, xlim([0,3]), axis square, grid on
title('Reconstruction'), xlabel('Time (s)'), ylabel('Amplitude (a.u.)'), legend('Orig','Pred')

% Plot test accuracy
subplot(2,2,4), bar(1,loom.rmax), hold on, bar(2,mean(loom.test.r)), hold off
set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')

% EPS
% Plot CV accuracy
figure
subplot(2,2,1), errorbar(1:numel(lambda),mean(eps.cv.r,[1,3]),std(eps.cv.r,0,[1,3])/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',-6:2:6), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')

% Plot CV error
subplot(2,2,2), errorbar(1:numel(lambda),mean(eps.cv.err,[1,3]),std(eps.cv.err,0,[1,3])/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',-6:2:6), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')

% Plot reconstruction
subplot(2,2,3), plot((1:length(eps.stest))/fs,eps.stest(:,1),'linewidth',2), hold on
plot((1:length(eps.pred))/fs,mean(eps.pred,2),'linewidth',2), hold off, xlim([0,3]), axis square, grid on
title('Reconstruction'), xlabel('Time (s)'), ylabel('Amplitude (a.u.)'), legend('Orig','Pred')

% Plot test accuracy
subplot(2,2,4), bar(1,eps.rmax), hold on, bar(2,mean(eps.test.r)), hold off
set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')

% Rec
% Plot CV accuracy
figure
subplot(2,2,1), errorbar(1:numel(lambda),mean(rec.cv.r,[1,3]),std(rec.cv.r,0,[1,3])/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',-6:2:6), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')

% Plot CV error
subplot(2,2,2), errorbar(1:numel(lambda),mean(rec.cv.err,[1,3]),std(rec.cv.err,0,[1,3])/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',-6:2:6), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')

% Plot reconstruction
subplot(2,2,3), plot((1:length(rec.stest))/fs,rec.stest(:,2),'linewidth',2), hold on
plot((1:length(rec.pred(:,2)))/fs,mean(rec.pred,2),'linewidth',2), hold off, xlim([0,3]), axis square, grid on
title('Reconstruction'), xlabel('Time (s)'), ylabel('Amplitude (a.u.)'), legend('Orig','Pred')

% Plot test accuracy
subplot(2,2,4), bar(1,rec.rmax), hold on, bar(2,mean(rec.test.r)), hold off
set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')

% figure;
% subplot(1,2,1)
% mTRFplot(model1000,'trf');
% title('TRF weights')
% grid on
%
% % Plot GFP
% subplot(1,2,2)
% mTRFplot(model1000,'gfp',[],'all');
% title('TRF Global Field Power')
% grid on
%
% figure; plot(model1000.w(:,:,60)');