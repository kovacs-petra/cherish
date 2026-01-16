function [loom,rec,pps,eps] = runTRF

% Set number of observations
NSub = 22;
% allSub = 3:24;
allSub = 3:13;
NTrialsPerCond = 200;
saveSub = [];

% Manage paths
dir_name = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\data\wp-1\EEG\07_interpolated\';
path_eeglab = '\\kfs.oeaw.ac.at\Fileserver\Projektdaten\CherISH\tools\eeglab2025.0.0';
path_trf = '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\tools\mTRF-Toolbox\mtrf';
dir_file = dir(dir_name);
addpath(dir_name,path_eeglab,path_trf);
[ALLEEG, EEG, CURRENTSET] = eeglab;

dDim = 1; % col in which distance data is stored
aDim = 2; % col in which azi data is stored

% Init stimulus structures
loom        = struct;
loom.names  = {'distance','azimuth'};

rec         = struct;
rec.names   = {'distance','azimuth'};
rec.allModels = struct;
    rec.allModels.w = cell(NSub,1);
    rec.allModels.b = cell(NSub,1);
    rec.allModels.t = cell(NSub,1);

pps         = struct;
pps.names   = {'distance','azimuth'};
pps.allModels = struct;
    pps.allModels.w = cell(NSub,1);
    pps.allModels.b = cell(NSub,1);
    pps.allModels.t = cell(NSub,1);

eps         = struct;
eps.names   = {'distance','azimuth'};
eps.allModels = struct;
    eps.allModels.w = cell(NSub,1);
    eps.allModels.b = cell(NSub,1);
    eps.allModels.t = cell(NSub,1);

allModels = struct;
    allModels.loom.w = cell(NSub,1);
    allModels.loom.b = cell(NSub,1);
    allModels.loom.t = cell(NSub,1);

    allModels.rec.w = cell(NSub,1);
    allModels.rec.b = cell(NSub,1);
    allModels.rec.t = cell(NSub,1);

    allModels.pps.w = cell(NSub,1);
    allModels.pps.b = cell(NSub,1);
    allModels.pps.t = cell(NSub,1);

    allModels.eps.w = cell(NSub,1);
    allModels.eps.b = cell(NSub,1);
    allModels.eps.t = cell(NSub,1);

% Init saving the data
loom_cvR    = []; rec_cvR    = []; pps_cvR    = []; eps_cvR    = [];
loom_cvErr  = []; rec_cvErr  = []; pps_cvErr  = []; eps_cvErr  = [];
loom_origD  = []; rec_origD  = []; pps_origA  = []; eps_origA  = [];
loom_predD  = []; rec_predD  = []; pps_predA  = []; eps_predA  = []; 
loom_rmax   = []; rec_rmax   = []; pps_rmax   = []; eps_rmax   = []; 
loom_testR  = []; rec_testR  = []; pps_testR  = []; eps_testR  = []; 
loom_trfWgt = []; rec_trfWgt = []; pps_trfWgt = []; eps_trfWgt = []; 

for i = 3:round(size(dir_file,1)/2) 
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
        [loom,rec,pps,eps,nfold...
            ...,resp_loom,resp_rec,resp_pps,resp_eps...
            ] = ...
        prepStim(EEG,NTrialsPerCond,loom,rec,pps,eps,dDim,aDim...
            ...,resp_loom,resp_rec,resp_pps,resp_eps...
            );

        % Save sampling rate and channel locations
        loom.fs    = EEG.srate;
        rec.fs     = EEG.srate;
        pps.fs     = EEG.srate;
        eps.fs     = EEG.srate;

        % Clear memory
        ALLEEG = []; EEG = []; CURRENTSET = [];

        % Optimize the encoder's ability to predict EEG features from new stimulus
        % data: tune the regularization parameter using an efficient leave-one-out
        % cross-validation (CV) procedure
        fs = loom.fs;
        Dir = -1;
        tmin = -100;
        tmax = 500;
        lambdavals = 0:2:10;
        % lambda = 10.^(-6:2:6);
        % lambda = 10.^(-2:2:12);
        lambda = 10.^lambdavals;
        loom.cv = mTRFcrossval(loom.strain,loom.rtrain,fs,Dir,tmin,tmax,lambda,'zeropad',1,'fast',0);
        rec.cv = mTRFcrossval(rec.strain,rec.rtrain,fs,Dir,tmin,tmax,lambda,'zeropad',1,'fast',0);
        pps.cv = mTRFcrossval(pps.strain,pps.rtrain,fs,Dir,tmin,tmax,lambda,'zeropad',1,'fast',0);
        eps.cv = mTRFcrossval(eps.strain,eps.rtrain,fs,Dir,tmin,tmax,lambda,'zeropad',1,'fast',0);

        % Get the optimal regularization value from the cv structure
        [loom.rmax,loom.idx] = max(mean(loom.cv.r(:,:,dDim)));
        [rec.rmax,rec.idx]  = max(mean(rec.cv.r(:,:,dDim)));
        [pps.rmax,pps.idx]  = max(mean(pps.cv.r(:,:,aDim)));
        [eps.rmax,eps.idx]  = max(mean(eps.cv.r(:,:,aDim)));

        % Train the model using this lambda value
        loom.model = mTRFtrain(loom.strain,loom.rtrain,fs,Dir,tmin,tmax,lambda(loom.idx),'zeropad',1);
        rec.model  = mTRFtrain(rec.strain,rec.rtrain,fs,Dir,tmin,tmax,lambda(rec.idx),'zeropad',1);
        pps.model  = mTRFtrain(pps.strain,pps.rtrain,fs,Dir,tmin,tmax,lambda(pps.idx),'zeropad',1);
        eps.model  = mTRFtrain(eps.strain,eps.rtrain,fs,Dir,tmin,tmax,lambda(eps.idx),'zeropad',1);

        % Test the model
        [loom.pred,loom.test] = mTRFpredict(loom.stest,loom.rtest,loom.model,'zeropad',1);
        [rec.pred,rec.test]   = mTRFpredict(rec.stest,rec.rtest,rec.model,'zeropad',1);
        [pps.pred,pps.test]   = mTRFpredict(pps.stest,pps.rtest,pps.model,'zeropad',1);
        [eps.pred,eps.test]   = mTRFpredict(eps.stest,eps.rtest,eps.model,'zeropad',1);

        % Save data
        % Cross-validation accuracy (Pearson's r)
        loom_cvR    = [loom_cvR; loom.cv.r];
        rec_cvR     = [rec_cvR; rec.cv.r]; 
        pps_cvR     = [pps_cvR; pps.cv.r]; 
        eps_cvR     = [eps_cvR; eps.cv.r];

        % Cross-validation error
        loom_cvErr  = [loom_cvErr; loom.cv.err]; 
        rec_cvErr   = [rec_cvErr; rec.cv.err]; 
        pps_cvErr   = [pps_cvErr; pps.cv.err]; 
        eps_cvErr   = [eps_cvErr; eps.cv.err];

        % Original spatial location (Distance or Azi)
        loom_origD  = [loom_origD, loom.stest(:,dDim)]; 
        rec_origD   = [rec_origD, rec.stest(:,dDim)]; 
        pps_origA   = [pps_origA, pps.stest(:,aDim)]; 
        eps_origA   = [eps_origA, eps.stest(:,aDim)];

        % Predicted spatial location (Distance or Azi)
        loom_predD  = [loom_predD, loom.pred(:,dDim)]; 
        rec_predD   = [rec_predD, rec.pred(:,dDim)]; 
        pps_predA   = [pps_predA, pps.pred(:,aDim)]; 
        eps_predA   = [eps_predA, eps.pred(:,aDim)];

        % Model performance: largest r on the training set
        loom_rmax   = [loom_rmax, loom.rmax]; 
        rec_rmax    = [rec_rmax, rec.rmax]; 
        pps_rmax    = [pps_rmax, pps.rmax]; 
        eps_rmax    = [eps_rmax, eps.rmax];

        % Model performance: r on the test set
        loom_testR  = [loom_testR, loom.test.r(:,dDim)]; 
        rec_testR   = [rec_testR, rec.test.r(:,dDim)]; 
        pps_testR   = [pps_testR, pps.test.r(:,aDim)]; 
        eps_testR   = [eps_testR, eps.test.r(:,aDim)];

        % TRF weights: what spatial location the brain represents
        w_chanavg_loom = mean(loom.model.w,2);
        w_chanavg_rec  = mean(rec.model.w,2);
        w_chanavg_pps  = mean(pps.model.w,2);
        w_chanavg_eps  = mean(eps.model.w,2);
        loom_trfWgt = [loom_trfWgt, w_chanavg_loom(:,dDim)]; 
        rec_trfWgt  = [rec_trfWgt, w_chanavg_rec(:,dDim)]; 
        pps_trfWgt  = [pps_trfWgt, w_chanavg_pps(:,aDim)]; 
        eps_trfWgt  = [eps_trfWgt, w_chanavg_eps(:,aDim)];

        % Save all models
        allModels.loom.w{sub,1} = loom.model.w;
        allModels.loom.b{sub,1} = loom.model.b;
        allModels.loom.t{sub,1} = loom.model.t;

        allModels.rec.w{sub,1} = rec.model.w;
        allModels.rec.b{sub,1} = rec.model.b;
        allModels.rec.t{sub,1} = rec.model.t;

        allModels.pps.w{sub,1} = pps.model.w;
        allModels.pps.b{sub,1} = pps.model.b;
        allModels.pps.t{sub,1} = pps.model.t;

        allModels.eps.w{sub,1} = eps.model.w;
        allModels.eps.b{sub,1} = eps.model.b;
        allModels.eps.t{sub,1} = eps.model.t;

        % Save subject numbers
        saveSub = [saveSub,sub];

        % Save results so far
        save('results2.mat',"loom_cvR","rec_cvR","pps_cvR","eps_cvR",...
            "loom_cvErr","rec_cvErr","pps_cvErr","eps_cvErr",...
            "loom_origD","rec_origD","pps_origA","eps_origA",...
            "loom_predD","rec_predD","pps_predA","eps_predA",...
            "loom_rmax","rec_rmax","pps_rmax","eps_rmax",...
            "loom_testR","rec_testR","pps_testR","eps_testR",...
            "loom_trfWgt","rec_trfWgt","pps_trfWgt","eps_trfWgt",...
            "allModels","saveSub","lambda","lambdavals");
    end
end

% % Average the data across subjects
% % Cross-validation accuracy
% mean_cvR_loom = mean(loom_cvR(:,:,dDim));
% mean_cvR_rec = mean(rec_cvR(:,:,dDim));
% mean_cvR_pps = mean(pps_cvR(:,:,aDim));
% mean_cvR_eps = mean(eps_cvR(:,:,aDim));
% 
% std_cvR_loom
% std_cvR_rec
% std_cvR_pps
% std_cvR_eps
% 
% % Cross-validation error
% mean_cvErr_loom = mean(loom_cvErr(:,:,dDim));
% mean_cvErr_rec = mean(rec_cvErr(:,:,dDim));
% mean_cvErr_pps = mean(pps_cvErr(:,:,aDim));
% mean_cvErr_eps = mean(eps_cvErr(:,:,aDim));
% 
% % Original spatial location (Distance or Azi)
% mean_origD_loom = mean(loom_origD,2);
% mean_origD_rec = mean(rec_origD,2);
% mean_origA_pps = mean(pps_origA,2);
% mean_origA_eps = mean(eps_origA,2);
% 
% % Model performance: largest r on the training set
% mean_rmax_loom = mean(loom_rmax);
% mean_rmax_rec = mean(rec_rmax);
% mean_rmax_pps = mean(pps_rmax);
% mean_rmax_eps = mean(eps_rmax);
% 
% % Model performance: r on the test set
% mean_testR_loom = mean(loom_testR);
% mean_testR_rec = mean(rec_testR);
% mean_testR_pps = mean(pps_testR);
% mean_testR_eps = mean(eps_testR);
% 
% % TRF weights
% mean_trfWgt_loom = mean(loom_trfWgt,2);
% mean_trfWgt_rec = mean(rec_trfWgt,2);
% mean_trfWgt_pps = mean(pps_trfWgt,2);
% mean_trfWgt_eps = mean(eps_trfWgt,2);
% 
% save('TRFresults.mat',"mean_cvR_loom","mean_cvR_rec","mean_cvR_pps","mean_cvR_eps",...
%     "mean_cvErr_loom","mean_cvErr_rec","mean_cvErr_pps","mean_cvErr_eps",...
%     "mean_origD_loom","mean_origD_rec","mean_origA_pps","mean_origA_eps",...
%     "mean_rmax_loom","mean_rmax_rec","mean_rmax_pps","mean_rmax_eps",...
%     "mean_testR_loom","mean_testR_rec","mean_testR_pps","mean_testR_eps",...
%     "mean_trfWgt_loom","mean_trfWgt_rec","mean_trfWgt_pps","mean_trfWgt_eps", ...
%     "lambdavals","lambda");

toPlot = 0;
if toPlot
% % Plotting the data % %
plotTRFres(dDim,aDim,fs,nfold)

% % Looming
% figure
% subplot(2,2,1), errorbar(1:numel(lambda),mean(loom.cv.r(:,:,dDim)),std(loom.cv.r(:,:,dDim))/sqrt(nfold-1),'linewidth',2)
% set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
% title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')
% 
% % Plot CV error
% subplot(2,2,2), errorbar(1:numel(lambda),mean(loom.cv.err(:,:,dDim)),std(loom.cv.err(:,:,dDim))/sqrt(nfold-1),'linewidth',2)
% set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
% title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')
% 
% % Plot reconstruction
% subplot(2,2,3), plot((1:length(loom.stest))/fs,loom.stest(:,dDim),'linewidth',2), hold on
% plot((1:length(loom.pred(:,dDim)))/fs,mean(loom.pred(:,dDim),2),'linewidth',2), hold off, xlim([0,3]), axis square, grid on
% title('Reconstruction'), xlabel('Time (s)'), ylabel('Amplitude (a.u.)'), legend('Orig','Pred')
% 
% % Plot test accuracy
% subplot(2,2,4), bar(1,loom.rmax), hold on, bar(2,mean(loom.test.r(dDim))), hold off
% set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
% title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
% 
% % EPS
% % Plot CV accuracy
% figure
% subplot(2,2,1), errorbar(1:numel(lambda),mean(eps.cv.r(:,:,aDim)),std(eps.cv.r(:,:,aDim))/sqrt(nfold-1),'linewidth',2)
% set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
% title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')
% 
% % Plot CV error
% subplot(2,2,2), errorbar(1:numel(lambda),mean(eps.cv.err(:,:,aDim)),std(eps.cv.err(:,:,aDim))/sqrt(nfold-1),'linewidth',2)
% set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
% title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')
% 
% % Plot reconstruction
% subplot(2,2,3), plot((1:length(eps.stest))/fs,eps.stest(:,aDim),'linewidth',2), hold on
% plot((1:length(eps.pred(:,aDim)))/fs,mean(eps.pred(:,aDim),2),'linewidth',2), hold off, xlim([0,3]), axis square, grid on
% title('Reconstruction'), xlabel('Time (s)'), ylabel('Amplitude (a.u.)'), legend('Orig','Pred')
% 
% % Plot test accuracy
% subplot(2,2,4), bar(1,eps.rmax), hold on, bar(2,mean(eps.test.r(aDim))), hold off
% set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
% title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
% 
% % Rec
% figure
% subplot(2,2,1), errorbar(1:numel(lambda),mean(rec.cv.r(:,:,dDim)),std(rec.cv.r(:,:,dDim))/sqrt(nfold-1),'linewidth',2)
% set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
% title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')
% 
% % Plot CV error
% subplot(2,2,2), errorbar(1:numel(lambda),mean(rec.cv.err(:,:,dDim)),std(rec.cv.err(:,:,dDim))/sqrt(nfold-1),'linewidth',2)
% set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
% title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')
% 
% % Plot reconstruction
% subplot(2,2,3), plot((1:length(rec.stest))/fs,rec.stest(:,dDim),'linewidth',2), hold on
% plot((1:length(rec.pred(:,dDim)))/fs,mean(rec.pred(:,dDim),2),'linewidth',2), hold off, xlim([0,3]), axis square, grid on
% title('Reconstruction'), xlabel('Time (s)'), ylabel('Amplitude (a.u.)'), legend('Orig','Pred')
% 
% % Plot test accuracy
% subplot(2,2,4), bar(1,rec.rmax), hold on, bar(2,mean(rec.test.r(dDim))), hold off
% set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
% title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
% 
% % PPS
% figure
% subplot(2,2,1), errorbar(1:numel(lambda),mean(pps.cv.r(:,:,aDim)),std(pps.cv.r(:,:,aDim))/sqrt(nfold-1),'linewidth',2)
% set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
% title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')
% 
% % Plot CV error
% subplot(2,2,2), errorbar(1:numel(lambda),mean(pps.cv.err(:,:,aDim)),std(pps.cv.err(:,:,aDim))/sqrt(nfold-1),'linewidth',2)
% set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
% title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')
% 
% % Plot reconstruction
% subplot(2,2,3), plot((1:length(pps.stest))/fs,pps.stest(:,aDim),'linewidth',2), hold on
% plot((1:length(pps.pred(:,aDim)))/fs,mean(pps.pred(:,aDim),2),'linewidth',2), hold off, xlim([0,3]), axis square, grid on
% title('Reconstruction'), xlabel('Time (s)'), ylabel('Amplitude (a.u.)'), legend('Orig','Pred')
% 
% % Plot test accuracy
% subplot(2,2,4), bar(1,pps.rmax), hold on, bar(2,mean(pps.test.r(aDim))), hold off
% set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
% title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
% 
% 
% figure;
% subplot(1,2,1)
% mTRFplot(loom.model,'trf');
% title('TRF weights')
% grid on
% 
% % Plot GFP
% subplot(1,2,2)
% mTRFplot(loom.model,'gfp',[],'all');
% title('TRF Global Field Power')
% grid on

end