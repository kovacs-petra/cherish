function plotTRFres(dDim,aDim,fs,nfold,lambda,lambdavals,NSub)

%% Aggregate data files
% Subject 1, subject 2-11, and subject 12-24 run and saved separately (for
% efficiency) - concatenate
load("results_s03-10_20260116.mat");

    loom_cvR_s1    = loom_cvR;
    rec_cvR_s1     = rec_cvR; 
    pps_cvR_s1     = pps_cvR; 
    eps_cvR_s1     = eps_cvR;

    % Cross-validation error
    loom_cvErr_s1  = loom_cvErr; 
    rec_cvErr_s1   = rec_cvErr; 
    pps_cvErr_s1   = pps_cvErr; 
    eps_cvErr_s1   = eps_cvErr;

    % Original spatial location (Distance or Azi)
    loom_origD_s1  = loom_origD; 
    rec_origD_s1   = rec_origD; 
    pps_origA_s1   = pps_origA; 
    eps_origA_s1   = eps_origA;

    % Predicted spatial location (Distance or Azi)
    loom_predD_s1  = loom_predD; 
    rec_predD_s1   = rec_predD; 
    pps_predA_s1   = pps_predA; 
    eps_predA_s1   = eps_predA;

    % Model performance: largest r on the training set
    loom_rmax_s1   = loom_rmax; 
    rec_rmax_s1    = rec_rmax; 
    pps_rmax_s1    = pps_rmax; 
    eps_rmax_s1    = eps_rmax;

    % Model performance: r on the test set
    loom_testR_s1  = loom_testR; 
    rec_testR_s1   = rec_testR; 
    pps_testR_s1   = pps_testR; 
    eps_testR_s1   = eps_testR;

    % TRF weights: what spatial location the brain represents
    % loom_trfWgt_s1 = loom_trfWgt; 
    % rec_trfWgt_s1  = rec_trfWgt; 
    % pps_trfWgt_s1  = pps_trfWgt; 
    % eps_trfWgt_s1  = eps_trfWgt;
    allModels_s1 = allModels;

load("results_s11-13_20260116.mat"); 

    loom_cvR_s2    = loom_cvR;
    rec_cvR_s2     = rec_cvR; 
    pps_cvR_s2     = pps_cvR; 
    eps_cvR_s2     = eps_cvR;

    % Cross-validation error
    loom_cvErr_s2  = loom_cvErr; 
    rec_cvErr_s2   = rec_cvErr; 
    pps_cvErr_s2   = pps_cvErr; 
    eps_cvErr_s2   = eps_cvErr;

    % Original spatial location (Distance or Azi)
    loom_origD_s2  = loom_origD; 
    rec_origD_s2   = rec_origD; 
    pps_origA_s2   = pps_origA; 
    eps_origA_s2   = eps_origA;

    % Predicted spatial location (Distance or Azi)
    loom_predD_s2  = loom_predD; 
    rec_predD_s2   = rec_predD; 
    pps_predA_s2   = pps_predA; 
    eps_predA_s2   = eps_predA;

    % Model performance: largest r on the training set
    loom_rmax_s2   = loom_rmax; 
    rec_rmax_s2    = rec_rmax; 
    pps_rmax_s2    = pps_rmax; 
    eps_rmax_s2    = eps_rmax;

    % Model performance: r on the test set
    loom_testR_s2  = loom_testR; 
    rec_testR_s2   = rec_testR; 
    pps_testR_s2   = pps_testR; 
    eps_testR_s2   = eps_testR;

    % TRF weights: what spatial location the brain represents
    % loom_trfWgt_s2 = loom_trfWgt; 
    % rec_trfWgt_s2  = rec_trfWgt; 
    % pps_trfWgt_s2  = pps_trfWgt; 
    % eps_trfWgt_s2  = eps_trfWgt;
    allModels_s2 = allModels;

load("results_s14-24_20260116.mat"); % s12-24

    loom_cvR_s3    = loom_cvR;
    rec_cvR_s3     = rec_cvR; 
    pps_cvR_s3     = pps_cvR; 
    eps_cvR_s3     = eps_cvR;

    % Cross-validation error
    loom_cvErr_s3  = loom_cvErr; 
    rec_cvErr_s3   = rec_cvErr; 
    pps_cvErr_s3   = pps_cvErr; 
    eps_cvErr_s3   = eps_cvErr;

    % Original spatial location (Distance or Azi)
    loom_origD_s3  = loom_origD; 
    rec_origD_s3   = rec_origD; 
    pps_origA_s3   = pps_origA; 
    eps_origA_s3   = eps_origA;

    % Predicted spatial location (Distance or Azi)
    loom_predD_s3  = loom_predD; 
    rec_predD_s3   = rec_predD; 
    pps_predA_s3   = pps_predA; 
    eps_predA_s3   = eps_predA;

    % Model performance: largest r on the training set
    loom_rmax_s3   = loom_rmax; 
    rec_rmax_s3    = rec_rmax; 
    pps_rmax_s3    = pps_rmax; 
    eps_rmax_s3    = eps_rmax;

    % Model performance: r on the test set
    loom_testR_s3  = loom_testR; 
    rec_testR_s3   = rec_testR; 
    pps_testR_s3   = pps_testR; 
    eps_testR_s3   = eps_testR;

    % TRF weights: what spatial location the brain represents
    % loom_trfWgt_s3 = loom_trfWgt; 
    % rec_trfWgt_s3  = rec_trfWgt; 
    % pps_trfWgt_s3  = pps_trfWgt; 
    % eps_trfWgt_s3  = eps_trfWgt;
    allModels_s3 = allModels;

%% Master variables
    loom_cvR = [loom_cvR_s1;loom_cvR_s2;loom_cvR_s3];
    rec_cvR  = [rec_cvR_s1;rec_cvR_s2;rec_cvR_s3];
    pps_cvR  = [pps_cvR_s1;pps_cvR_s2;pps_cvR_s3]; 
    eps_cvR  = [eps_cvR_s1;eps_cvR_s2;eps_cvR_s3];

    % Cross-validation error
    loom_cvErr  = [loom_cvErr_s1;loom_cvErr_s2;loom_cvErr_s3];
    rec_cvErr   = [rec_cvErr_s1;rec_cvErr_s2;rec_cvErr_s3];
    pps_cvErr   = [pps_cvErr_s1;pps_cvErr_s2;pps_cvErr_s3]; 
    eps_cvErr   = [eps_cvErr_s1;eps_cvErr_s2;eps_cvErr_s3];

    % Original spatial location (Distance or Azi)
    loom_origD  = [loom_origD_s1,loom_origD_s2,loom_origD_s3]; 
    rec_origD   = [rec_origD_s1,rec_origD_s2,rec_origD_s3]; 
    pps_origA   = [pps_origA_s1,pps_origA_s2,pps_origA_s3]; 
    eps_origA   = [eps_origA_s1,eps_origA_s2,eps_origA_s3];

    % Predicted spatial location (Distance or Azi)
    loom_predD  = [loom_predD_s1,loom_predD_s2,loom_predD_s3]; 
    rec_predD   = [rec_predD_s1,rec_predD_s2,rec_predD_s3]; 
    pps_predA   = [pps_predA_s1,pps_predA_s2,pps_predA_s3]; 
    eps_predA   = [eps_predA_s1,eps_predA_s2,eps_predA_s3];

    % Model performance: largest r on the training set
    loom_rmax   = [loom_rmax_s1,loom_rmax_s2,loom_rmax_s3]; 
    rec_rmax    = [rec_rmax_s1,rec_rmax_s2,rec_rmax_s3]; 
    pps_rmax    = [pps_rmax_s1,pps_rmax_s2,pps_rmax_s3]; 
    eps_rmax    = [eps_rmax_s1,eps_rmax_s2,eps_rmax_s3]; 

    % Model performance: r on the test set
    loom_testR  = [loom_testR_s1,loom_testR_s2,loom_testR_s3]; 
    rec_testR   = [rec_testR_s1,rec_testR_s2,rec_testR_s3]; 
    pps_testR   = [pps_testR_s1,pps_testR_s2,pps_testR_s3];
    eps_testR   = [eps_testR_s1,eps_testR_s2,eps_testR_s3];

    % TRF models
    allModels.loom.w = {allModels_s1.loom.w{1:8,1},allModels_s2.loom.w{9:11,1},...
        allModels_s3.loom.w{12:end,1}}';
    allModels.loom.b = {allModels_s1.loom.b{1:8,1},allModels_s2.loom.b{9:11,1},...
        allModels_s3.loom.b{12:end,1}}';

    allModels.rec.w = {allModels_s1.rec.w{1:8,1},allModels_s2.rec.w{9:11,1},...
        allModels_s3.rec.w{12:end,1}}';
    allModels.rec.b = {allModels_s1.rec.b{1:8,1},allModels_s2.rec.b{9:11,1},...
        allModels_s3.rec.b{12:end,1}}';

    allModels.pps.w = {allModels_s1.pps.w{1:8,1},allModels_s2.pps.w{9:11,1},...
        allModels_s3.pps.w{12:end,1}}';
    allModels.pps.b = {allModels_s1.pps.b{1:8,1},allModels_s2.pps.b{9:11,1},...
        allModels_s3.pps.b{12:end,1}}';

    allModels.eps.w = {allModels_s1.eps.w{1:8,1},allModels_s2.eps.w{9:11,1},...
        allModels_s3.eps.w{12:end,1}}';
    allModels.eps.b = {allModels_s1.eps.b{1:8,1},allModels_s2.eps.b{9:11,1},...
        allModels_s3.eps.b{12:end,1}}';

    % TRF weights
    % Transform model structs to have subject as the 4th dimension
    for ss = 1:NSub
        model_loom.wi(:,:,:,ss) = allModels.loom.w{ss,1};
        model_loom.bi(:,:,ss) = allModels.loom.b{ss,1};

        model_rec.wi(:,:,:,ss) = allModels.rec.w{ss,1};
        model_rec.bi(:,:,ss) = allModels.rec.b{ss,1};

        model_pps.wi(:,:,:,ss) = allModels.pps.w{ss,1};
        model_pps.bi(:,:,ss) = allModels.pps.b{ss,1};

        model_eps.wi(:,:,:,ss) = allModels.eps.w{ss,1};
        model_eps.bi(:,:,ss) = allModels.eps.b{ss,1};
    end

    % Take the avg of w and b by subject (4th and 3rd dim, respectively)
    model_loom.w = mean(model_loom.wi,4);
    model_loom.b = mean(model_loom.bi,3);
    model_loom.t = allModels.loom.t{12,1}; % it's the same for everyone
    model_loom.fs = 100;
    model_loom.Dir = -1;
    model_loom.type = 'multi';

    model_rec.w = mean(model_rec.wi,4);
    model_rec.b = mean(model_rec.bi,3);
    model_rec.t = allModels.rec.t{12,1}; % it's the same for everyone
    model_rec.fs = 100;
    model_rec.Dir = -1;
    model_rec.type = 'multi';

    model_pps.w = mean(model_pps.wi,4);
    model_pps.b = mean(model_pps.bi,3);
    model_pps.t = allModels.pps.t{12,1}; % it's the same for everyone
    model_pps.fs = 100;
    model_pps.Dir = -1;
    model_pps.type = 'multi';

    model_eps.w = mean(model_eps.wi,4);
    model_eps.b = mean(model_eps.bi,3);
    model_eps.t = allModels.eps.t{12,1}; % it's the same for everyone
    model_eps.fs = 100;
    model_eps.Dir = -1;
    model_eps.type = 'multi';

%% Grand avgs
% Cross-validation accuracy
mean_cvR_loom = mean(loom_cvR(:,:,dDim));
mean_cvR_rec = mean(rec_cvR(:,:,dDim));
mean_cvR_pps = mean(pps_cvR(:,:,aDim));
mean_cvR_eps = mean(eps_cvR(:,:,aDim));

std_cvR_loom = std(loom_cvR(:,:,dDim));
std_cvR_rec = std(rec_cvR(:,:,dDim));
std_cvR_pps = std(pps_cvR(:,:,aDim));
std_cvR_eps = std(eps_cvR(:,:,aDim));

% Cross-validation error
mean_cvErr_loom = mean(loom_cvErr(:,:,dDim));
mean_cvErr_rec = mean(rec_cvErr(:,:,dDim));
mean_cvErr_pps = mean(pps_cvErr(:,:,aDim));
mean_cvErr_eps = mean(eps_cvErr(:,:,aDim));

std_cvErr_loom = std(loom_cvErr(:,:,dDim));
std_cvErr_rec = std(rec_cvErr(:,:,dDim));
std_cvErr_pps = std(pps_cvErr(:,:,aDim));
std_cvErr_eps = std(eps_cvErr(:,:,aDim));

% Original spatial location (Distance or Azi)
mean_origD_loom = mean(loom_origD,2);
mean_origD_rec = mean(rec_origD,2);
mean_origA_pps = mean(pps_origA,2);
mean_origA_eps = mean(eps_origA,2);

% Predicted spatial location (Distance or Azi)
mean_predD_loom = mean(loom_predD,2);
mean_predD_rec = mean(rec_predD,2);
mean_predA_pps = mean(pps_predA,2);
mean_predA_eps = mean(eps_predA,2);

% Model performance: largest r on the training set
mean_rmax_loom = mean(loom_rmax);
mean_rmax_rec = mean(rec_rmax);
mean_rmax_pps = mean(pps_rmax);
mean_rmax_eps = mean(eps_rmax);

% Model performance: r on the test set
mean_testR_loom = mean(loom_testR);
mean_testR_rec = mean(rec_testR);
mean_testR_pps = mean(pps_testR);
mean_testR_eps = mean(eps_testR);

% Save avgs
save("TRFresults_20260119.mat","mean_cvR_loom","mean_cvR_rec","mean_cvR_pps","mean_cvR_eps",...
    "std_cvR_loom","std_cvR_rec","std_cvR_pps","std_cvR_eps",...
    "mean_cvErr_loom","mean_cvErr_rec","mean_cvErr_pps","mean_cvErr_eps",...
    "std_cvErr_loom","std_cvErr_rec","std_cvErr_pps","std_cvErr_eps",...
    "mean_origD_loom","mean_origD_rec","mean_origA_pps","mean_origA_eps",...
    "mean_predD_loom","mean_predD_rec","mean_predA_pps","mean_predA_eps",...
    "mean_rmax_loom","mean_rmax_rec","mean_rmax_pps","mean_rmax_eps",...
    "mean_testR_loom","mean_testR_rec","mean_testR_pps","mean_testR_eps",...
    "model_loom","model_rec","model_pps","model_eps",...
    "lambda","lambdavals");

%% Looming
% Plot CV accuracy
figure(Name='Looming decoding results')
subplot(2,2,1), errorbar(1:numel(lambda),mean_cvR_loom,std_cvR_loom/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')

% Plot CV error
subplot(2,2,2), errorbar(1:numel(lambda),mean_cvErr_loom,std_cvErr_loom/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')

% Plot reconstruction
subplot(2,2,3), plot((1:length(mean_origD_loom))/fs,mean_origD_loom,'linewidth',2), hold on
plot((1:length(mean_predD_loom))/fs,mean_predD_loom,'linewidth',2), hold off, xlim([0,3]), axis square, grid on
title('Distance reconstruction'), xlabel('Time (s)'), ylabel('Distance (a.u.)'), legend('Orig','Pred')

% Plot test accuracy
subplot(2,2,4), bar(1,mean_rmax_loom), hold on, bar(2,mean_testR_loom), hold off
set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
sgtitle('Decoding results: LOOMING')

%% Receding
% Plot CV accuracy
figure(Name='Receding decoding results')
subplot(2,2,1), errorbar(1:numel(lambda),mean_cvR_rec,std_cvR_rec/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')

% Plot CV error
subplot(2,2,2), errorbar(1:numel(lambda),mean_cvErr_rec,std_cvErr_rec/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')

% Plot reconstruction
subplot(2,2,3), plot((1:length(mean_origD_rec))/fs,mean_origD_rec,'linewidth',2), hold on
plot((1:length(mean_predD_rec))/fs,mean_predD_rec,'linewidth',2), hold off, xlim([0,3]), axis square, grid on
title('Distance reconstruction'), xlabel('Time (s)'), ylabel('Distance (a.u.)'), legend('Orig','Pred')

% Plot test accuracy
subplot(2,2,4), bar(1,mean_rmax_rec), hold on, bar(2,mean_testR_rec), hold off
set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
sgtitle('Decoding results: RECEDING')

%% PPS
% Plot CV accuracy
figure(Name='PPS decoding results')
subplot(2,2,1), errorbar(1:numel(lambda),mean_cvR_pps,std_cvR_pps/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')

% Plot CV error
subplot(2,2,2), errorbar(1:numel(lambda),mean_cvErr_pps,std_cvErr_pps/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')

% Plot reconstruction
subplot(2,2,3), plot((1:length(mean_origA_pps))/fs,mean_origA_pps,'linewidth',2), hold on
plot((1:length(mean_predA_pps))/fs,mean_predA_pps,'linewidth',2), hold off, xlim([0,3]), axis square, grid on
title('Azimuth reconstruction'), xlabel('Time (s)'), ylabel('Azimuth (a.u.)'), legend('Orig','Pred')

% Plot test accuracy
subplot(2,2,4), bar(1,mean_rmax_pps), hold on, bar(2,mean_testR_pps), hold off
set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
sgtitle('Decoding results: ROTATING IN PPS')

%% EPS
% Plot CV accuracy
figure(Name='EPS decoding results')
subplot(2,2,1), errorbar(1:numel(lambda),mean_cvR_eps,std_cvR_eps/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')

% Plot CV error
subplot(2,2,2), errorbar(1:numel(lambda),mean_cvErr_eps,std_cvErr_eps/sqrt(nfold-1),'linewidth',2)
set(gca,'xtick',1:length(lambda),'xticklabel',lambdavals), xlim([0,numel(lambda)+1]), axis square, grid on
title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')

% Plot reconstruction
subplot(2,2,3), plot((1:length(mean_origA_eps))/fs,mean_origA_eps,'linewidth',2), hold on
plot((1:length(mean_predA_eps))/fs,mean_predA_eps,'linewidth',2), hold off, xlim([0,3]), axis square, grid on
title('Azimuth reconstruction'), xlabel('Time (s)'), ylabel('Azimuth (a.u.)'), legend('Orig','Pred')

% Plot test accuracy
subplot(2,2,4), bar(1,mean_rmax_eps), hold on, bar(2,mean_testR_eps), hold off
set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
sgtitle('Decoding results: ROTATING IN EPS')

%% TRF filter plots
figure(Name='TRF filters')
subplot(2,2,1)
mTRFplot(loom.model,'trf',dDim,"all"); % specify channels
grid on
title('Looming TRF filter')

subplot(2,2,2)
mTRFplot(rec.model,'trf',dDim,21); % specify channels
grid on
title('Receding TRF filter')

subplot(2,2,3)
mTRFplot(pps.model,'trf',aDim,21); % specify channels
grid on
title('Rotating in PPS TRF filter')

subplot(2,2,4)
mTRFplot(eps.model,'trf',aDim,21); % specify channels
grid on
title('Rotating in EPS TRF filter')
sgtitle('TRF filters, channels: Cz')

% subplot(2,2,1), plot((1:length(mean_trfWgt_loom))/fs,mean_trfWgt_loom,'linewidth',2), grid on
% title('Looming TRF filter'), xlabel('Time (s)'), ylabel('TRF weights'), xlim([0 0.61])
% 
% subplot(2,2,2), plot((1:length(mean_trfWgt_rec))/fs,mean_trfWgt_rec,'linewidth',2), grid on
% title('Receding TRF filter'), xlabel('Time (s)'), ylabel('TRF weights'), xlim([0 0.61])
% 
% subplot(2,2,3), plot((1:length(mean_trfWgt_pps))/fs,mean_trfWgt_pps,'linewidth',2), grid on
% title('Rotating in PPS TRF filter'), xlabel('Time (s)'), ylabel('TRF weights'), xlim([0 0.61])
% 
% subplot(2,2,4), plot((1:length(mean_trfWgt_eps))/fs,mean_trfWgt_eps,'linewidth',2), grid on
% title('Rotating in EPS TRF filter'), xlabel('Time (s)'), ylabel('TRF weights'), xlim([0 0.61])



