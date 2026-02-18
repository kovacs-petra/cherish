% Single lag decoder analysis
% Set number of observations
NSub = 22;
allSub = 3:24;
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

pps         = struct;
pps.names   = {'distance','azimuth'};

eps         = struct;
eps.names   = {'distance','azimuth'};

% Init saving the data
loom_cvR    = []; rec_cvR    = []; pps_cvR    = []; eps_cvR    = [];
loom_cvErr  = []; rec_cvErr  = []; pps_cvErr  = []; eps_cvErr  = [];
loom_cvT    = []; rec_cvT    = []; pps_cvT    = []; eps_cvT    = [];

for i = 3:round(size(dir_file,1)) 
    filename = dir_file(i).name;
    if contains(filename,'.set')
        sub = str2double(filename(7:8));

    if ismember(sub,allSub)

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
        tmin = 100;
        tmax = 500;
        % Lambda chosen based on the previous decoding analyses (R values)
        lambda_loom = 10.^4;
        lambda_rec = 10.^4;
        lambda_pps = 10.^2;
        lambda_eps = 10.^0;
        [loom.cv,loom.cv_t] = mTRFcrossval(loom.strain,loom.rtrain,fs,Dir,tmin,tmax,lambda_loom,'type','single','zeropad',0);
        [rec.cv,rec.cv_t] = mTRFcrossval(rec.strain,rec.rtrain,fs,Dir,tmin,tmax,lambda_rec,'type','single','zeropad',0);
        [pps.cv,pps.cv_t] = mTRFcrossval(pps.strain,pps.rtrain,fs,Dir,tmin,tmax,lambda_pps,'type','single','zeropad',0);
        [eps.cv,eps.cv_t] = mTRFcrossval(eps.strain,eps.rtrain,fs,Dir,tmin,tmax,lambda_eps,'type','single','zeropad',0);

        % Save data
        % Cross-validation accuracy (Pearson's r)
        loom_cvR    = [loom_cvR; loom.cv.r(:,:,dDim,:)];
        rec_cvR     = [rec_cvR; rec.cv.r(:,:,dDim,:)]; 
        pps_cvR     = [pps_cvR; pps.cv.r(:,:,aDim,:)]; 
        eps_cvR     = [eps_cvR; eps.cv.r(:,:,aDim,:)];

        % Cross-validation error
        loom_cvErr  = [loom_cvErr; loom.cv.err(:,:,dDim,:)]; 
        rec_cvErr   = [rec_cvErr; rec.cv.err(:,:,dDim,:)]; 
        pps_cvErr   = [pps_cvErr; pps.cv.err(:,:,aDim,:)]; 
        eps_cvErr   = [eps_cvErr; eps.cv.err(:,:,aDim,:)];

        % T
        loom_cvT    = [loom_cvT; loom.cv_t];
        rec_cvT    = [rec_cvT; rec.cv_t];
        pps_cvT    = [pps_cvT; pps.cv_t];
        eps_cvT    = [eps_cvT; eps.cv_t];
        
        % Save subject numbers
        saveSub = [saveSub,sub];
        chanlocs = loom.chanlocs;

        % Save results so far
        save('results_singleLag.mat',"loom_cvR","rec_cvR","pps_cvR","eps_cvR",...
            "loom_cvErr","rec_cvErr","pps_cvErr","eps_cvErr",...
            "loom_cvT","rec_cvT","pps_cvT","eps_cvT",...
            "saveSub","chanlocs");
    end
    end
end

% Compute mean and variance
macc_loom = squeeze(mean(loom_cvR))'; vacc_loom = squeeze(var(loom_cvR))';
merr_loom = squeeze(mean(loom_cvErr))'; verr_loom = squeeze(var(loom_cvErr))';

macc_rec = squeeze(mean(rec_cvR))'; vacc_rec = squeeze(var(rec_cvR))';
merr_rec = squeeze(mean(rec_cvErr))'; verr_rec = squeeze(var(rec_cvErr))';

macc_pps = squeeze(mean(pps_cvR))'; vacc_pps = squeeze(var(pps_cvR))';
merr_pps = squeeze(mean(pps_cvErr))'; verr_pps = squeeze(var(pps_cvErr))';

macc_eps = squeeze(mean(eps_cvR))'; vacc_eps = squeeze(var(eps_cvR))';
merr_eps = squeeze(mean(eps_cvErr))'; verr_eps = squeeze(var(eps_cvErr))';

% Compute variance bound
xacc_loom = [-fliplr(loom_cvT),-loom_cvT]; yacc_loom = [fliplr(macc_loom-sqrt(vacc_loom/nfold)),macc_loom+sqrt(vacc_loom/nfold)];
xerr_loom = [-fliplr(loom_cvT),-loom_cvT]; yerr_loom = [fliplr(merr_loom-sqrt(verr_loom/nfold)),merr_loom+sqrt(verr_loom/nfold)];

xacc_rec = [-fliplr(rec_cvT),-rec_cvT]; yacc_rec = [fliplr(macc_rec-sqrt(vacc_rec/nfold)),macc_rec+sqrt(vacc_rec/nfold)];
xerr_rec = [-fliplr(rec_cvT),-rec_cvT]; yerr_rec = [fliplr(merr_rec-sqrt(verr_rec/nfold)),merr_rec+sqrt(verr_rec/nfold)];

xacc_pps = [-fliplr(pps_cvT),-pps_cvT]; yacc_pps = [fliplr(macc_pps-sqrt(vacc_pps/nfold)),macc_pps+sqrt(vacc_pps/nfold)];
xerr_pps = [-fliplr(pps_cvT),-pps_cvT]; yerr_pps = [fliplr(merr_pps-sqrt(verr_pps/nfold)),merr_pps+sqrt(verr_pps/nfold)];

xacc_eps = [-fliplr(eps_cvT),-eps_cvT]; yacc_eps = [fliplr(macc_eps-sqrt(vacc_eps/nfold)),macc_eps+sqrt(vacc_eps/nfold)];
xerr_eps = [-fliplr(eps_cvT),-eps_cvT]; yerr_eps = [fliplr(merr_eps-sqrt(verr_eps/nfold)),merr_eps+sqrt(verr_eps/nfold)];

% Plotting the conditions separately
%% Loom
% Plot accuracy
figure(Name='Looming decoding accuracy')
subplot(1,2,1), h = fill(xacc_loom,yacc_loom,'b','edgecolor','none'); hold on
set(h,'facealpha',0.01), xlim([tmin,tmax]), axis square, grid on
plot(-fliplr(loom_cvT),fliplr(macc_loom),'linewidth',2), hold off
title('Reconstruction Accuracy'), xlabel('Time lag (ms)'), ylabel('Correlation')

% Plot error
subplot(1,2,2)
h = fill(xerr_loom,yerr_loom,'b','edgecolor','none'); hold on
set(h,'facealpha',0.01), xlim([tmin,tmax]), axis square, grid on
plot(-fliplr(loom_cvT),fliplr(merr_loom),'linewidth',2), hold off
title('Reconstruction Error'), xlabel('Time lag (ms)'), ylabel('MSE')
sgtitle('Looming decoding accuracy')
%% Rec
% Plot accuracy
figure(Name='Receding decoding accuracy')
subplot(1,2,1), h = fill(xacc_rec,yacc_rec,'b','edgecolor','none'); hold on
set(h,'facealpha',0.01), xlim([tmin,tmax]), axis square, grid on
plot(-fliplr(rec_cvT),fliplr(macc_rec),'linewidth',2), hold off
title('Reconstruction Accuracy'), xlabel('Time lag (ms)'), ylabel('Correlation')

% Plot error
subplot(1,2,2)
h = fill(xerr_rec,yerr_rec,'b','edgecolor','none'); hold on
set(h,'facealpha',0.01), xlim([tmin,tmax]), axis square, grid on
plot(-fliplr(rec_cvT),fliplr(merr_rec),'linewidth',2), hold off
title('Reconstruction Error'), xlabel('Time lag (ms)'), ylabel('MSE')
sgtitle('Receding decoding accuracy')

%% PPS
% Plot accuracy
figure(Name='PPS decoding accuracy')
subplot(1,2,1), h = fill(xacc_pps,yacc_pps,'b','edgecolor','none'); hold on
set(h,'facealpha',0.01), xlim([tmin,tmax]), axis square, grid on
plot(-fliplr(pps_cvT),fliplr(macc_pps),'linewidth',2), hold off
title('Reconstruction Accuracy'), xlabel('Time lag (ms)'), ylabel('Correlation')

% Plot error
subplot(1,2,2)
h = fill(xerr_pps,yerr_pps,'b','edgecolor','none'); hold on
set(h,'facealpha',0.01), xlim([tmin,tmax]), axis square, grid on
plot(-fliplr(pps_cvT),fliplr(merr_pps),'linewidth',2), hold off
title('Reconstruction Error'), xlabel('Time lag (ms)'), ylabel('MSE')
sgtitle('PPS decoding accuracy')

%% EPS
% Plot accuracy
figure(Name='EPS decoding accuracy')
subplot(1,2,1), h = fill(xacc_eps,yacc_eps,'b','edgecolor','none'); hold on
set(h,'facealpha',0.01), xlim([tmin,tmax]), axis square, grid on
plot(-fliplr(eps_cvT),fliplr(macc_eps),'linewidth',2), hold off
title('Reconstruction Accuracy'), xlabel('Time lag (ms)'), ylabel('Correlation')

% Plot error
subplot(1,2,2)
h = fill(xerr_eps,yerr_eps,'b','edgecolor','none'); hold on
set(h,'facealpha',0.01), xlim([tmin,tmax]), axis square, grid on
plot(-fliplr(eps_cvT),fliplr(merr_eps),'linewidth',2), hold off
title('Reconstruction Error'), xlabel('Time lag (ms)'), ylabel('MSE')
sgtitle('EPS decoding accuracy')

%% Plotting the conditions together
% Accuracy
figure(Name='Decoding accuracy')
subplot(1,2,1)
    l = fill(xacc_loom,yacc_loom,'r','edgecolor','none'); hold on
    r = fill(xacc_rec,yacc_rec,'b','edgecolor','none'); 
    p = fill(xacc_pps,yacc_pps,'m','edgecolor','none'); 
    e = fill(xacc_eps,yacc_eps,'g','edgecolor','none'); 
set([l,r,p,e],'facealpha',0.01), xlim([tmin,tmax]), axis square, grid on
    ll = plot(-fliplr(loom_cvT),fliplr(macc_loom),'r','linewidth',2);
    rr = plot(-fliplr(rec_cvT),fliplr(macc_rec),'b','linewidth',2);
    pp = plot(-fliplr(pps_cvT),fliplr(macc_pps),'m','linewidth',2);
    ee = plot(-fliplr(eps_cvT),fliplr(macc_eps),'g','linewidth',2); hold off
title('Reconstruction accuracy'), xlabel('Time lag (ms)'), ylabel('Correlation')
sgtitle('Decoding accuracy')
% legend('Looming','Receding','Rotating in PPS','Rotating in EPS');

% Error
subplot(1,2,2)
    l = fill(xerr_loom,yerr_loom,'r','edgecolor','none'); hold on
    r = fill(xerr_rec,yerr_rec,'b','edgecolor','none'); 
    p = fill(xerr_pps,yerr_pps,'m','edgecolor','none'); 
    e = fill(xerr_eps,yerr_eps,'g','edgecolor','none'); 
set([l,r,p,e],'facealpha',0.01), xlim([tmin,tmax]), axis square, grid on
    plot(-fliplr(loom_cvT),fliplr(merr_loom),'r','linewidth',2);
    plot(-fliplr(rec_cvT),fliplr(merr_rec),'b','linewidth',2);
    plot(-fliplr(pps_cvT),fliplr(merr_pps),'m','linewidth',2);
    plot(-fliplr(eps_cvT),fliplr(merr_eps),'g','linewidth',2), hold off
title('Reconstruction error'), xlabel('Time lag (ms)'), ylabel('MSE')
sgtitle('Decoding error')
