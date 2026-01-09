function [pred_loom,test_loom,pred_rec,test_rec,...
    pred_pps,test_pps,pred_eps,test_eps] = runTRF(stim_loom,stim_rec,stim_pps,stim_eps,...
         resp_loom,resp_rec,resp_pps,resp_eps)

% addpath '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\tools\mTRF-Toolbox\mtrf';

NSub = 22;

% Prepare stimulus and eeg structures if not already inputted
if ~nargin
    [stim_loom,stim_rec,stim_pps,stim_eps,...
        resp_loom,resp_rec,resp_pps,resp_eps] = prepTRF;
end

% Partition data into training and test sets
nfold = 10;
testTrial = 1;

% Init output for partitioning
% strain_loom = zeros(1,NSub);
% rtrain_loom = zeros(1,NSub);
% strain_rec = zeros(1,NSub);
% rtrain_rec = zeros(1,NSub);
% strain_pps = zeros(1,NSub);
% rtrain_pps = zeros(1,NSub);
% strain_eps = zeros(1,NSub);
% rtrain_eps = zeros(1,NSub);
% 
% strain_loom = zeros(1,NSub);
% rtrain_loom = zeros(1,NSub);
% strain_rec = zeros(1,NSub);
% rtrain_rec = zeros(1,NSub);
% strain_pps = zeros(1,NSub);
% rtrain_pps = zeros(1,NSub);
% strain_eps = zeros(1,NSub);
% rtrain_eps = zeros(1,NSub);

for i = 1:NSub
[strain_loom{i},rtrain_loom{i},stest_loom{i},rtest_loom{i}] = ...
    mTRFpartition(stim_loom.data(i,:),resp_loom.data(i,:),nfold,testTrial);
[strain_rec{i},rtrain_rec{i},stest_rec{i},rtest_rec{i}] = ...
    mTRFpartition(stim_rec.data(i,:),resp_rec.data(i,:),nfold,testTrial);
[strain_pps{i},rtrain_pps{i},stest_pps{i},rtest_pps{i}] = ...
    mTRFpartition(stim_pps.data(i,:),resp_pps.data(i,:),nfold,testTrial);
[strain_eps{i},rtrain_eps{i},stest_eps{i},rtest_eps{i}] = ...
    mTRFpartition(stim_eps.data(i,:),resp_eps.data(i,:),nfold,testTrial);


% Optimize the encoder's ability to predict EEG features from new stimulus 
% data: tune the regularization parameter using an efficient leave-one-out 
% cross-validation (CV) procedure
fs = stim_loom.fs;
Dir = 1;
tmin = 0;
tmax = 500;
lambda = 10.^(-6:2:6);
cv_loom = mTRFcrossval(strain_loom(i),rtrain_loom(i),fs,Dir,tmin,tmax,lambda,'zeropad',0,'fast',1,'dim',2);
cv_rec = mTRFcrossval(strain_rec(i),rtrain_rec(i),fs,Dir,tmin,tmax,lambda,'zeropad',0,'fast',1,'dim',2);
cv_pps = mTRFcrossval(strain_pps(i),rtrain_pps(i),fs,Dir,tmin,tmax,lambda,'zeropad',0,'fast',1,'dim',2);
cv_eps = mTRFcrossval(strain_eps(i),rtrain_eps(i),fs,Dir,tmin,tmax,lambda,'zeropad',0,'fast',1,'dim',2);

% Get the optimal regularization value from the cv structure
[~,idx_loom(i)] = max(mean(cv_loom.r));
[~,idx_rec(i)] = max(mean(cv_rec.r));
[~,idx_pps(i)] = max(mean(cv_pps.r));
[~,idx_eps(i)] = max(mean(cv_eps.r));

% Train the model using this lambda value
model_loom(i) = mTRFtrain(strain_loom(i),rtrain_loom(i),fs,Dir,tmin,tmax,lambda(idx_loom(i)),'zeropad',0);
model_rec(i) = mTRFtrain(strain_rec(i),rtrain_rec(i),fs,Dir,tmin,tmax,lambda(idx_rec(i)),'zeropad',0);
model_pps(i) = mTRFtrain(strain_pps(i),rtrain_pps(i),fs,Dir,tmin,tmax,lambda(idx_pps(i)),'zeropad',0);
model_eps(i) = mTRFtrain(strain_eps(i),rtrain_eps(i),fs,Dir,tmin,tmax,lambda(idx_eps(i)),'zeropad',0);

% Test the model
[pred_loom(i),test_loom(i)] = mTRFpredict(stest_loom(i),rtest_loom(i),model_loom(i),'zeropad',0);
[pred_rec(i),test_rec(i)] = mTRFpredict(stest_rec(i),rtest_rec(i),model_rec(i),'zeropad',0);
[pred_pps(i),test_pps(i)] = mTRFpredict(stest_pps(i),rtest_pps(i),model_pps(i),'zeropad',0);
[pred_eps(i),test_eps(i)] = mTRFpredict(stest_eps(i),rtest_eps(i),model_eps(i),'zeropad',0);


end

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