% Note: I'm averageing loom and rec conditions to try out, but in
% reality these should be values from the same condition but different
% subjects

% CV accuracy: Concatenate subjects' correlation values by lambda and
% calculate a mean
x = [loom.cv.r;rec.cv.r];
mean_agg = mean(x(:,:,dDim)); % Distance dimension for loom and rec, azi dimension for rest

% CV error: Do the same with the error values
y = [loom.cv.err;rec.cv.err];
mean_err = mean(y(:,:,dDim));

% Reconstruction of the trajectory: Time dimension needs to stay!
% Original values: Collect Distance values (loom, rec) in each col and average each col
% row by row. Do the same with Azi values for pps, eps.
z = [loom.stest(:,dDim),rec.stest(:,dDim)];
mean_d = mean(z,2);
% Predicted values: the same but in a different vector

% Model performance: Collect rmax values and avg them. Collect test.r
% values in the given dimension and avg them
rmaxes = [loom.rmax, rec.rmax];
mean_rmax = mean(rmaxes);

testRs = [loom.test.r(:,dDim),rec.test.r(:,dDim)];
mean_testR = mean(testRs);

% TRF filter: Average the model weights keeping the time (x) dimension, using
% the corresponding distance or azi (z) dimension. Channels (y) can be averaged, 
% as we don't have expectations about specific channels encoding the
% stimulus
avg_chan_model_l = mean(loom.model.w,2);
avg_chan_model_r = mean(rec.model.w,2);
weights = [avg_chan_model_l(:,dDim),avg_chan_model_r(:,dDim)];
mean_weights = mean(weights,2);

figure;
mTRFplot(loom.model,'trf',dDim,'all',[-500 100],0,1);

% T = table(x,'VariableNames',{'loomCVr'},'RowNames',{'s1','s2'});
% mean_agg = mean(T.loomCVr{:,:,dDim});