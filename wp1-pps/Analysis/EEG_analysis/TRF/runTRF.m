addpath '\\kfs.oeaw.ac.at\fileserver\Projektdaten\CherISH\tools\mTRF-Toolbox\mtrf';
model2000=mTRFtrain(stim.data,resp.data,stim.fs,1,-100,2000,0.1);

figure;
subplot(1,2,1)
mTRFplot(model1000,'trf');
title('TRF weights')
grid on

% Plot GFP
subplot(1,2,2)
mTRFplot(model1000,'gfp',[],'all');
title('TRF Global Field Power')
grid on

figure; plot(model1000.w(:,:,60)');