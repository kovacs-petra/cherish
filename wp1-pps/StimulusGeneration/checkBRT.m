% Plot some stimulus info to check how precise the stimulus
% generation turned out with BRT

% Check delay between desired and actual radius (for radial movt)
time = linspace(0,totalDur,length(rC));
figure;
plot(time,rC);
[~,~,radius]=cart2sph(EmitterPosition(1,1,:), EmitterPosition(1,2,:),EmitterPosition(1,3,:));
rActual = squeeze(radius);
hold on
plot(M,rActual)

% Check the stimulus
figure; plot(Data.Receiver')

% PPS vs. receding timing
figure; subplot(2,1,1); plot(linspace(0,max(cueParams.timeGoal),length(out)),out')
hold on; plot([0.4 0.4],[-1.5 1.5],'Color','k')
subplot(2,1,2); plot(linspace(0,max(cueParams.timeGoal),length(out)),out')
hold on; plot([0.4 0.4],[-1.5 1.5],'Color','k')

%% SPL levels in different stimuli
% Values from 10-Mar-2025-1445
load("spl_EPS.mat");
figure; sgtitle('1 m distance, -90 to +90 deg');
subplot(4,3,1); plot(spl_01'); ylim([-100,0]); title('310 Hz');
subplot(4,3,2); plot(spl_02'); ylim([-100,0]); title('340 Hz');
subplot(4,3,3); plot(spl_03'); ylim([-100,0]); title('380 Hz');
subplot(4,3,4); plot(spl_04'); ylim([-100,0]); title('410 Hz');
subplot(4,3,5); plot(spl_05'); ylim([-100,0]); title('450 Hz');
subplot(4,3,6); plot(spl_06'); ylim([-100,0]); title('480 Hz');
subplot(4,3,7); plot(spl_07'); ylim([-100,0]); title('520 Hz');
subplot(4,3,8); plot(spl_08'); ylim([-100,0]); title('550 Hz');
subplot(4,3,9); plot(spl_09'); ylim([-100,0]); title('590 Hz');
subplot(4,3,11); plot(spl_10'); ylim([-100,0]); title('620 Hz');

% save("spl_EPS.mat","spl_01","spl_02","spl_03","spl_04",...
%     "spl_05","spl_06","spl_07","spl_08","spl_09","spl_10");

load("spl_EPS_lr.mat")
figure; sgtitle('1 m distance, +90 to -90 deg');
subplot(4,3,1); plot(spl_01_rl'); ylim([-100,0]); title('310 Hz');
subplot(4,3,2); plot(spl_02_rl'); ylim([-100,0]); title('340 Hz');
subplot(4,3,3); plot(spl_03_rl'); ylim([-100,0]); title('380 Hz');
subplot(4,3,4); plot(spl_04_rl'); ylim([-100,0]); title('410 Hz');
subplot(4,3,5); plot(spl_05_rl'); ylim([-100,0]); title('450 Hz');
subplot(4,3,6); plot(spl_06_rl'); ylim([-100,0]); title('480 Hz');
subplot(4,3,7); plot(spl_07_rl'); ylim([-100,0]); title('520 Hz');
subplot(4,3,8); plot(spl_08_rl'); ylim([-100,0]); title('550 Hz');
subplot(4,3,9); plot(spl_09_rl'); ylim([-100,0]); title('590 Hz');
subplot(4,3,11); plot(spl_10_rl'); ylim([-100,0]); title('620 Hz');

% Plot ILD
figure; plot(spl_03(1,:)'-spl_03(2,:)')





