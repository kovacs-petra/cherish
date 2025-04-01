% Plot some stimulus info to check how precise the stimulus
% generation turned out with BRT

%% Check delay between desired and actual radius (for radial motion)
time = linspace(0,dur,length(rC));
figure;
plot(time,rC);
[~,~,radius]=cart2sph(EmitterPosition(1,1,:), EmitterPosition(1,2,:),EmitterPosition(1,3,:));
rActual = squeeze(radius);
hold on
plot(M,rActual)

%% Plot ILDs for stimuli
% (Loaded one by one from folder named after date of creation)
% (Figure saved into ild.png)
[ild310,time310] = extractILD(out',fs,0.1);
[ild340,time340] = extractILD(out',fs,0.1);
[ild380,time380] = extractILD(out',fs,0.1);
[ild410,time410] = extractILD(out',fs,0.1);
[ild450,time450] = extractILD(out',fs,0.1);
[ild480,time480] = extractILD(out',fs,0.1);
[ild520,time520] = extractILD(out',fs,0.1);
[ild550,time550] = extractILD(out',fs,0.1);
[ild590,time590] = extractILD(out',fs,0.1);
[ild620,time620] = extractILD(out',fs,0.1);

figure; sgtitle('Rotating in EPS, default sadie 48 kHz HRTF, BRIR disabled');
subplot(4,3,1); plot(time310',ild310); ylim([-14,20]); title('310 Hz'); ylabel('ILD (dB)'); xlabel('Time (s)');
hold on; yline(0,'k--');
subplot(4,3,2); plot(time340',ild340); ylim([-14,20]); title('340 Hz'); ylabel('ILD (dB)'); xlabel('Time (s)');
hold on; yline(0,'k--');
subplot(4,3,3); plot(time380',ild380); ylim([-14,20]); title('380 Hz'); ylabel('ILD (dB)'); xlabel('Time (s)');
hold on; yline(0,'k--');
subplot(4,3,4); plot(time410',ild410); ylim([-14,20]); title('410 Hz'); ylabel('ILD (dB)'); xlabel('Time (s)');
hold on; yline(0,'k--');
subplot(4,3,5); plot(time450',ild450); ylim([-14,20]); title('450 Hz'); ylabel('ILD (dB)'); xlabel('Time (s)');
hold on; yline(0,'k--');
subplot(4,3,6); plot(time480',ild480); ylim([-14,20]); title('480 Hz'); ylabel('ILD (dB)'); xlabel('Time (s)');
hold on; yline(0,'k--');
subplot(4,3,7); plot(time520',ild520); ylim([-14,20]); title('520 Hz'); ylabel('ILD (dB)'); xlabel('Time (s)');
hold on; yline(0,'k--');
subplot(4,3,8); plot(time550',ild550); ylim([-14,20]); title('550 Hz'); ylabel('ILD (dB)'); xlabel('Time (s)');
hold on; yline(0,'k--');
subplot(4,3,9); plot(time590',ild590); ylim([-14,20]); title('590 Hz'); ylabel('ILD (dB)'); xlabel('Time (s)');
hold on; yline(0,'k--');
subplot(4,3,11); plot(time620',ild620); ylim([-14,20]); title('620 Hz'); ylabel('ILD (dB)'); xlabel('Time (s)');
hold on; yline(0,'k--');


