% Load HRTF set to obtain sampling rate

SOFAdbURL('http://www.sofacoustics.org/data')

database = 'scut';       
HRTFfilename = 'SCUT_KEMAR_radius_all.sofa';
fullfn = fullfile(SOFAdbPath, 'database', database, HRTFfilename);
X = SOFAload(fullfn);
SOFAinfo(X)
%% Prepare HRTF

X=SOFAremoveVariable(X,'MeasurementSourceAudioChannel');
X=SOFAremoveVariable(X,'MeasurementAudioLatency');
if isfield(X,'GLOBAL__NCProperties'), X=rmfield(X,'GLOBAL__NCProperties'); end

% Remove elevation
Xhor=X;
idx = X.SourcePosition(:,2)==0;
Xhor.Data.IR=Xhor.Data.IR(idx,:,:);
Xhor.SourcePosition=Xhor.SourcePosition(idx,:);
Xhor = SOFAupdateDimensions(Xhor);

%% Interpolate HRTF
ele = 0;

% Interpolate HRTF horizontally
azires = 1;
r = min(X.SourcePosition(:,3));
[XhorintPPS, TOAcheckhor] = LocaDyn_InterpolateHRTFs_horPlane(Xhor,azires,ele,r,2,1);
r = max(X.SourcePosition(:,3));
% r = 2; % EPS 2 m
[XhorintEPS, TOAcheckhor] = LocaDyn_InterpolateHRTFs_horPlane(Xhor,azires,ele,r,2,1);

% Interpolate HRTFs radially
radires = .01;
azi = 90;
[XradintLeft, TOAcheck] = LocaDyn_InterpolateHRTFs_distance(Xhor,radires,azi,ele,2,0);
azi = 270;
[XradintRight, TOAcheck] = LocaDyn_InterpolateHRTFs_distance(Xhor,radires,azi,ele,2,0);

%% Save the four HRTFs separately
SOFAsave("HRTF_PPS.sofa", XhorintPPS);
SOFAsave("HRTF_EPS.sofa", XhorintEPS);
SOFAsave("HRTF_left.sofa", XradintLeft);
SOFAsave("HRTF_right.sofa", XradintRight);

% %% Combine interpolated HRTFs -> SOFAspat seems to have problems with the combined HRTF set -> keep them separate!
% Xint = X;
% Xint.Data.IR= [XhorintPPS.Data.IR;XhorintEPS.Data.IR; XradintLeft.Data.IR; XradintRight.Data.IR];
% Xint.SourcePosition=[XhorintPPS.SourcePosition;XhorintEPS.SourcePosition; XradintLeft.SourcePosition; XradintRight.SourcePosition];
% Xint = SOFAupdateDimensions(Xint);

%% Plot original and interpolated HRTFs
Xall = [Xhor,Xint];
for xx = 1: length(Xall)
for R = 1:2
  fs = Xall(xx).Data.SamplingRate;
  pos = Xall(xx).SourcePosition;
  idx = pos(:,1) == azi & pos(:,2) == ele;
  hM=double(squeeze(Xall(xx).Data.IR(idx,R,:)));

  [r,i]=sort(pos(idx,3));
  hM = hM(i,:);

  M=(20*log10(abs(fft(hM')')));
  M=M(:,1:floor(size(M,2)/2));  % only positive frequencies

  Mt=(20*log10(abs(hM)));

  figure('Name',num2str(R))
  subplot(1,2,1) 
  time = 0:1/fs*1000:(size(Mt,2)-1)/fs*1000;
  h=surface(time,r,Mt(:,:));
  shading flat
  xlabel('Time (ms)');
  ylabel('Radius (m)');

  subplot(1,2,2) 
  freq = 0:fs/size(hM,2):(size(M,2)-1)*fs/size(hM,2);
  h=surface(freq,r,M(:,:));
  shading flat
  xlabel('Frequency (Hz)');
  ylabel('Radius (m)');
end
end

%%
[x,y,z]= sph2cart(deg2rad(Xint.SourcePosition(:,1)),deg2rad(Xint.SourcePosition(:,2)),Xint.SourcePosition(:,3)); 
figure; 
scatter3(x,y,z)

%% Listening examples
sig = sig_triwave(400,44100,3);

%% Looming sound
traj.r = linspace(1,0.2,length(sig));
traj.azi = 90*ones(length(traj.r),1);
traj.ele = ele*ones(length(traj.r),1);
[out, aziActual, eleActual, rActual, idx] = SOFAspat((sig./traj.r)',XradintLeft,traj.azi,traj.ele,traj.r);
figure; plot(out)
soundsc(out,fs)

%% Clockwise rotating sound, PPS
traj.r = linspace(0.2,0.2,length(sig));
traj.azi = linspace(-90,90,length(sig));
traj.ele = ele*ones(length(traj.r),1);
[out, aziActual, eleActual, rActual, idx] = SOFAspat((sig./traj.r)',XhorintPPS,traj.azi,traj.ele,traj.r);
figure; plot(out)
soundsc(out,fs)

%% Clockwise rotating sound, EPS
traj.r = linspace(2,2,length(sig));
traj.azi = linspace(-90,90,length(sig));
traj.ele = ele*ones(length(traj.r),1);
[out, aziActual, eleActual, rActual, idx] = SOFAspat((sig./traj.r)',XhorintEPS,traj.azi,traj.ele,traj.r);
figure; plot(out)
soundsc(out,fs)