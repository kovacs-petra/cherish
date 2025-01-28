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
% Interpolate HRTF horizontally
[Xhorint, TOAcheckhor] = LocaDyn_InterpolateHRTFs_horPlane(Xhor,1,2,1);

% Interpolate HRTFs radially
azi = [90 0 -90];
ele = 0;
radiusresolution = .01;
[Xradint, TOAcheck] = LocaDyn_InterpolateHRTFs_distance(Xhor,radiusresolution,azi,ele,2,0);

%% Combine interpolated HRTFs
Xint = X;
Xint.Data.IR= [Xhorint.Data.IR; Xradint.Data.IR];
Xint.SourcePosition=[Xhorint.SourcePosition; Xradint.SourcePosition];
Xint = SOFAupdateDimensions(Xint);

%% Plot original and interpolated HRTFs
% Xall = [Xhor,Xint];
% for xx = 1: length(Xall)
% for R = 1:2
%   fs = Xall(xx).Data.SamplingRate;
%   pos = Xall(xx).SourcePosition;
%   idx = pos(:,1) == azi & pos(:,2) == ele;
%   hM=double(squeeze(Xall(xx).Data.IR(idx,R,:)));
% 
%   [r,i]=sort(pos(idx,3));
%   hM = hM(i,:);
% 
%   M=(20*log10(abs(fft(hM')')));
%   M=M(:,1:floor(size(M,2)/2));  % only positive frequencies
% 
%   Mt=(20*log10(abs(hM)));
% 
%   figure('Name',num2str(R))
%   subplot(1,2,1) 
%   time = 0:1/fs*1000:(size(Mt,2)-1)/fs*1000;
%   h=surface(time,r,Mt(:,:));
%   shading flat
%   xlabel('Time (ms)');
%   ylabel('Radius (m)');
% 
%   subplot(1,2,2) 
%   freq = 0:fs/size(hM,2):(size(M,2)-1)*fs/size(hM,2);
%   h=surface(freq,r,M(:,:));
%   shading flat
%   xlabel('Frequency (Hz)');
%   ylabel('Radius (m)');
% end
% end

%% Listening example
sig = sig_triwave(400,fs,3);
% traj.r = linspace(0.2,0.2,length(sig));
traj.r = [linspace(0.2,0.2,length(traj.r)/3),...
    linspace(0.2,2,length(traj.r)/3), linspace(2,2,length(traj.r)/3)];
traj.azi = -90*ones(length(traj.r),1);
% traj.azi = [linspace(1,1,length(traj.r)/3),...
%     linspace(1,89,length(traj.r)/3), linspace(89,89,length(traj.r)/3)];
traj.ele = ele*ones(length(traj.r),1);
[out, aziActual, eleActual, rActual, idx] = SOFAspat((sig./traj.r)',Xint,traj.azi,traj.ele,traj.r);
soundsc(out,fs)

%% Save
SOFAsave("SCUT_KEMAR_radius_all_interp2.sofa", Xint);