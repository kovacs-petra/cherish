function [Xrec, TOAcheck] = LocaDyn_InterpolateHRTFs_horPlane(X,azires,TOAremoval,fig)
% function Xrec = LocaDyn_InterpolateHRTFs_horPlane(X,TOAremoval,fig)
%
%  interpolating HRTFs to a denser grid on the horizontal plane.
%  TOA estimation % Based on TOA Model (Ziegelwanger and Majdak, 2014)
%
%  input:	X: reference SOFA object
%         azires: angular resolution along azimuth
%         TOAremoval: 0 .. no removal; 1: .. circshift; 2: .. minimum phase
%         fig: plot figures (true/false)
%
%  output:	Xrec: modified SOFA object
%  
%
%  #Author: Piotr Majdak, original function name LocaDyn_InterpolateHRTFs (01.2022)
%  #Author: Robert Baumgartner: only horizontal plane for Dynamtes (11.2022)

 if ~exist('fig','var')
     % third parameter does not exist, so default it to something
      fig=false;
 end
fs=X.Data.SamplingRate; 

%% Plot Original HRTFs
if fig
  figT = figure;%('Position',[ 73       386.14       1164.6       333.71]);
  subplot(2,3,1); 
  SOFAplotHRTF(X,'etchorizontal',1);
  title('Original HRIRs, left');
  subplot(2,3,2); 
  SOFAplotHRTF(X,'etchorizontal',2);
  title('Original HRIRs, right');
  subplot(2,3,3); 
  SOFAplotHRTF(X,'itdhorizontal');
  title('Original HRIRs');

  figF = figure;%('Position',[ 73       386.14       1164.6       333.71]);
  subplot(2,2,1); 
  SOFAplotHRTF(X,'maghorizontal',1); %ylim([-90 270]);
  title('Original HRTFs, left');
  subplot(2,2,2); 
  SOFAplotHRTF(X,'maghorizontal',2); %ylim([-90 270]);
  title('Original HRTFs, right');
end


%% Fit the TOA model

% [~,res]=LocaDyn_ziegelwanger2014_adapt(X,4,1,1);
[~,res]=ziegelwanger2014(X,4,1,1); % use thr method on lp-filtered IRs, AMT 1.1.1
p=res.p_offaxis(:,1);
txt = 'head radius=%0.3g cm; head center=[%5.3g, %5.3g, %5.3g] cm; constant delay=%5.3g ms; ear azimuth=%5.3g°; ear elevation=%5.3g°\n';
leftEar = sprintf([' Left:  ' txt], ...
   p(1)*100, p(2)*100, p(3)*100, p(4)*100, p(5)*1000, rad2deg(p(6)), rad2deg(p(7)));
p=res.p_offaxis(:,2);
rightEar = sprintf([' Right: ' txt], ...
   p(1)*100, p(2)*100, p(3)*100, p(4)*100, p(5)*1000, rad2deg(p(6)), rad2deg(p(7)));

TOAcheck = ['TOA model sanity check:' 10 leftEar rightEar];
disp(TOAcheck);

X = SOFAaddVariable(X,'TOAModel','NR',[res.p_offaxis; zeros(X.API.N-size(res.p_offaxis,1),X.API.R)]);

%% remove the TOA: Xnd
switch TOAremoval, 
  case 0 % No removal
    Xnd=X;
  case 1 % Circular shift
    Xnd=X;
    toa=zeros(X.API.M, X.API.R); 
    for r=1:X.API.R
      toa(:,r)=ziegelwanger2014_offaxis(res.p_offaxis(:,r),[deg2rad(X.SourcePosition(:,1)) deg2rad(X.SourcePosition(:,2)) deg2rad(X.SourcePosition(:,3))]);
      for ii=1:X.API.M
        Xnd.Data.IR(ii,r,:)=circshift(squeeze(X.Data.IR(ii,r,:)),round((-toa(ii,r)+res.p_offaxis(5,r))*fs)); % res.p_offaxis(5,r) == tau
      end
    end
  case 2 % Minimum-phase
    Xnd=X;
    Xnd.Data.IR=shiftdim(LocaDyn_MinimalPhase(shiftdim(Xnd.Data.IR,2)),1);
    toa=zeros(X.API.M, X.API.R); 
    for r=1:X.API.R
      toa(:,r)=ziegelwanger2014_offaxis(res.p_offaxis(:,r),[deg2rad(X.SourcePosition(:,1)) deg2rad(X.SourcePosition(:,2)) deg2rad(X.SourcePosition(:,3))]);
      for ii=1:X.API.M
        Xnd.Data.IR(ii,r,:)=circshift(squeeze(Xnd.Data.IR(ii,r,:)),round((res.p_offaxis(5,r))*fs)); % res.p_offaxis(5,r) == tau
      end
    end
  otherwise
    disp('Selected TOAremoval does not exist');
    return
end
Xnd=SOFAupdateDimensions(Xnd);

% if fig
%   figure('Position',[ 73       386.14       1164.6       333.71]);
%   subplot(1,2,1); 
%   SOFAplotHRTF(Xnd,'etchorizontal',1);
%   title('No-Delay HRIRs, left');
%   subplot(1,2,2); 
%   SOFAplotHRTF(Xnd,'etchorizontal',2);
%   title('No-Delay HRIRs, right');
% end

%% Add a silent pole at the bottom
Xsp=Xnd; 
Xsp.Data.IR(end+1,:,:)=zeros(1,X.API.R,X.API.N);
Xsp.SourcePosition(end+1,:)=[0,-90,Xsp.SourcePosition(1,3)];
Xsp.API.M=X.API.M+1;
Xsp=SOFAupdateDimensions(Xsp);

%% Define to be interpolated positions -> horizontal plane
azi=Xsp.SourcePosition(:,1);
ele=Xsp.SourcePosition(:,2);
r=Xsp.SourcePosition(:,3);
[x,y,z]=sph2cart(deg2rad(azi),deg2rad(ele),r);
rmean=mean(r);
P=[x y z];
T = freeBoundary(delaunayTriangulation(P));
% if display
%   figure;
%   trisurf(T, P(:,1), P(:,2), P(:,3), 'FaceColor', 'w','EdgeColor','b', 'FaceAlpha',0.5);
%   title('Subdivision: original'); 
%   axis equal
%   view(75,10);
% end
% new.azi = 0:azires:360-azires;
% new.ele = zeros(size(new.azi));
% new.r = r(1).*ones(size(new.azi));
distances = [0.2,2];
new.azi = repmat(0:azires:360-azires,1,length(distances));
new.ele = zeros(size(new.azi));
new.r = repmat(distances,1,length(new.azi)/length(distances));
[new.x,new.y,new.z]=sph2cart(deg2rad(new.azi),deg2rad(new.ele),new.r);
Pnew=[new.x;new.y;new.z]';
% [a,e,r]=cart2sph(Pnew(:,1),Pnew(:,2),Pnew(:,3)); to double check the transformation


%% Define target grid
Xnew=Xsp;
Xnew.Data.IR=zeros(length(Pnew),Xsp.API.R,Xsp.API.N); 
Xnew.SourcePosition=[new.azi(:), new.ele(:), new.r(:)];
MInspect=NaN; % stop at MInspect measurement for inspection
tic
for ii=1:length(Pnew)  
  [g,Tidx] = VBAPgain(Pnew(ii,:), P, T);
  if isnan(g), error('face for VBAPing not found'); end
  for r=1:Xsp.API.R
    A=squeeze(Xsp.Data.IR(T(Tidx,:),r,:))'; % A: three IRs from the face around the new position
    Xnew.Data.IR(ii,r,:)=A*g'; % matrix-vector product
  end
  if mod(ii,500)==0, disp([num2str(ii) '/' num2str(length(Pnew))]); end
  if ii==MInspect
    S1=Xsp.SourcePosition(T(Tidx,1),:);
    S2=Xsp.SourcePosition(T(Tidx,2),:);
    S3=Xsp.SourcePosition(T(Tidx,3),:);
    disp(S1); disp(S2); disp(S3); 
    disp(Xnew.SourcePosition(ii,:)); 
    disp('---'); 
    figure(1);
    plot(A); hold on;
    plot(squeeze(Xnew.Data.IR(ii,r,:)),'R');
    legend({['Gain: ' num2str(g(1))], ...
            ['Gain: ' num2str(g(2))], ...
            ['Gain: ' num2str(g(3))], ...
            'New'});
    title(['#' num2str(ii) ' of ' num2str(length(Pnew))]);
    waitforbuttonpress;
    close(1);
  end
end
toc
Xnew.Data.Delay=zeros(1,Xsp.API.R);
Xnew=SOFAupdateDimensions(Xnew);

% if fig
%   figure('Position',[ 73       386.14       1164.6       333.71]);
%   subplot(1,2,1); 
%   SOFAplotHRTF(Xnew,'etchorizontal',1);
%   title('Interporalted HRIRs, left');
%   subplot(1,2,2); 
%   SOFAplotHRTF(Xnew,'etchorizontal',2);
%   title('Interporalted HRIRs, right');
% end

%% reconstruct the TOA in the HRIRs

Xrec=Xnew;
if TOAremoval,
  toa=zeros(Xrec.API.M, Xrec.API.R); 
  for r=1:X.API.R
    toa(:,r)=ziegelwanger2014_offaxis(Xrec.TOAModel(:,r), [deg2rad(Xrec.SourcePosition(:,1)) deg2rad(Xrec.SourcePosition(:,2)) Xrec.SourcePosition(:,3)]);
    for ii=1:Xrec.API.M
      Xrec.Data.IR(ii,r,:)=circshift(squeeze(Xrec.Data.IR(ii,r,:)),round(-(-toa(ii,r)+Xrec.TOAModel(5,r))*fs));
    end
  end
end
Xrec=SOFAremoveVariable(Xrec,'TOAModel');
Xrec=SOFAupdateDimensions(Xrec); 

if fig
  figure(figT);%'Position',[ 73       386.14       1164.6       333.71]);
  subplot(2,3,4); 
  SOFAplotHRTF(Xrec,'etchorizontal',1);
  title('Interpolated HRIRs, left');
  subplot(2,3,5); 
  SOFAplotHRTF(Xrec,'etchorizontal',2);
  title('Interpolated HRIRs, right');
  subplot(2,3,6); 
  SOFAplotHRTF(Xrec,'itdhorizontal');
  title('Interpolated HRIRs');

  figure(figF);%'Position',[ 73       386.14       1164.6       333.71]);
  subplot(2,2,3); 
  SOFAplotHRTF(Xrec,'maghorizontal',1); %ylim([-90 270]);
  title('Interpolated HRTFs, left');
  subplot(2,2,4); 
  SOFAplotHRTF(Xrec,'maghorizontal',2); %ylim([-90 270]);
  title('Interpolated HRTFs, right');
end
