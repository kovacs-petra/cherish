function [g,TRidx,par,Pi] = VBAPgain(P, X, TR)
% [g,TRidx,par,Pi] = VBAPgain(P, X, TR)
%   P: requested position (cartesian)
%   X: coordinates of the available grid
%   TR: faces of the grid (can be obtained by TR=freeBoundary(delaunayTriangulation(LS)); 
%   
%   g: vector with gains
%   TRidx: index of the face found for the calculations
%   par: vector with the plane parameters 
%   Pi: intersection of (O,P) with the plane (P1,P2,P3)
%
%   #Author: Piotr Majdak (2017): Initial implementation for spatial reverb in the LAS
%   #Author: Piotr Majdak (2021): Error robustness improved for LocaDyn


for ii=1:size(TR,1)
    P1=X(TR(ii,1),:);
    if prod(prod(P==P1)), g=[1 0 0]; TRidx=ii; par=g; Pi=P1; return; end
    P2=X(TR(ii,2),:);
    if prod(prod(P==P2)), g=[0 1 0]; TRidx=ii; par=g; Pi=P2; return; end
    P3=X(TR(ii,3),:);
    if prod(prod(P==P3)), g=[0 0 1]; TRidx=ii; par=g; Pi=P3; return; end
    if cond([P; P2-P1; P3-P1])>1e10,
%         warning('eeee');
    else
      a=[P; P2-P1; P3-P1]'\(P-P1)'; % https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
      if a(2)>=-2*eps && a(2)<=1+eps && a(3)>=-2*eps && a(3)<=1+eps && (a(2)+a(3))<=1+eps && a(1)<=1+eps, 
        TRidx=ii;
        par=a';
        Pi=P-par(1)*P; % Pi: intersection of (O,P) with the plane (P1,P2,P3)
        g=[1-par(3)-par(2) par(2) par(3)]; % loudspeaker gains
        return;
      end
    end
end
% warning on
warning('Face surrounding the requested direction not found!');
g=NaN; TRidx=NaN; par=NaN; Pi=NaN; 
