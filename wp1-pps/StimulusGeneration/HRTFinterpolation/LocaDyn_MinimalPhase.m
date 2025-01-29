function hM2=LocaDyn_MinimalPhase(hM,opt)

% AA_MinimalPhase               - Replace the phase by minimal phase
%
% hM2=AA_MinimalPhase(hM);
%
% Replace the phase spectrum of hM by the minimal phase. The minimal phase
% is calculated using hilbert transformation.
%
% Input:
%  hM: data matrix with impulse respnoses (IR): 
%      dim 1: time in samples
%      dim 2: each IR
%      dim 3: each record channel
%  opt: IR (0, default), Amplitudespectrum (1)
%
% Output:
%  hM2: hM with minimal phase
%
% Piotr Majdak, 26.5.2006

if ~exist('opt','var')
    opt=0;
end

n=size(hM,1);
itnr=size(hM,2);
rec=size(hM,3);
hM2=zeros(size(hM));

for jj=1:rec
    for ii=1:itnr
        if ~opt
            h=squeeze(hM(:,ii,jj));
            % decompose signal
            amp1=abs(fft(h));
        else
            amp1=squeeze(hM(:,ii,jj));
        end

        % transform
        amp2=amp1;
        an2u=-imag(hilbert(log(amp1))); % minimal phase

        % reconstruct signal from amp2 and an2u
        % build a symmetrical phase 
        an2u=an2u(1:floor(n/2)+1);
        an2u=[an2u; -flipud(an2u(2:end+mod(n,2)-1))];
        an2=an2u-round(an2u/2/pi)*2*pi;  % wrap around +/-pi: wrap(x)=x-round(x/2/pi)*2*pi
        % amplitude
        amp2=amp2(1:floor(n/2)+1);
        amp2=[amp2; flipud(amp2(2:end+mod(n,2)-1))];
        % back to time domain
        h2=real(ifft(amp2.*exp(1i*an2)));
        hM2(:,ii,jj)=h2;
    end
end