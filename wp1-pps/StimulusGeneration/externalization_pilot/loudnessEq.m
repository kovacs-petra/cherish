% Loudness equalize white noise and triangle wave stimuli for CherISH wp1
% externalization pilot round 3

% Initialize AMT
if ~exist("amt_start.m","file")
    amtPath = '\\KFS\Fileserver\ProjektDaten\CherISH\code\amtoolbox-full-1.4.0';
    addpath(amtPath);
    amt_start;
end

addpath('D:\TTK Projects\02 Figure-ground\Pilot_baba\stimulus\');
fs = 48e3;

[f0,d,a,t] = size(stimStruct.stim);
stimStructEq = struct;
stimStructEq.fs = fs;
stimStructEq.stim = cell(size(stimStruct.stim));

for ff = 1:f0
    for dd = 1:d
        for aa = 1:a
            r = stimStruct.stim{ff,dd,aa,2}{2}; % get distance

            white = stimStruct.stim{ff,dd,aa,2}{1};
            triangle = stimStruct.stim{ff,dd,aa,1}{1};

            % FFT
            N = length(white);
            f = linspace(0, fs/2, floor(N/2)+1);
            white_fft = abs(fft(white)).^2;
            tri_fft = abs(fft(triangle)).^2;

            % Take positive frequencies only
            white_psd = white_fft(1:length(f));
            tri_psd = tri_fft(1:length(f));

            % Get ISO 226 SPL contour at 60 phon (can be adjusted)
            [isoSPL, isoF] = iso226(60, f, true);  % in dB

            % Normalize to 0 dB at 1 kHz
            isoSPL = isoSPL - isoSPL(f == 1000);

            % Convert dB to linear scale
            isoWeight = 10.^(isoSPL/20);

            % Apply perceptual weighting
            white_weighted = white_psd .* isoWeight;
            tri_weighted = tri_psd .* isoWeight;

            % Estimate perceptual loudness as RMS of weighted spectrum
            white_loud = sqrt(mean(white_weighted));
            tri_loud = sqrt(mean(tri_weighted));

            % Match triangle wave to white noise loudness
            scaleFactor = white_loud / tri_loud;
            triangle_eq = triangle * scaleFactor;

            % Scale to 1
            scaleto = 0.2/r;
            scaledbW = 20 * log10(scaleto / max(white,[],"all"));
            scaledbT = 20 * log10(scaleto / max(triangle_eq,[],"all"));

            scaled1W = scaletodbspl(white(:,1),dbspl(white(:,1))+scaledbW);
            scaled2W = scaletodbspl(white(:,2),dbspl(white(:,2))+scaledbW);
            white_scaled = [scaled1W scaled2W];

            scaled1T = scaletodbspl(triangle_eq(:,1),dbspl(triangle_eq(:,1))+scaledbT);
            scaled2T = scaletodbspl(triangle_eq(:,2),dbspl(triangle_eq(:,2))+scaledbT);
            triangle_eq_scaled = [scaled1T scaled2T];

            stimStructEq.stim{ff,dd,aa,1}{1} = triangle_eq_scaled;
            stimStructEq.stim{ff,dd,aa,1}{2} = stimStruct.stim{ff,dd,aa,1}{2};
            stimStructEq.stim{ff,dd,aa,2}{1} = white_scaled;
            stimStructEq.stim{ff,dd,aa,2}{2} = stimStruct.stim{ff,dd,aa,2}{2};
        end
    end
end

save("stimStructEq.mat","stimStructEq");

