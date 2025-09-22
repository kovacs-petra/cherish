Externalization pilot for CherISH wp1

- Stimulus generation script
- Stimulus presentation script

stimStruct.mat:
- fs: sampling freq
- stim: 4-D structure -- f0 (310-620 Hz) x distance (0.2, 2, 0.4, 0.9, 1.3, 1.8 m) x azimuth (90 deg, -90 deg) x stimulus type (triangle wave, white noise)

stimStructEq.mat: the same as stimStruct.mat, but loudness equalized using the script loudnessEq.m