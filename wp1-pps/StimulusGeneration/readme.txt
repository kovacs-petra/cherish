Folders named after the date of creation:
	- Spatialized stimuli: .mat files
	- Nonspatialzed stimuli (input to BRT): .wav files
	- Stimuli parameters: .csv file
	- ILD trajectory for each stimulus: .png file
	% Each folder contains 10 unique stimuli in one trajectory (loom, rec, PPS, or EPS) and one side (left or right)
	% This will result in 20 stimuli per side and condition (i.e., 20 unique stimuli per condition). Each will be played 5 times in the experiment.

Folder bugreport: report about ILD problem in the BRT Renderer App. The solution was to disable reverb and use the default BRT HRTF set.

Folder HRTFinterpolation: attempts to interpolate the SCUT distance-dependent HRTF set and thus remove some artefacts. Set not used in the final version.

Files:
- generateStimuli.m: main stimulus generation script
	- oscsend.m: sends osc messages from MATLAB
	- BRTspat.m: spatializes the cues using the BRT Renderer App
- getStimuliArray.m: organizes all stimuli into one big array including all parameters and audio, and containing each stimulus 5 times
- checkBRT.m: contains snippets to visualize any artefacts, e.g the ILDs
	- extractILD.m: extracts ILD info from a signal
- rampAndFilter.m: uses SOFAspat for spatialization rather than BRT; not used in the final version
- sig_triwave.m: generates triangle wave examples 


