% Difference Wave Simulation Helper Function: Simulate One Subject's ERP Data File

% Purpose: This function simulates an ERP data file for one subject. The file
% contains trial-level ERP waveforms that have been generated for 100 trials
% consisting of 2 emotion conditions (condition A/Neutral and condition B/Happy), 
% 5 different 'actors', and 10 presentations of each emotion condition/actor combination.

% Notes: 
    % - The amplitude of trial-level waveforms are modulated by emotion, 
    %   trial presentation number, actor, subject, maternal sensitivity, and 
    %   age. Some modulations are condition-specific so that they will 
    %   be detectable in the trial-level *difference* wave (e.g., subject 
    %   intercept is added to the peak amplitude for Condition B only so
    %   that there are subject-level variations for the B-minus-A
    %   difference wave).
    % - The peak latency of trial-level waveforms are modulated by subject. 
    % - Note that peak amplitude values for each variable (e.g., emotion 
    %   condition) were chosen to produce a specific mean amplitude value 
    %   when extracted from C4 over a 300-500 ms time window. For example,
    %   condition A (Neutral) has a peak amplitude of 581 µV, which corresponds 
    %   to a mean amplitude of -10.002 µV over a 300-500 ms time window at C4. 
    % - Real noise from the infant noise profile are added in a subsequent
    %   step in simulateOneSample_forDiffWave. 

% ***See Section3_SimulatedData README.md available on the LME_MixedEffectsERPDifferenceWave
% GitHub for pipeline details: https://github.com/basclab/LME_MixedEffectsERPDifferenceWave/tree/main/Section3_SimulatedData

% Format:
    % ERP = simulateOneSubject_forDiffWave(sampleEmotionPeakAmp_condB_oneSubject, decayRate, sampleActorPeakAmpArray, leadField, sourceLocs)

% Inputs:
    % - sampleEmotionPeakAmp_condB_oneSubject: Peak amplitude value for this 
    %   subject's assigned maternal sensitivity and age group. 
    % - decayRate: Decay rate type used to specify whether a 'different' or 
    %   'same' decay rate is simulated across the two emotion conditions. 
    % - sampleActorPeakAmpArray: Array of actor intercepts randomly generated in
    %   the simulateOneSample function. The first value corresponds to the
    %   intercept for actor 01, the second value corresponds to actor 02, etc.
    % - leadField: Data structure created in
    %   DiffWaveSim_02_SimulateERPData.m for specifying the lead field, 
    %   electrode montage, electrodes of interest, and dipole orientation.
    % - sourceLocs: Index used to identify the dipole location from the
    %   leadField structure. This variable was created in
    %   DiffWaveSim_02_SimulateERPData.m.
    
% Other Requirements: 
    % - Needs MATLAB R2019a, EEGLAB v 2019_0, ERPLAB v 8.01, SEREEGA v 1.1.0
        % - For more information on EEGLAB, see: Delorme, A. & Makeig, S. (2004).
        %   EEGLAB: An open source toolbox for analysis of single-trial EEG dynamics
        %   including independent component analysis. https://sccn.ucsd.edu/eeglab/index.php
        % - For more information on ERPLAB, see: Lopez-Calderon, J., & Luck, S. J.
        %   (2014). ERPLAB: An open-source toolbox for the analysis of event-related
        %   potentials. https://erpinfo.org/erplab/    
        % - For more information on SEREEGA, see: Krol, L. R., Pawlitzki, J., Lotte, F.,
        %   Gramann, K., & Zander, T. O. (2018). SEREEGA: Simulating event-related EEG
        %   activity. https://github.com/lrkrol/SEREEGA
    % - Pediatric Head Atlas release v 1.1, Atlas 1 (0-2 years old) files
        % - Atlas files should be requested and downloaded from: https://www.pedeheadmod.net/pediatric-head-atlases/
        % - The folder containing the atlas files should then be added
        %   to the MATLAB path (via "Home" > "Set Path" > "Add with Subfolders").
    % - Filepath to the following file used during processing:
        % - binDescriptorFilename: File specifying each bin's number, label, 
        %   and 5-digit event markers. This file is created by the
        %   DiffWaveSim_01_CreateBinDescriptorFile.m script.
        
% Function Steps:
    % 1. Define simulation parameters
    % 2. Define function for generating ERP signal class
    % 3. Generate trial-level waveforms with parameters and function from steps 1-2
    % 4. Update event marker preceding codes with trial presentation number
    % 5. Extract bin-based epochs
    % 6. Calculate trial-level ERPs
    % 7. Delete extra simulated electrodes

% Outputs:
    % - ERP: ERP file containing the subject's trial-level waveforms.

% Usage Example:
    % >> sampleEmotionPeakAmp_condB_oneSubject = 581; 
    % >> decayRate = 'different';
    % >> sampleActorPeakAmpArray = [736.9235, 821.8267, -655.0656, -39.9697, -488.5581];
    % >> leadField = lf_generate_frompha('0to2','128','labels',{'E36','E104'});
    % >> sourceLocs = lf_get_source_nearest(leadField, [10 46 18]);
    % >> leadField.orientation(sourceLocs,:) = [0.57 -0.70 -0.01];
    % >> ERP = simulateOneSubject_forDiffWave(sampleEmotionPeakAmp_condB_oneSubject, decayRate, sampleActorPeakAmpArray, leadField, sourceLocs)

% Copyright 2024 Megan J. Heise, Serena K. Mon, Lindsay C. Bowman
% Brain and Social Cognition Lab, University of California Davis, Davis, CA, USA.

% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function ERP = simulateOneSubject_forDiffWave(sampleEmotionPeakAmp_condB_oneSubject, decayRate, sampleActorPeakAmpArray, leadField, sourceLocs)
%% 1. DEFINE SIMULATION PARAMETERS 

    % Specify filepath of bin descriptor file
    binDescriptorFilename = 'C:\Users\basclab\Desktop\Section3_SimulatedData\LMESimulation_BinDescriptorFile.txt';
        
    % Specify each stimulus' 3-digit preceding code used for creating event  
    % markers. Each preceding code corresponds to one emotion condition and actor. 
    % The first digit (e.g., 6) corresponds to the emotion condition (e.g., A)
    % and the last two digits correspond to the actor ID. For more
    % information, see the LMESimulation_EventMarkerMappingKey.xlsx spreadsheet
    emotionPrecCodes = ["601", "602", "603", "604", "605", "301", "302", "303", "304", "305"];
    % NOTE: The above array is ordered based on emotion condition such that all
    % markers starting with "6" belong to emotion condition A and are listed
    % first. Markers starting with "3" belong to condition B and are listed
    % second. This order is important for ensuring that the correct emotion 
    % condition population mean is assigned in step 3.

    % Define parameters for simulated epoch window 
    preStimPeriod = 200; % Pre-stimulus/baseline period (ms)
    postStimPeriod = 600; % Post-stimulus period (ms)
    samplingRate = 1000; % Sampling rate (Hz)

    % Define parameters for presentation of each emotion condition/actor
    presentN = 10; % Number of presentations of each stimulus (emotion condition/actor)
    presentNArray = repmat(pad(string(1:presentN),2,'left','0')',length(emotionPrecCodes),1); % Create a string array with repeating values of 1, 2, ..., presentN for each emotion condition/actor combination

    % Create epoch structure based on above parameters: each stimulus is
    % presented 10 times, data are sampled at 1000 Hz, and the epoch consists of
    % a 200 ms pre-stimulus period and 600 ms post-stimulus period.
    epochs = struct('n', presentN, 'srate', samplingRate, 'length', preStimPeriod+postStimPeriod,'prestim',preStimPeriod);

    % Define parameters for each emotion condition (A and B)   
    emotionPeakAmp_condA = 581;
    emotionPeakAmpArray = [emotionPeakAmp_condA, sampleEmotionPeakAmp_condB_oneSubject]; % Peak amplitude population mean for each emotion condition distribution (one distribution per condition; value for condition B will depend on the subject's assigned maternal sensitivity and age group) 
    emotionPeakAmpDv = 0; % Peak amplitude deviation (this value is 0 because we introduce waveform noise via real noise profiles in a subsequent step)
    
    % Peak amplitude decay rate/slope (i.e., change in amplitude with each successive presentation) 
    if strcmp(decayRate, 'different') % Specify condition A slope depending on decayRate
        emotionPeakAmpSlope_condA = -261; % Value corresponds to a mean amplitude increase of 0.499 µV (different slope from condition B)
    elseif strcmp(decayRate, 'same')
        emotionPeakAmpSlope_condA = -784; % Value corresponds to a mean amplitude increase of 1.499 µV (same slope as condition B)
    else
        error('decayRate must be either ''different'' or ''same''.')
    end
    emotionPeakAmpSlope_condB = -784; % Value corresponds to a mean amplitude increase of 1.499 µV for condition B

    emotionPeakLatency = preStimPeriod+400; % Peak latency (400 ms post-stimulus)
    emotionPeakWindow = 200; % Window length for simulated peak (i.e., a peak spanning 300-500 ms has a window length of 200 ms)

    % For the subject's peak latency shift, one value is drawn from the
    % following array using a uniform sampling distribution and used to
    % shift the subject's peak from the emotionPeakLatency specified above. 
    % This latency shift is constant across all trials and corresponds to a
    % minimal change in mean amplitude between subjects (ranging from 
    % 0.054 to 0 µV). 
    subjectPeakLatencyShiftArray = [-20 20]; 
    
    % Define parameters for subject-specific intercept added to peak amplitude
    % for condition B only
    subjectPeakAmpMean_condB = 0;
    subjectPeakAmpSD_condB = 290;
    
    % Specify number of unique actors based on input argument
    actorN = length(sampleActorPeakAmpArray);

    % Specify electrode to keep 
    keepERPElec = 2; 

    % Specify unneeded channels as a string for ERPLAB input.
    % In this example, we simulate 2 channels and the channel we want (C4)
    % is channel 2 so we delete the first channel. This step is important
    % because the real noise file contains data for one channel (C4)
    deleteERPChans = 'delerpchan( 1)';
    
%% 2. DEFINE FUNCTION FOR GENERATING ERP SIGNAL CLASS

    % Anonymous function to generate an ERP "signal class" with the specified
    % peak amplitude, latency, and slope. The signal class is then used by the
    % generate_scalpdata function to simulate trial-level waveforms.  
    % - Inputs:
        % - overallPeakAmp: Peak amplitude summed over the peak amplitude values
        %   for emotion condition population mean, actor intercept, subject intercept, 
        %   and maternal sensitivity/age intercept.  
        % - overallPeakLatency: Peak latency summed over the default peak
        %   latency for this component (i.e., emotionPeakLatency) and the subject's
        %   peak latency shift (i.e., a value drawn from the subjectPeakLatencyShiftArray).
	    % - emotionPeakAmpSlope: Decay rate for the peak amplitude of the specified
	    %   emotion condition. 
        % - Other variables (e.g., emotionPeakWindow) are defined above and held
        %   constant for all simulated subjects. 
    % - Output: 
        % - ERP signal class with the specified parameters. 
    erp = @(overallPeakAmp, overallPeakLatency, emotionPeakAmpSlope) ...
        utl_check_class(struct( ...
        'type', 'erp', ...
        'peakLatency', overallPeakLatency, ...
        'peakWidth', emotionPeakWindow, ...
        'peakAmplitude', overallPeakAmp, ...
        'peakAmplitudeSlope', emotionPeakAmpSlope, ...
        'peakAmplitudeDv', emotionPeakAmpDv, ...
        'peakLatencyDv', 0, ...
        'peakLatencyShift', 0));
   
    
%% 3. GENERATE TRIAL-LEVEL WAVEFORMS WITH PARAMETERS AND FUNCTION FROM STEPS 1-2

    % Generate subject intercept for condition B from a normal
    % distribution with above parameters. There is no subject intercept
    % added to condition A.
    subjectPeakAmpIncrement_condA = 0;
    subjectPeakAmpIncrement_condB = normrnd(subjectPeakAmpMean_condB,subjectPeakAmpSD_condB);
    subjectPeakAmpIncrementArray = [subjectPeakAmpIncrement_condA, subjectPeakAmpIncrement_condB];
    
    % Generate peak amplitude array with one value for each unique stimulus
    % (emotion condition/actor). Each value consists of the sum of the emotion
    % condition population mean, actor intercept, subject intercept, and maternal
    % sensitivity/age intercept. (NOTE: While the emotion condition and actor intercept
    % varies by stimulus, the subject and maternal sensitivity/age intercepts do not.) 
        % - For example, the first value of overallPeakAmpArray is the sum
        %   of the peak amplitudes for emotion condition A + actor 01 + subject
        %   intercept + maternal sensitivity/age intercept. The second value is the sum of 
        %   condition A + actor 02 + subject intercept + maternal sensitivity/age intercept, and so on. 
        % - Variables such as emotionPeakAmpSlope and emotionPeakAmpDv are taken into
        %   account separately with the erp function.
    overallPeakAmpArray = repelem(emotionPeakAmpArray,actorN) + repelem(subjectPeakAmpIncrementArray,actorN) + repmat(sampleActorPeakAmpArray, 1, length(emotionPeakAmpArray));
    
    % Randomly select this subject's peak latency shift from a uniform 
    % distribution with the array specified above and add it to the 
    % emotionPeakLatency variable. The overallPeakLatency variable is used to
    % define the peak latency of this subject's waveform.
    overallPeakLatency = emotionPeakLatency + randi(subjectPeakLatencyShiftArray,1,1); 

    % Simulate the first stimulus' 10 trial-level waveforms (corresponding to
    % 10 presentations/stimulus). The waveform's amplitude reduces with 
    % each unique presentation based on the emotionPeakAmpSlope_condA variable. 
    v = 1; % Counter variable indexing the overallPeakAmpArray variable
    % Define a component consisting of a neural source, dipole orientation
    % and ERP signal
    componentTemp1 = struct('source', sourceLocs, ... 
        'signal', {{erp(overallPeakAmpArray(v),overallPeakLatency,emotionPeakAmpSlope_condA)}});
    % Validate the component structure
    componentTemp1 = utl_check_component(componentTemp1, leadField);
    % Simulate the trial-level waveforms with the specified component,
    % leadfield, and epoch structure
    dataTemp1 = generate_scalpdata(componentTemp1, leadField, epochs, 'showprogress', 0);
    % Format the simulated data as an EEGLAB .set file and add the
    % corresponding event marker preceding code (e.g., 601)
    EEGTemp1 = utl_create_eeglabdataset(dataTemp1, epochs, leadField, ...
        'marker', convertStringsToChars(emotionPrecCodes(v)));

    % Repeat this process for the next unique stimulus until trial-level
    % waveforms for all 10 unique stimuli have been simulated
    v = 2;
    componentTemp2 = struct('source', sourceLocs, ...
        'signal', {{erp(overallPeakAmpArray(v),overallPeakLatency,emotionPeakAmpSlope_condA)}});
    componentTemp2 = utl_check_component(componentTemp2, leadField);
    dataTemp2 = generate_scalpdata(componentTemp2, leadField, epochs, 'showprogress', 0);
    EEGTemp2 = utl_create_eeglabdataset(dataTemp2, epochs, leadField, ...
        'marker', convertStringsToChars(emotionPrecCodes(v)));
    EEG1 = pop_mergeset(EEGTemp1,EEGTemp2); % Merge the two EEG .set files together (note that pop_mergeset only accepts two input arguments at a time)

    v = 3;
    componentTemp3 = struct('source', sourceLocs, ...
        'signal', {{erp(overallPeakAmpArray(v),overallPeakLatency,emotionPeakAmpSlope_condA)}});
    componentTemp3 = utl_check_component(componentTemp3, leadField);
    dataTemp3 = generate_scalpdata(componentTemp3, leadField, epochs, 'showprogress', 0);
    EEGTemp3 = utl_create_eeglabdataset(dataTemp3, epochs, leadField, ...
        'marker', convertStringsToChars(emotionPrecCodes(v)));
    EEG2 = pop_mergeset(EEG1,EEGTemp3);

    v = 4;
    componentTemp4 = struct('source', sourceLocs, ...
        'signal', {{erp(overallPeakAmpArray(v),overallPeakLatency,emotionPeakAmpSlope_condA)}});
    componentTemp4 = utl_check_component(componentTemp4, leadField);
    dataTemp4 = generate_scalpdata(componentTemp4, leadField, epochs, 'showprogress', 0);
    EEGTemp4 = utl_create_eeglabdataset(dataTemp4, epochs, leadField, ...
        'marker', convertStringsToChars(emotionPrecCodes(v)));
    EEG3 = pop_mergeset(EEG2,EEGTemp4);

    v = 5;
    componentTemp5 = struct('source', sourceLocs, ...
        'signal', {{erp(overallPeakAmpArray(v),overallPeakLatency,emotionPeakAmpSlope_condA)}});
    componentTemp5 = utl_check_component(componentTemp5, leadField);
    dataTemp5 = generate_scalpdata(componentTemp5, leadField, epochs, 'showprogress', 0);
    EEGTemp5 = utl_create_eeglabdataset(dataTemp5, epochs, leadField, ...
        'marker', convertStringsToChars(emotionPrecCodes(v)));
    EEG4 = pop_mergeset(EEG3,EEGTemp5);

    % The following trials are simulated for emotion condition B so we
    % specify the slope as emotionPeakAmpSlope_condB
    v = 6;
    componentTemp6 = struct('source', sourceLocs, ...
        'signal', {{erp(overallPeakAmpArray(v),overallPeakLatency,emotionPeakAmpSlope_condB)}});
    componentTemp6 = utl_check_component(componentTemp6, leadField);
    dataTemp6 = generate_scalpdata(componentTemp6, leadField, epochs, 'showprogress', 0);
    EEGTemp6 = utl_create_eeglabdataset(dataTemp6, epochs, leadField, ...
        'marker', convertStringsToChars(emotionPrecCodes(v)));
    EEG5 = pop_mergeset(EEG4,EEGTemp6);

    v = 7;
    componentTemp7 = struct('source', sourceLocs, ...
        'signal', {{erp(overallPeakAmpArray(v),overallPeakLatency,emotionPeakAmpSlope_condB)}});
    componentTemp7 = utl_check_component(componentTemp7, leadField);
    dataTemp7 = generate_scalpdata(componentTemp7, leadField, epochs, 'showprogress', 0);
    EEGTemp7 = utl_create_eeglabdataset(dataTemp7, epochs, leadField, ...
        'marker', convertStringsToChars(emotionPrecCodes(v)));
    EEG6 = pop_mergeset(EEG5,EEGTemp7);

    v = 8;
    componentTemp8 = struct('source', sourceLocs, ...
        'signal', {{erp(overallPeakAmpArray(v),overallPeakLatency,emotionPeakAmpSlope_condB)}});
    componentTemp8 = utl_check_component(componentTemp8, leadField);
    dataTemp8 = generate_scalpdata(componentTemp8, leadField, epochs, 'showprogress', 0);
    EEGTemp8 = utl_create_eeglabdataset(dataTemp8, epochs, leadField, ...
        'marker', convertStringsToChars(emotionPrecCodes(v)));
    EEG7 = pop_mergeset(EEG6,EEGTemp8);

    v = 9;
    componentTemp9 = struct('source', sourceLocs, ...
        'signal', {{erp(overallPeakAmpArray(v),overallPeakLatency,emotionPeakAmpSlope_condB)}});
    componentTemp9 = utl_check_component(componentTemp9, leadField);
    dataTemp9 = generate_scalpdata(componentTemp9, leadField, epochs, 'showprogress', 0);
    EEGTemp9 = utl_create_eeglabdataset(dataTemp9, epochs, leadField, ...
        'marker', convertStringsToChars(emotionPrecCodes(v)));
    EEG8 = pop_mergeset(EEG7,EEGTemp9);

    v = 10;
    componentTemp10 = struct('source', sourceLocs, ...
        'signal', {{erp(overallPeakAmpArray(v),overallPeakLatency,emotionPeakAmpSlope_condB)}});
    componentTemp10 = utl_check_component(componentTemp10, leadField);
    dataTemp10 = generate_scalpdata(componentTemp10, leadField, epochs, 'showprogress', 0);
    EEGTemp10 = utl_create_eeglabdataset(dataTemp10, epochs, leadField, ...
        'marker', convertStringsToChars(emotionPrecCodes(v)));
    EEG = pop_mergeset(EEG8,EEGTemp10); 

    EEG = epoch2continuous(EEG); % Concatenate the epoched dataset into a continuous dataset

%% 4. UPDATE EVENT MARKER PRECEDING CODES WITH TRIAL PRESENTATION NUMBER

    allEventTypes = {EEG.event.type}'; % Extract the subject's event markers (consisting of 3-digit preceding codes) 
    nonBoundaryEventIdx = ~strcmp(allEventTypes, "boundary"); % Locate non-boundary event markers (i.e., event markers from the emotionPrecCodes array)
    % Update the event marker name with the trial presentation number (e.g., convert
    % the first instance of 601 to 60101; convert the second instance of 601 to
    % 60102, etc.)
    allEventTypes(nonBoundaryEventIdx) = strcat(allEventTypes(nonBoundaryEventIdx), cellstr(presentNArray));
    [EEG.event.type] = deal(allEventTypes{:}); % Update the subject's event marker array with the final 5-digit event markers
    EEG = eeg_checkset(EEG, 'eventconsistency'); % Check EEG event array for inconsistencies (e.g., event markers out of order) 

%% 5. EXTRACT BIN-BASED EPOCHS

    % Create EventList
    EEG  = pop_creabasiceventlist(EEG, 'AlphanumericCleaning', 'on', 'BoundaryNumeric', {-99}, 'BoundaryString', {'boundary'});

    % Assign events to bins
    EEG  = pop_binlister(EEG, 'BDF', binDescriptorFilename, 'IndexEL',  1, 'SendEL2', 'EEG', 'UpdateEEG', 'off', 'Voutput', 'EEG','Report','off');

    % Extract bin-based epochs and baseline correct
    EEG = pop_epochbin(EEG, [-preStimPeriod  postStimPeriod],  'pre');

%% 6. CALCULATE TRIAL-LEVEL ERPS

    ERP = pop_averager(EEG , 'Criterion', 'good', 'DQ_flag', 1, 'ExcludeBoundary', 'off', 'SEM', 'on');

%% 7. DELETE EXTRA SIMULATED ELECTRODES 

    % We delete C3 (extra simulated electrode) from the data file so that 
    % the data file has only C4. This makes it easier to append the C4 data
    % when we add the simulated ERPset with the real noise file.
    ERP = pop_erpchanoperator( ERP, { deleteERPChans} , 'ErrorMsg', 'popup', 'KeepLocations',  0, 'Warning', 'on' );
    
    % Save binerror for final bin channel only
    ERP.binerror = ERP.binerror(keepERPElec,:,:); 
end