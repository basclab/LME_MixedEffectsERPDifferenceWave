% Difference Wave Simulation Helper Function: Simulate One Data Sample

% Purpose: This function simulates a trial-level ERP data file for each subject
% in a sample. Then, all subjects' ERP files are concatenated into one file
% and the mean amplitude for the electrode of interest (C4) and time window 
% (300-500 ms) is extracted for each subject and trial-specific bin. In
% addition, a subject data log file is created, which lists each subject's
% assigned maternal sensitivity and age group. 

% Notes:
% - Actor intercepts are added to the peak amplitude for both emotion conditions 
%   (condition A/Neutral and condition B/Happy).
% - Subjects are assigned to a maternal sensitivity group (low or high) and
%   an age group (younger or older). These two group assignments determine a 
%   subject's intercept that is added to their peak amplitude for condition B (Happy). 
% - Peak amplitude values were chosen to produce a specific mean amplitude 
%   value when extracted from C4 over a 300-500 ms time window. For example, 
%   the actor peak amplitude intercept of 581 �V corresponds to a mean amplitude of 
%   -10.002 �V over a 300-500 ms time window at C4. 

% ***See Section3_SimulatedData README.md available on the LME_MixedEffectsERPDifferenceWave
% GitHub for pipeline details: https://github.com/basclab/LME_MixedEffectsERPDifferenceWave/tree/main/Section3_SimulatedData

% Format:
    % simulateOneSample_forDiffWave(sample, subjectN, decayRate, saveFolder, leadField, sourceLocs, noiseERPArray)

% Inputs:
    % - sample: Sample ID for the current simulated sample (e.g., if sample = 1,
    %   the current sample will have an ID of #1).
    % - subjectN: Number of simulated subjects per sample.
    % - decayRate: Decay rate type used to specify whether a 'different' or 
    %   'same' decay rate is simulated across the two emotion conditions. 
    % - saveFolder: Folder for saving simulated data output files.
    %   This parent folder has two subfolders: 01_MeanAmpOutput_PreMerge
    %   and 01_SubjectDataLog, which are used for saving the corresponding
    %   mean amplitude output files and subject data logs (see
    %   Outputs section below for more information).  
    % - leadField: Data structure created in 
    %   DiffWaveSim_02_SimulateERPData.m for specifying the lead field, 
    %   electrode montage, electrodes of interest, and dipole orientation.
    % - sourceLocs: Index used to identify the dipole location from the
    %   leadField structure. This variable was created in
    %   DiffWaveSim_02_SimulateERPData.m.
    % - noiseERPArray: EEGLAB array of .erp files containing resting state
    %   trials extracted from real participants. The epoch window length
    %   must be identical to the length of the simulated data in order to
    %   combine them together.
    
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
        
% Function Steps:
    % 1. Specify output filenames and other variables
    % 2. Define simulation parameters
    % 3. Simulate each subject's ERP file
    % 4. Export this sample's mean amplitude output file 
    % 5. Export this sample's subject data log 
    % 6. (Optional) Export this sample's .erp file

% Outputs:
    % - Mean amplitude .txt files with one mean amplitude value per bin/
    %   electrode/subject. There is one file for each simulated sample.
    %   Each file contains the following columns:
        % - value: The mean amplitude for this simulated subject�s specified
        %   bin and electrode. This value is extracted from C4 (corresponding 
        %   to E104 of the EGI HydroCel GSN 128-electrode montage) over a 
        %   300-500 ms time window. 
        % - chindex: The electrode number (e.g., electrode #1 corresponds to 
        %   E104 in this simulation).
        % - chlabel: The label for this electrode (e.g., E104).
        % - bini: A number generated by ERPLAB when concatenating ERP
        %   data files across each sample's subjects. NOTE: This number
        %   is NOT used for subsequent processing. 
        % - binlabel: A label composed of [SUBJECTID]_:_[trial-specific bin label].
        %   The trial-specific bin label corresponds to the labels from the
        %   bin descriptor file (created in DiffWaveSim_01_CreateBinDescriptorFile.m).
        %   The binlabel column is used to identify subject ID and stimuli-related  
        %   information in DiffWaveSim_03_OrganizeDataFiles.R. 
        % - ERPset: This column is intentionally empty because it is NOT needed 
        %   for identifying subject ID (see binlabel column above). 
    % - Subject data log .txt files with one row corresponding to each subject.
    %   There is one file for each simulated sample. Each file contains the 
    %   following columns:
        % - SUBJECTID: Simulated subject ID (e.g., 01, 02, �).
        % - mSens_age: String composed of a subject's assigned [maternal sensitivity group]_[age group]
        %   (e.g., lowMSens_youngerAgeGroup). 
    % - (Optional) .erp files containing the trial-level waveforms for all
    %   subjects in a sample. There is one file for each simulated sample
    %   and they are saved directly in the saveFolder specified above (not
    %   in a subfolder). These files are useful for visualizing waveforms
    %   or troubleshooting.  

% Usage Example:
    % >> sample = 1;
    % >> subjectN = 48;
    % >> decayRate = 'different';
    % >> saveFolder = 'C:\Users\basclab\Desktop\Section3_SimulatedData';
    % >> leadField = lf_generate_frompha('0to2','128','labels',{'E36','E104'});
    % >> sourceLocs = lf_get_source_nearest(leadField, [10 46 18]);
    % >> leadField.orientation(sourceLocs,:) = [0.57 -0.70 -0.01];
    % >> [tempERP noiseERPArray] = pop_loaderp( 'filename', {'InfantNoise01.erp','InfantNoise02.erp','InfantNoise03.erp','InfantNoise04.erp','InfantNoise05.erp','InfantNoise06.erp','InfantNoise07.erp','InfantNoise08.erp','InfantNoise09.erp','InfantNoise10.erp'}, 'filepath', 'C:\Users\basclab\Desktop\Section3_SimulatedData\00_InfantNoiseRepository\');
    % >> simulateOneSample_forDiffWave(sample, subjectN, decayRate, saveFolder, leadField, sourceLocs, noiseERPArray);
    
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

function simulateOneSample_forDiffWave(sample, subjectN, decayRate, saveFolder, leadField, sourceLocs, noiseERPArray)
%% 1. SPECIFY OUTPUT FILENAMES AND OTHER VARIABLES

    % Specify filenames for output files: mean amplitude file, subject data
    % log, and (optional) .erp file. The mean amplitude filename includes 
    % 'PreMerge' to signify that the file does not contain age group
    % information yet (this information is added by the
    % DiffWaveSim_03_OrganizeDataFiles.R script). 
    saveSampleFilename = strcat('Sample', sprintf('%04d',sample)); % Format sample ID with leading zeros
    saveSampleMeanAmpFilename = strcat(saveFolder, '\01_NCMeanAmpOutput_PreMerge\', saveSampleFilename, '-NCMeanAmpOutput_PreMerge.txt'); 
    saveSampleSubjectDataLogFilename = strcat(saveFolder, '\01_SubjectDataLog\', saveSampleFilename, '-SubjectDataLog.txt');
    %saveSampleERPFilename = strcat(saveSampleFilename, '.erp');
    
    % Create ALLERP structure for storing all subjects' ERP files
    ALLERP = buildERPstruct([]);

    % Define variables for mean amplitude extraction: 
    % In this simulation, the NC mean amplitude is extracted over a 300-500 ms
    % time window for C4 (electrode #1).
    timeWindowArray = [300 500];
    electrodeArray = [1];

   %% 2. DEFINE SIMULATION PARAMETERS 

    % Define parameters for actor intercepts:
    % Each actor intercept is drawn from a normal distribution with one of
    % the following mean and standard deviation values (e.g., actor 01's 
    % normal distribution has a mean of 581 and standard deviation of 290;
    % actor 02's has a mean of 290 and standard deviation of 290, etc.).
    actorPeakAmpArray = [581, 290, 0, -290, -581]; % Peak amplitude population mean for each actor intercept distribution (one distribution per actor; values correspond to mean amplitudes of [-10.002, -4.993, 0, 4.993, 10.002], respectively, in units of �V)
    actorPeakAmpSD = 290; % Peak amplitude population standard deviation for all actor intercept distributions (value corresponds to mean amplitude of 4.993 �V)
    
    % Generate actor intercepts for this sample by randomly selecting one value
    % per actor. Each subject in this sample will have the same 5 actor intercepts.
    % These values are added to the subject's peak amplitude for both 
    % emotion conditions in simulateOneSubject_forDiffWave function below. 
    sampleActorPeakAmpArray = [normrnd(actorPeakAmpArray(1),actorPeakAmpSD), normrnd(actorPeakAmpArray(2),actorPeakAmpSD), ...
        normrnd(actorPeakAmpArray(3),actorPeakAmpSD), normrnd(actorPeakAmpArray(4),actorPeakAmpSD), ...
        normrnd(actorPeakAmpArray(5),actorPeakAmpSD)];

    % Define parameters for maternal sensitivity/age group intercept:
    % - Each subject is assigned to a maternal sensitivity group (i.e., lowMSens 
    %   or highMSens) and an age group (i.e., youngerAgeGroup or
    %   olderAgeGroup) such that each group has an approximately 
    %   equal number of subjects (see the subjectNSubgroup variable below).
    % - These group intercept values are constant for all simulated samples 
    %   (e.g., a subject in the lowMSens_olderAgeGroup in sample 1 will have 
    %   the same age intercept as a lowMSens_olderAgeGroup subject in sample 3). 
    % - Each subject's group intercept is added to their peak amplitude for 
    %   condition B in simulateOneSubject_forDiffWave function below. 
    mSensAgeLabelArray = ["lowMSens_olderAgeGroup"; "higMSens_olderAgeGroup"; "lowMSens_youngerAgeGroup"; "higMSens_youngerAgeGroup"]; % Array of age group labels
    emotionPeakAmpArray_condB = [581, 813, 349, 581]; % Peak amplitude group intercepts (value corresponds to mean amplitudes of [-10.002, -13.992, -6.008, -9.998], respectively, in units of �V)

    subjectNSubgroup = ceil(subjectN/length(mSensAgeLabelArray)); % Number of subjects assigned to each age group (this number is rounded up if subjectN is not evenly divisible by the number of age groups)
    sampleMSensAgeLabelArray_original = repelem(mSensAgeLabelArray, subjectNSubgroup); % Create group intercept array (NOTE: This array is not randomized yet and may have a length greater than subjectN)
    sampleMSensAgeLabelArray_random = datasample(sampleMSensAgeLabelArray_original, subjectN, 'Replace', false); % Extract randomized age intercept array of length subjectN 
    
    % Assign a real noise file to subject. If there are more simulated subjects  
    % than noise files, multiple subjects may be assigned the same noise file 
    % but the trial-level pairing of simulated and real noise trials will 
    % differ across subjects due to randomization below.
    noiseFileN = length(noiseERPArray);
    sampleNoiseNSubgroup = ceil(subjectN/noiseFileN); % Number of subjects assigned to each noise file (this number is rounded up if subjectN is not evenly divisible by the number of noise files)
    sampleNoiseIndexArray_original = repelem(1:noiseFileN, sampleNoiseNSubgroup); % Create noise filename array (NOTE: This array is not randomized yet and may have a length greater than subjectN)
    sampleNoiseIndexArray_random = datasample(sampleNoiseIndexArray_original, subjectN, 'Replace', false); % Extract randomized noise filename array of length subjectN 
    
%% 3. SIMULATE EACH SUBJECT'S ERP FILE

    for subject = 1:subjectN % Loop through each simulated subject
        subjectMSensAgeIndex = find(mSensAgeLabelArray == sampleMSensAgeLabelArray_random(subject));
        sampleEmotionPeakAmp_condB_oneSubject = emotionPeakAmpArray_condB(subjectMSensAgeIndex); % Extract this subject's peak amplitude for condition B 
        
        subjectALLERP = buildERPstruct([]); % Create empty ALLERP structure to store subject's simulated ERP and real noise ERP
        subjectALLERP(1) = simulateOneSubject_forDiffWave(sampleEmotionPeakAmp_condB_oneSubject, decayRate, sampleActorPeakAmpArray, leadField, sourceLocs); % Generate this subject's simulated ERP data file 
 
        % Load this subject's assigned real noise file
        subjectALLERP(2) = noiseERPArray(sampleNoiseIndexArray_random(subject));

        % Append the two ERP files
        subjectERP_append = pop_appenderp( subjectALLERP , 'Erpsets', [ 1 2]);
        % Extract simulated bin numbers
        binArray_simul = find(~contains(subjectERP_append.bindescr, 'Real noise'));
        binArray_noise = find(contains(subjectERP_append.bindescr, 'Real noise'));
        
        % Randomly sort binArray_noise
        binArray_noiseRandom = datasample(binArray_noise, length(binArray_noise),'Replace', false);
        % Extract bin labels from simulated data file
        binLabelArray = subjectALLERP(1).bindescr;
        % Create bin equations with randomly paired simulated trials and
        % real noise trials
        subjectBinEquations = strcat("nb", string(binArray_simul), " = b", string(binArray_simul), " + b", string(binArray_noiseRandom), " label ", binLabelArray);
        subjectERP_final = pop_binoperator( subjectERP_append, cellstr(subjectBinEquations));
     
        ALLERP(subject) = subjectERP_final; % Store this subject's final ERP file in the ALLERP structure
    end

%% 4. EXPORT THIS SAMPLE'S MEAN AMPLITUDE OUTPUT FILE

    % Concatenate all subjects' ERP data files into one ERP file
    ERP = pop_appenderp(ALLERP, 'Erpsets', [1:subjectN],'Prefixes',cellstr(pad(string(1:subjectN),2,'left','0')));
    
    % Calculate the mean amplitude for each bin over the time window and
    % electrodes specified above and save in .txt file
    ERP_meanAmp = pop_geterpvalues(ERP, timeWindowArray, [1:ERP.nbin], electrodeArray , ...
        'Baseline', 'pre', 'Binlabel', 'on', 'FileFormat', 'long', ...
        'Filename', saveSampleMeanAmpFilename, 'Fracreplace', 'NaN', 'InterpFactor',  1, ...
        'Measure', 'meanbl', 'PeakOnset',  1, 'Resolution',  3);
    
%% 5. EXPORT THIS SAMPLE'S SUBJECT DATA LOG

    % Create string array of all subject IDs
    SUBJECTID = pad(string(1:subjectN),2,'left','0')';
    
    subjectDataLog = table(SUBJECTID, sampleMSensAgeLabelArray_random); % Combine subject IDs and assigned groups into one table
    subjectDataLog.Properties.VariableNames([2]) = {'mSens_age'}; % Update group column name 
    writetable(subjectDataLog,saveSampleSubjectDataLogFilename); % Save table as .txt file
    
%% 6. (OPTIONAL) EXPORT THIS SAMPLE'S .ERP FILE
    
    % Uncomment the line below to save the concatenated .erp file containing
    % all subjects' trial-level waveforms    
    %ERP = pop_savemyerp(ERP, 'erpname', saveSampleERPFilename, 'filename', saveSampleERPFilename, 'filepath', saveFolder);

end