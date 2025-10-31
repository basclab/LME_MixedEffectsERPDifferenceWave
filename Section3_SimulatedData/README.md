This folder contains pipeline scripts used to simulate realistic Negative Central (NC) waveforms for Section 3 of Heise et al. (2025). We simulated two decay rate conditions: 'different decay' and 'same decay'. For each decay condition, we simulated a total of 1,000 datasets (samples). For each dataset, we simulated 48 subjects, and each simulated subject saw 100 trials (2 emotion conditions and 5 different ‘actors’ each presented 10 times per emotion condition). We include example output files for 2 datasets from the different decay rate condition.

## Table of Contents  
* [Script overview](#script-overview) 
* [Script requirements](#script-requirements)
  
## Script overview

*MATLAB scripts for simulating trial-level ERP data*
* **DiffWaveSim_01_CreateBinDescriptorFile.m**: Creates a bin descriptor file using the LMESimulation_EventMarkerMappingKey.xlsx file. The bin descriptor file is used to extract trial-level bins (i.e., each bin corresponds to a presentation of a specific emotion condition and actor) in subsequent scripts (DiffWaveSim_02_SimulateERPData.m, simulateOneSubject_forDiffWave.m). 
* **DiffWaveSim_02_SimulateERPData.m**: Simulates trial-level ERP data using two helper functions (simulateOneSample_forDiffWave.m, and simulateOneSubject_forDiffWave.m). More information about the simulation parameters are included in Appendix B of Heise et al. (2025).

*R scripts for organizing data files and fitting LME and ANOVA models*
* **DiffWaveSim_03_OrganizeDataFiles.R**: Prepares simulated data files for analysis by extracting trial and subject-related information.
* **DiffWaveSim_04_ExtractModelOutput.R**: Induces missing trials and subjects and analyzes data using seven difference wave approaches using two helper functions (DiffWaveSim_04_funs.R and setNSFunctions.R).
* **DiffWaveSim_05_CompareModelPerformance.R**: Calculates the number of usable trials and subjects, statistical power to detect the maternal sensitivity effect, and error and bias of the maternal sensitivity estimates for each difference wave approach. 

## Script requirements
* MATLAB R2019a: [https://www.mathworks.com/](https://www.mathworks.com/)
* EEGLAB v. 2019_0: [https://sccn.ucsd.edu/eeglab/index.php](https://sccn.ucsd.edu/eeglab/index.php)
* ERPLAB v. 8.01: [https://erpinfo.org/erplab/](https://erpinfo.org/erplab/)
* SEREEGA v. 1.1.0: [https://github.com/lrkrol/SEREEGA](https://github.com/lrkrol/SEREEGA)
* Pediatric Head Atlas release v. 1.1, Atlas 1 (0-2 years old): [https://www.pedeheadmod.net/pediatric-head-atlases/](https://www.pedeheadmod.net/pediatric-head-atlases/)
  * After requesting and downloading the atlas, the atlas folder should be added to the MATLAB path (via "Home" > "Set Path" > "Add with Subfolders").
* Real infant noise files: [https://osf.io/b53wj/](https://osf.io/b53wj/) > "Files" > "InfantNoise" > "ForNCSimulations"
* R v. 3.6.1: [https://www.r-project.org/](https://www.r-project.org/)
* lme4 v. 1.1-25: [https://cran.r-project.org/web/packages/lme4/index.html](https://cran.r-project.org/web/packages/lme4/index.html)
* lmerTest v. 3.1-3: [https://cran.r-project.org/web/packages/lmerTest/index.html](https://cran.r-project.org/web/packages/lmerTest/index.html)
* MatchIt v. 4.5.0: [https://cran.r-project.org/web/packages/MatchIt/index.html](https://cran.r-project.org/web/packages/MatchIt/index.html)
* afex v. 0.28-1: [https://cran.r-project.org/web/packages/afex/index.html](https://cran.r-project.org/web/packages/afex/index.html)
* emmeans v. 1.5.3: [https://cran.r-project.org/web/packages/emmeans/index.html](https://cran.r-project.org/web/packages/emmeans/index.html)
* See scripts for other required variables and packages
* See GitHub repository folder structure (00_InfantNoiseRepository, 01_ERPFiles, etc.) for required folders for saving output files
