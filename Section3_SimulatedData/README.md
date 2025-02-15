This folder contains pipeline scripts used to simulate data for Section 3 of Heise et al. (submitted). We simulate a total of 1,000 datasets (samples). For each dataset, we simulate 48 subjects, and each simulated subject saw 100 trials (2 emotion conditions and 5 different ‘actors’ with 10 presentations each). 

## Table of Contents  
* [Script overview](#script-overview) 
* [Script requirements](#script-requirements)
  
## Script overview

*MATLAB scripts for simulating trial-level ERP data*
* **DiffWaveSim_01_CreateBinDescriptorFile.m**: Creates a bin descriptor file used to extract trial-level bins (i.e., each bin corresponds to a presentation of a specific stimulus) in subsequent scripts (DiffWaveSim_02_SimulateERPData.m, simulateOneSubject_forDiffWave.m). 
* **DiffWaveSim_02_SimulateERPData.m**: Simulates trial-level ERP data using two helper functions (simulateOneSample_forDiffWave.m, and simulateOneSubject_forDiffWave.m).

*R scripts for organizing data files and fitting LME and ANOVA models*
* **DiffWaveSim_03_OrganizeDataFiles.R**: Prepares simulated data files for analysis by extracting trial and subject-related information.
* **DiffWaveSim_04_ExtractModelOutput.R**: Induces missing trials/subjects and analyzes data using seven difference wave approaches using two helper functions (DiffWaveSim_04_funs.R and setNSFunctions.R).
* **DiffWaveSim_05_CompareModelPerformance.R**: Calculates the number of usable trials and subjects, statistical power to detect the maternal sensitivity effect, and error and bias of the maternal sensitivity estimates for each difference wave approach. 

## Script requirements
* MATLAB R2019a: [https://www.mathworks.com/](https://www.mathworks.com/)
* EEGLAB v. 2019_0: [https://sccn.ucsd.edu/eeglab/index.php](https://sccn.ucsd.edu/eeglab/index.php)
* ERPLAB v. 8.01: [https://erpinfo.org/erplab/](https://erpinfo.org/erplab/)
* SEREEGA v. 1.1.0: [https://github.com/lrkrol/SEREEGA](https://github.com/lrkrol/SEREEGA)
* Pediatric Head Atlas release v. 1.1, Atlas 1 (0-2 years old): [https://www.pedeheadmod.net/pediatric-head-atlases/](https://www.pedeheadmod.net/pediatric-head-atlases/)
  * After requesting and downloading the atlas, the atlas folder should be added to the MATLAB path (via "Home" > "Set Path" > "Add with Subfolders").
* Real infant noise files extracted from a resting state task (see Appendix A of Heise et al., submitted)
* R v. 3.6.1: [https://www.r-project.org/](https://www.r-project.org/)
* lme4 v. 1.1-25: [https://cran.r-project.org/web/packages/lme4/index.html](https://cran.r-project.org/web/packages/lme4/index.html)
* lmerTest v. 3.1-3: [https://cran.r-project.org/web/packages/lmerTest/index.html](https://cran.r-project.org/web/packages/lmerTest/index.html)
* MatchIt v. 4.5.0: [https://cran.r-project.org/web/packages/MatchIt/index.html](https://cran.r-project.org/web/packages/MatchIt/index.html)
* afex v. 0.28-1: [https://cran.r-project.org/web/packages/afex/index.html](https://cran.r-project.org/web/packages/afex/index.html)
* emmeans v. 1.5.3: [https://cran.r-project.org/web/packages/emmeans/index.html](https://cran.r-project.org/web/packages/emmeans/index.html)
* See R scripts for other required packages 
