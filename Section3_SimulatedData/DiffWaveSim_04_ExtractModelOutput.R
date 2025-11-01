# Difference Wave Simulation Script: 4. Extract Model Output

# This script imports each simulated mean amplitude file, induces missing 
# trials and subjects, and analyzes data using seven difference wave approaches:
  # - Six LME approaches: Interaction, Exact Match, three types of Nearest
  #   Neighbor, Random Permutation
  # - One conventional ANOVA approach using a common trial threshold (10 trials
  #   per condition)

# ***See Section3_SimulatedData README.md available on the LME_MixedEffectsERPDifferenceWave
# GitHub for pipeline details: https://github.com/basclab/LME_MixedEffectsERPDifferenceWave/tree/main/Section3_SimulatedData

# Requirements: 
  # - Needs R Version 3.6.1 and packages listed below
  # - importFolder: Folder containing the mean amplitude files created by 
  #   DiffWaveSim_03_OrganizeDataFiles.R
  # - See instructions in data environment and steps 1 and 2 for specifying other
  #   variables
  # - DiffWaveSim_04_funs.R and setNSFunctions.R are saved in the same working directory

# Script Functions:
  # 1. Specify missingness pattern 
  # 2. Specify other simulation parameters
  # 3. Load each simulated data file
  # 4. Fit difference wave approaches to population data (i.e., no missing data)
  # 5. Induce missing data
  # 6. Fit difference wave approaches after inducing missing data
  # 7. Save model output across all simulated datasets

# Outputs: 
  # - One .csv file containing the estimated marginal means for each maternal 
  #   sensitivity group and other model predictors across all models and simulated
  #   datasets for one type of decay and missingness pattern (e.g., different decay/
  #   Missingness Pattern #1). See ModelOutput_DataDictionary.xlsx file for further
  #   information.

# Copyright 2024 Megan J. Heise, Serena K. Mon, Lindsay C. Bowman
# Brain and Social Cognition Lab, University of California Davis, Davis, CA, USA.

# Permission is hereby granted, free of charge, to any person obtaining a 
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

library(data.table) # V.1.13.2; fread function
library(stringr) # V.1.4.0; str_sub function
library(car) # V.3.0-10; contr.sum function
library(lme4) # V.1.1-25; LME models
library(lmerTest) # V.3.1-3; p-values for LME
library(MatchIt) # V.4.5.0; nearest neighbor trial pairing  
library(afex) # V.0.28-1; ANOVA models
library(emmeans) # V.1.5.3; estimated marginal means calculation
library(performance) # V.0.6.1; check_convergence function
library(foreach) # V.1.5.2; foreach function
library(doParallel) # V.1.0.17; %dopar% operator
library(doRNG) # V.1.8.3; reproducible parallel results
RNGkind("L'Ecuyer-CMRG")

#-------------------------------------------------------------------------------
# DATA ENVIRONMENT

# Set up parallel processing to reduce processing time
numCore <- 7
clust <- makeCluster(numCore) 
registerDoParallel(clust)
showConnections()

# Add functions used for fitting LME and ANOVA models 
source("DiffWaveSim_04_funs.R") 

# Specify folder location of data files for analysis 
importFolder <- 'C:/Users/basclab/Desktop/Section3_SimulatedData/02_NCMeanAmpOutput_Final'

# Make directory of all .csv files in importFolder
fileDir <- list.files(path = importFolder, pattern = ".csv", full.names = TRUE, recursive = FALSE)
sampleN <- length(fileDir) # Number of simulated datasets (samples)
rpIter <- 10000 # Number of Random Permutation iterations for each sample

# Specify folder location to save output file
saveFolder <- 'C:/Users/basclab/Desktop/Section3_SimulatedData/03_ModelOutput'
missType <- 'Miss1' # Used for saveModelOutputName variable only
decayRate <- 'different' # Used for saveModelOutputName variable only

#-------------------------------------------------------------------------------
# 1. SPECIFY MISSINGNESS PATTERN 

# Specify probability weight distribution for trial presentation numbers 6-10. 
  # - The weight for numbers 1-5 is equal to 1 - presentNumberWeight6to10.
  # - For example, if presentNumberWeight6to10 = 0.7, then 70% of missing trials 
  #   are from trial presentation numbers 6-10 and 30% of trials are from trial
  #   presentation numbers 1-5.
  # - If both weight variables are equal to 0.5, an equal number of missing trials 
  #   are drawn from each trial presentation number (MCAR).
presentNumberWeight6to10 <- 0.7

# Specify probability weight distribution for Younger age group. 
  # - The weight for older age is equal to 1 - ageWeightyounger.
  # - For example, if ageWeightYounger = 0.7, then 70% of subjects selected for 
  #   more missing trials and subsequent casewise deletion are from the Younger 
  #   age group. 
ageWeightYounger <- 0.7

#-------------------------------------------------------------------------------
# 2. SPECIFY OTHER SIMULATION PARAMETERS

emotionTrialN <- 5*10 # Total number of trials per emotion condition (= number of actors*number of trial presentations)
emotionLabel <- c("A", "B") # Name of each emotion condition (used for induceMissingTrials_v2 function)
emotionN <- length(emotionLabel) # Number of emotion conditions
subjectN <- 48 # Total number of subjects 
presentNumberRef <- 5.5 # Average presentation number simulated in data (used for emmeans calculations)
trialMissingThreshold <- emotionTrialN - 10 # Max number of trials before subject is considered to have a low trial count

# Specify percent of subjects with low trial-count who will be casewise deleted
# for ANOVA model (e.g., 6 corresponds to 6% of subjects induced to have less 
# than 10 trials/condition remaining). 
caseDeletionPctArray <-  c(0, 6, 11, 32)
caseDeletionNArray <- ceiling((caseDeletionPctArray/100)*subjectN) # Corresponding number of subjects with low trial-count

# Specify LME Interaction model formula 
formulaLME_interact = meanAmpNC ~ emotion + presentNumber + mSens + age + mSens:emotion + age:emotion + presentNumber:emotion + (1|SUBJECTID) + (1|ACTOR)
# Specify formula for other LME models (Exact Match, Nearest Neighbor, Random Permutation)
formulaLME_withPresentNumberAvg = meanAmpNC_BMinusA ~ presentNumberAvg + mSens + age + (1|SUBJECTID)

# Load functions used to speed up model fitting
source("setNSFunctions.R")
origFuns <- setNSFunctions()

# Used to track processing duration
start_time <- Sys.time()

#-------------------------------------------------------------------------------
# Loop through each simulated sample using foreach loop for parallel processing

set.seed(20230416) # Specify seed for reproducible results
modelOutput <- foreach (i=1:sampleN, .packages=c('tidyr', 'data.table', 'performance', 'emmeans', 'stringr', 'car', 'MatchIt', 'lme4', 'lmerTest', 'afex')) %dorng% { # Specify R packages used within foreach loop
  setDTthreads(1) # Used for fread function 
  
  # Preallocate empty list to store model output. We specify the length as the 
  # number of analysis models fit to each sample (3 models fit at population level + 
  # 7 models fit for each percentage of low trial-count subjects) 
  modelOutput <- vector("list", 3+(length(caseDeletionPctArray)*7)) 

  #-----------------------------------------------------------------------------
  # 3. LOAD EACH SIMULATED DATA FILE 
  
  # Import mean amplitude file ("population" dataset)
  dfOriginal <- fread(fileDir[i]) 
  
  # Extract sample ID
  sampleID <- str_sub(fileDir[i],-24,-21)
  
  # Create probability weight columns for trial presentation number and age group 
  # based on values specified above
  dfOriginal$presentNumberWeight <- ifelse(dfOriginal$presentNumber > 5,
                                           (presentNumberWeight6to10/(emotionTrialN/2)), 
                                           ((1-presentNumberWeight6to10)/(emotionTrialN/2)))
  dfOriginal$ageWeight <- ifelse(dfOriginal$age == 'youngerAgeGroup', 
                                 ageWeightYounger/(subjectN/2), (1-ageWeightYounger)/(subjectN/2))
 
  # Extract age probability weights (one weight per subject) used for randomly 
  # selecting subjects with low trial-counts
  dfAgeWeight <- aggregate(ageWeight ~ SUBJECTID, dfOriginal, mean, 
                           na.action = na.omit)
  
  #-----------------------------------------------------------------------------
  # 4. FIT DIFFERENCE WAVE APPROACHES TO POPULATION DATA (I.E., NO MISSING DATA)
  # Note we do not fit LME Nearest Neighbor and Random Permutation models because 
  # all trials were available and data could be fit with the Exact Match model.

  # Fit population-level LME Interaction model
  modelOutput[[1]] <- allTrials_Int(dfOriginal, formulaLME_interact)
  modelOutput[[1]]$caseDeletionPct <- "Pop." # Record no missing data condition
  modelOutput[[1]]$matchMethod <- "Int" # Record analysis approach
  modelOutput[[1]]$sample <- sampleID # Record sample ID
  modelOutput[[1]]$sampleRP <- 0 # Placeholder variable (not used for models except for Random Permutation)
  
  # Fit population-level LME Exact Match model
  modelOutput[[2]] <- pairTrials_EM(dfOriginal, formulaLME_withPresentNumberAvg)
  modelOutput[[2]]$caseDeletionPct <- "Pop."
  modelOutput[[2]]$matchMethod <- "EM"
  modelOutput[[2]]$sample <- sampleID
  modelOutput[[2]]$sampleRP <- 0

  # Fit population-level ANOVA model
  modelOutput[[3]] <- fitANOVA(dfOriginal)
  modelOutput[[3]]$caseDeletionPct <- "Pop."
  modelOutput[[3]]$matchMethod <- "DW" # Record ANOVA 'difference wave' approach
  modelOutput[[3]]$sample <- sampleID
  modelOutput[[3]]$sampleRP <- 0
  
  #-----------------------------------------------------------------------------
  # For each percentage of subjects with low trial-count
  for (j in 1:length(caseDeletionNArray)) {
    #---------------------------------------------------------------------------
    # 5. INDUCE MISSING DATA
    
    # Randomly sample the subject IDs that will have low trial counts and (for 
    # ANOVA model only) will be casewise deleted
    subjectCaseDeletion <- sample(dfAgeWeight$SUBJECTID, caseDeletionNArray[j],
                                  replace = FALSE, prob=dfAgeWeight$ageWeight)
     
    # Use induceMissingTrials_v2 function to remove trials based on specified 
    # missingness (step 1) and subjectCaseDeletion variable
    dfMissing_noNA <- induceMissingTrials_v2(dfOriginal, subjectCaseDeletion)

    #---------------------------------------------------------------------------
    # 6. FIT DIFFERENCE WAVE APPROACHES AFTER INDUCING MISSING DATA 

    # Fit LME Interaction model 
    modelOutput[[(7*(j-1))+5]] <- allTrials_Int(dfMissing_noNA, formulaLME_interact)
    modelOutput[[(7*(j-1))+5]]$caseDeletionPct <- paste0(as.character(caseDeletionPctArray[j]),'%') # Record percentage of subjects with low trial-count
    modelOutput[[(7*(j-1))+5]]$matchMethod <-  "Int" 
    modelOutput[[(7*(j-1))+5]]$sample <- sampleID 
    modelOutput[[(7*(j-1))+5]]$sampleRP <- 0 
    
    # Fit LME Exact Match model
    modelOutput[[(7*(j-1))+6]] <- pairTrials_EM(dfMissing_noNA, formulaLME_withPresentNumberAvg)
    modelOutput[[(7*(j-1))+6]]$caseDeletionPct <- paste0(as.character(caseDeletionPctArray[j]),'%')
    modelOutput[[(7*(j-1))+6]]$matchMethod <-  "EM"
    modelOutput[[(7*(j-1))+6]]$sample <- sampleID
    modelOutput[[(7*(j-1))+6]]$sampleRP <- 0

    # Fit LME Nearest Neighbor: No Prioritization (NN:N) model
    modelOutput[[(7*(j-1))+7]] <- pairTrials_NN_N(dfMissing_noNA, formulaLME_withPresentNumberAvg)
    modelOutput[[(7*(j-1))+7]]$caseDeletionPct <- paste0(as.character(caseDeletionPctArray[j]),'%')
    modelOutput[[(7*(j-1))+7]]$matchMethod <-  "NN:N"
    modelOutput[[(7*(j-1))+7]]$sample <- sampleID
    modelOutput[[(7*(j-1))+7]]$sampleRP <- 0
    
    # Fit LME Nearest Neighbor: Presentation Number Prioritization (NN:P) model
    modelOutput[[(7*(j-1))+8]] <- pairTrials_NN_P(dfMissing_noNA, formulaLME_withPresentNumberAvg)
    modelOutput[[(7*(j-1))+8]]$caseDeletionPct <- paste0(as.character(caseDeletionPctArray[j]),'%')
    modelOutput[[(7*(j-1))+8]]$matchMethod <-  "NN:P"
    modelOutput[[(7*(j-1))+8]]$sample <- sampleID
    modelOutput[[(7*(j-1))+8]]$sampleRP <- 0
    
    # Fit LME Nearest Neighbor: Stimulus Feature Prioritization (NN:F) model
    modelOutput[[(7*(j-1))+9]] <- pairTrials_NN_F(dfMissing_noNA, formulaLME_withPresentNumberAvg)
    modelOutput[[(7*(j-1))+9]]$caseDeletionPct <- paste0(as.character(caseDeletionPctArray[j]),'%')
    modelOutput[[(7*(j-1))+9]]$matchMethod <-  "NN:F"
    modelOutput[[(7*(j-1))+9]]$sample <- sampleID
    modelOutput[[(7*(j-1))+9]]$sampleRP <- 0

    # Fit LME Random Permutation model
    modelOutput[[(7*(j-1))+10]] <- pairTrials_RP(dfMissing_noNA, formulaLME_withPresentNumberAvg, rpIter)
    modelOutput[[(7*(j-1))+10]]$caseDeletionPct <- paste0(as.character(caseDeletionPctArray[j]),'%')
    modelOutput[[(7*(j-1))+10]]$matchMethod <-  "RP"
    modelOutput[[(7*(j-1))+10]]$sample <- sampleID
    
    # Casewise delete subjects
    dfCaseDeletion <- dfMissing_noNA[!(dfMissing_noNA$SUBJECTID %in% subjectCaseDeletion),]
    
    # Fit ANOVA model to casewise deleted data
    modelOutput[[(7*(j-1))+11]] <- fitANOVA(dfCaseDeletion)
    modelOutput[[(7*(j-1))+11]]$caseDeletionPct <- paste0(as.character(caseDeletionPctArray[j]),'%')
    modelOutput[[(7*(j-1))+11]]$matchMethod <- "DW"
    modelOutput[[(7*(j-1))+11]]$sample <- sampleID
    modelOutput[[(7*(j-1))+11]]$sampleRP <- 0
  }
  modelOutput # This line ensures that the foreach loop will return the output variable 
}
modelOutputFinal <- do.call(rbind, do.call(Map, c(f = rbind, modelOutput))) # Bind model output across all samples into one dataframe

# Reset functions modified for speeding up model fitting
resetNSFunctions(origFuns)

# Calculate processing duration
end_time <- Sys.time()
end_time - start_time

#-------------------------------------------------------------------------------
# 7. SAVE MODEL OUTPUT ACROSS ALL SIMULATED DATASETS
saveModelOutputName <- paste0(saveFolder, '/ModelOutput_', missType, '_', decayRate, 'Decay_sampleN', sampleN, '_subN', subjectN, '_rpIter', rpIter, '.csv')
fwrite(modelOutputFinal, saveModelOutputName, row.names=FALSE)

# Close parallel cluster
stopCluster(clust)
showConnections()