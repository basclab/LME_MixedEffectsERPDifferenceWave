# Difference wave analysis of real preschooler ERP data

# This script was used to analyze the NC difference wave mean amplitudes presented 
# in Section 4 of Heise, Mon, and Bowman (submitted). The input files contain
# NC mean amplitudes for two emotion conditions: 
# Full-Intensity Angry and Reduced-Intensity Angry

# Our research question was: is the Full-minus-Reduced-Intensity Angry 
# NC difference wave modulated by a continuous predictor (perceptual sensitivity 
# subscale from the Early Childhood Behavioral Questionnaire-Short Form)?

# We compare several difference wave analysis approaches:
  # - Six LME approaches: Interaction, Exact Match, three types of Nearest
  #   Neighbor, Random Permutation
  # - Two conventional linear regression approaches using two common trial thresholds:
  #  10 trials per condition, 15 trials per condition
# See manuscript for more information about each analysis approach. 

# Requirements:
  # - Needs R Version 3.6.1 and packages listed below
  # - Input files containing trial-level and condition-level NC mean amplitudes
  #   (see comments below)
  # - setNSFunctions is saved in the same working directory 

# Script Functions:
  # 1. Define function for trial-level analysis 
  # 2. Load data files
  # 3. Format trial-level data for analysis
  # 4. Fit LME Interaction model
  # 5. Fit LME Exact Match model
  # 6. Fit LME Nearest Neighbor: No Prioritization model
  # 7. Fit LME Nearest Neighbor: Presentation Number Prioritization model
  # 8. Fit LME Nearest Neighbor: Stimulus Feature Prioritization model
  # 9. Fit LME: Random Permutation model
  # 10. Fit 10-Trial Regression model
  # 11. Fit 15-Trial Regression model
  # 12. Compare ECBQ perceptual sensitivity estimates across models

# Outputs: 
  # - Estimated ECBQ perceptual sensitivity predictor across models

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

# Load required packages
library(dplyr) # V.1.0.2; dataframe manipulation
library(car) # V.3.01-0; ontr.sum function 
library(rstatix) # V.0.6.0; regression assumption checks
library(tidyverse) # V.1.3.0; group_by function 
library(lme4) # V.1.1-25; lme models
library(lmerTest) # V.3.1-3; p-values for lme
library(MatchIt) # V.4.1.0; nearest neighbor trial pairing 
library(tidyr) # V.1.1.2; spread function
library(emmeans) # V.1.5.3; estimated marginal means calculation
library(performance) # V.0.6.1; check_convergence function
library(coxed) # V.0.3.3; bootstrapped CI calculation
library(foreach) # V.1.5.2; foreach function
library(doParallel) # V.1.0.17; %dopar%
library(doRNG) # V.1.8.3; reproducible parallel results
RNGkind("L'Ecuyer-CMRG") 

#-------------------------------------------------------------------------------
# 1. DEFINE FUNCTION FOR TRIAL-LEVEL ANALYSIS

# Load functions used for faster random permutation processing (same file used for
# simulation code)
source("setNSFunctions.R")

# makeWide_realData: Function adapted from LME simulation code to convert long dataframe 
# into wide to facilitate trial-level difference wave calculations. We changed 'A' labels to
# 'FullAngry' and 'B' labels to 'ReduAngry'.
# - Format:
#     dfOutput_wide <- makeWide_realData(dfInput_long)
# - Input: 
  # - dfInput_long: Trial-level dataframe in which rows have been sorted based on trial pairing ID
  #   and condition (e.g., row 1 is trial pair 1/Full Angry, row 2 is trial pair 1/Redu Angry,
  #   row 3 is trial pair 2/Full Angry, row 4 is trial pair 2/Redu Angry, etc.). 
# - Output: 
  # - dfOutput_wide: Trial-paired dataframe in wide format in which each condition
  #   has columns for mean amplitude and trial presentation number. 
makeWide_realData = function(dfInput_long) { 
  # Create array of indices for each pair of trials (extract 1st row of each pair)
  i <- seq(1, nrow(dfInput_long) - 1, by = 2)
  
  # For each 2nd row of each trial pair (corresponding to the second condition),
  # we extract mean amplitude and presentation number information and save as new columns
  # in the 1st row of each pair 
  dfInput_long[i, c("meanAmpNC.ReduAngry", "presentNumber.ReduAngry")] = dfInput_long[i+1, c("meanAmpNC", "presentNumber")]
  
  # Rename the original meanAmpNC and presentNumber columns to reflect the first condition 
  names(dfInput_long)[names(dfInput_long) %in% c("meanAmpNC")] <- "meanAmpNC.FullAngry"
  names(dfInput_long)[names(dfInput_long) %in% c("presentNumber")] <- "presentNumber.FullAngry"
  
  # Return the first row of each trial pair, which now contains the mean amplitude
  # and trial presentation number information for both conditions
  return(dfOutput_wide = dfInput_long[i, ])
}

#-------------------------------------------------------------------------------
# 2. LOAD DATA FILES 

# Data files are in long format:
# - For the trial-level data, each row corresponds to one participant, one electrode,
#   one trial presentation, and one condition.
# - For the condition-level data, each row corresponds to one participant and one 
#   condition (data has been averaged across electrode and trial presentations).
# Full-Intensity Angry is labelled as FullAngry and Reduced-Intensity Angry is
# labelled as ReduAngry.

setwd('C:/Users/basclab/Desktop/Section4_RealData')
dfTrial <- read.csv('nc_mAmp_trialLevel.csv')
dfCond_10TrialMin <- read.csv(nc_mAmp_condLevel_10TrialMin.csv)
dfCond_15TrialMin <- read.csv(nc_mAmp_condLevel_15TrialMin.csv)

#-------------------------------------------------------------------------------
# 3. FORMAT TRIAL-LEVEL DATA FOR ANALYSIS

# Format SUBJECTID and ACTOR columns as factors
dfTrial$SUBJECTID <- as.factor(dfTrial$SUBJECTID) 
dfTrial$ACTOR <- as.factor(dfTrial$ACTOR)

# Sum-code emotion column and specify factor ordering of FullAngry before ReduAngry 
# for subsequent processing with makeWide_realData function 
dfTrial$emotion <- factor(dfTrial$emotion, levels = c("FullAngry", "ReduAngry"))
contrasts(dfTrial$emotion)=contr.Sum(levels(dfTrial$emotion))

# Set Cz as the reference level for Interaction model analysis
dfTrial$electrode <- factor(dfTrial$electrode, levels = c("Cz","C3","C4"))

# Create column of 0's and 1's where FullAngry is set as treatment (1) for
# matchit function (Nearest Neighbor approaches)
dfTrial$emotionBinary <- ifelse(dfTrial$emotion == "FullAngry", 1, 0)

#-------------------------------------------------------------------------------
# 4. FIT LME INTERACTION MODEL

# Specify Interaction model formula where predictor of interest is
# emotion:ECBQ_per interaction
formulaLME_Int = meanAmpNC ~ emotion + presentNumber + electrode + ECBQ_per +  
  emotion:presentNumber + emotion:ECBQ_per + (1|SUBJECTID) + (1|ACTOR)
LME_Int <- lmer(formulaLME_Int, data=dfTrial, REML=TRUE)
summary(LME_Int)

# Check LME assumptions:
# - Linearity
plot(resid(LME_Int),dfTrial$meanAmpNC) 
# - Normal distribution of residuals
qqmath(LME_Int) 
# - Homoscedasticity (equal variance for each emotion condition)
leveneTest(residuals(LME_Int) ~ dfTrial$emotion) 

#-------------------------------------------------------------------------------
# 5. FIT LME EXACT MATCH MODEL

# Remove emotionBinary column from dataframe for Exact Match analysis
dfTrial_EM <- subset(dfTrial, select = -c(emotionBinary))
# Create trial pairs that are matched on all features using spread function
dfTrial_EM <- spread(dfTrial_EM, emotion, meanAmpNC)
# Calculate Full-minus-Reduced-Intensity Angry difference wave
dfTrial_EM$diffWaveMeanAmpNC <- dfTrial_EM$FullAngry - dfTrial_EM$ReduAngry

# Specify Exact Match model formula where predictor of interest is ECBQ_per
formulaLME_EM <- diffWaveMeanAmpNC ~ presentNumber + ECBQ_per + (1|SUBJECTID)
LME_EM <- lmer(formulaLME_EM, data=dfTrial_EM, REML=TRUE)
summary(LME_EM) 

# Check LME assumptions:
# - Linearity (Need to remove NA rows from dfTrial_EM in order to compare with model residuals)
plot(resid(LME_EM),dfTrial_EM$diffWaveMeanAmpNC[!is.na(dfTrial_EM$diffWaveMeanAmpNC)]) 
# - Normal distribution of residuals
qqmath(LME_EM) 
# - Homoscedasticity of emotion condition not applicable because only 1 difference wave condition 

#-------------------------------------------------------------------------------
# 6. FIT LME NEAREST NEIGHBOR: NO PRIORITIZATION MODEL 
# This Nearest Neighbor model is abbreviated as NN_N

# Identify trial pairs have an exact match based on SUBJECTID and electrode and a 
# greedy match based on presentation number and actor 
matchOutput_NN_N <- matchit(emotionBinary ~ presentNumber + ACTOR, data = dfTrial,
                            method = "nearest", distance = "mahalanobis", exact = c("SUBJECTID", "electrode"))
# Extract data with paired trials only (i.e., exclude unmatched trials)
dfTrial_NN_N <- match.data(matchOutput_NN_N) 
# Order rows from trial pairing ID (i.e., subclass) 1-s; and within each subclass: emotion FullAngry, then ReduAngry 
# (important for makeWide_realData function which assumes that each trial pair has been sorted by condition) 
# - Note that we don't need to sort based on SUBJECTID and electrode because that has
# already been taken into account by the trial pairing ID.
dfTrial_NN_N <- dfTrial_NN_N[order(dfTrial_NN_N$subclass, dfTrial_NN_N$emotion),]
# Make dataframe wide so we can calculate mean amp difference in next line
dfTrial_NN_N <- makeWide_realData(dfTrial_NN_N) 

# Calculate mean amplitude difference between trial pairs
dfTrial_NN_N$diffWaveMeanAmpNC <- dfTrial_NN_N$meanAmpNC.FullAngry - dfTrial_NN_N$meanAmpNC.ReduAngry

# Calculate average presentation number for each trial pair
dfTrial_NN_N$presentNumberAvg <- (dfTrial_NN_N$presentNumber.FullAngry + dfTrial_NN_N$presentNumber.ReduAngry)/2

# Specify Nearest Neighbor model formula where predictor of interest is ECBQ_per
formulaLME_NN <- diffWaveMeanAmpNC ~ presentNumberAvg + ECBQ_per + (1|SUBJECTID)
LME_NN_N <- lmer(formulaLME_NN, data=dfTrial_NN_N, REML=TRUE)
summary(LME_NN_N) 

# Check LME assumptions:
# - Linearity
plot(resid(LME_NN_N),dfTrial_NN_N$diffWaveMeanAmpNC)
# - Normal distribution of residuals
qqmath(LME_NN_N) 
# - Homoscedasticity of emotion condition not applicable because only 1 difference wave condition 

#-------------------------------------------------------------------------------
# 7. FIT LME NEAREST NEIGHBOR: PRESENTATION NUMBER PRIORITIZATION MODEL
# This Nearest Neighbor model is abbreviated as NN_P

# Identify trial pairs that have an exact match based on SUBJECTID, electrode, and
# presentation number and a greedy match based on actor 
matchOutput_NN_P <- matchit(emotionBinary ~ ACTOR, data = dfTrial,
                            method = "nearest", distance = "mahalanobis", exact = c("SUBJECTID", "electrode", "presentNumber"))
# Repeat same process of extracting trial-paired data and calculating difference wave 
# mean amplitude (see Section 6 for more detailed comments)
dfTrial_NN_P <- match.data(matchOutput_NN_P) 
dfTrial_NN_P <- dfTrial_NN_P[order(dfTrial_NN_P$subclass, dfTrial_NN_P$emotion),]
dfTrial_NN_P <- makeWide_realData(dfTrial_NN_P) 
dfTrial_NN_P$diffWaveMeanAmpNC <- dfTrial_NN_P$meanAmpNC.FullAngry - dfTrial_NN_P$meanAmpNC.ReduAngry
dfTrial_NN_P$presentNumberAvg <- (dfTrial_NN_P$presentNumber.FullAngry + dfTrial_NN_P$presentNumber.ReduAngry)/2
# Fit model using same Nearest Neighbor LME equation as in Section 6
LME_NN_P <- lmer(formulaLME_NN, data=dfTrial_NN_P, REML=TRUE)
summary(LME_NN_P) 

# Check LME assumptions:
# - Linearity
plot(resid(LME_NN_P),dfTrial_NN_P$diffWaveMeanAmpNC) 
# - Normal distribution of residuals
qqmath(LME_NN_P)
# - Homoscedasticity of emotion condition not applicable because only 1 difference wave condition 

#-------------------------------------------------------------------------------
# 8. FIT LME NEAREST NEIGHBOR: STIMULUS FEATURE PRIORITIZATION MODEL
# This Nearest Neighbor model is abbreviated as NN_F

# Identify trial pairs that have an exact match based on SUBJECTID, electrode, and
# actor (stimulus feature) and a greedy match based on presentation number
matchOutput_NN_F <- matchit(emotionBinary ~ presentNumber, data = dfTrial,
                            method = "nearest", distance = "mahalanobis", exact = c("SUBJECTID", "electrode", "ACTOR"))
# Repeat same process of extracting trial-paired data and calculating difference wave 
# mean amplitude (see Section 6 for more detailed comments)
dfTrial_NN_F <- match.data(matchOutput_NN_F) 
dfTrial_NN_F <- dfTrial_NN_F[order(dfTrial_NN_F$subclass, dfTrial_NN_F$emotion),]
dfTrial_NN_F <- makeWide_realData(dfTrial_NN_F) 
dfTrial_NN_F$diffWaveMeanAmpNC <- dfTrial_NN_F$meanAmpNC.FullAngry - dfTrial_NN_F$meanAmpNC.ReduAngry
dfTrial_NN_F$presentNumberAvg <- (dfTrial_NN_F$presentNumber.FullAngry + dfTrial_NN_F$presentNumber.ReduAngry)/2
# Fit model using same Nearest Neighbor LME equation as in Section 6
LME_NN_F <- lmer(formulaLME_NN, data=dfTrial_NN_F, REML=TRUE)
summary(LME_NN_F) 

# Check LME assumptions:
# - Linearity
plot(resid(LME_NN_F),dfTrial_NN_F$diffWaveMeanAmpNC) 
# - Normal distribution of residuals
qqmath(LME_NN_F) 
# - Homoscedasticity of emotion condition not applicable because only 1 difference wave condition 

#-------------------------------------------------------------------------------
# 9. FIT LME RANDOM PERMUTATION MODEL
# Note: We randomly pair trials WITHIN SUBJECTID and electrode

# Set up parallel processing for Random Permutation model
numCore <- 7
clust <- makeCluster(numCore) 
registerDoParallel(clust)
showConnections()

# Specify number of random permutation iterations
iterN <- 10000
origFuns <- setNSFunctions()

# Prepare dataframe for more efficient permutation by extracting number of trials per condition
# for each participant
dfTrial_order <- dfTrial[order(dfTrial$SUBJECTID, dfTrial$electrode, dfTrial$emotion), ] # Order rows from participant 1-n; and within each participant: order trials first by electrode, then by condition 
subjectEmotionNTable <- table(dfTrial$SUBJECTID, dfTrial$emotion) # Count number of trials for each participant and condition (note this is a multiple of 3 due to 3 electrodes)
subjectEmotionNTable <- subjectEmotionNTable[,c(1:2)] # Subset first two columns for subsequent processing
LME_output_RP_allIter <- vector("list", iterN)  # Temporary variable for storing output for each random perm iteration

# Specify Random Permutation model formula where predictor of interest is ECBQ_per 
# (same as Nearest Neighbor model formula)
formulaLME_RP <- diffWaveMeanAmpNC ~ presentNumberAvg + ECBQ_per + (1|SUBJECTID)

set.seed(20230317) # Set seed for reproducible results
LME_output_RP_allIter <- foreach (i=1:iterN, .packages=c('tidyverse', 'performance', 'emmeans', 'lme4', 'lmerTest', 'afex')) %dorng% { # Loop through each simulated data sample
  # For each row of subjectEmotionNTable (corresponding to the number of trials per condition for each subject):
  # 1) divide number of trials by 3 to get number of trials for one electrode
  # 2) randomly order integers 1-trialN_FullAngry; then randomly order integers 1-trialN_ReduAngry (no replacement);
  # 3) concatenate these two arrays into one array and then use unlist to save in one array;
  # 4) add this array as the subclass column in dfTrial_order
  # Note: The subclass pairing assumes that dfTrial_order has already been ordered based on participant, electrode, and condition! 
  dfTrial_order$subclass <- unlist(apply(subjectEmotionNTable, 1, 
                                            function(x) c(sample.int(x[1]/3), sample.int(x[2]/3), sample.int(x[1]/3), sample.int(x[2]/3), sample.int(x[1]/3), sample.int(x[2]/3))))
  
  # Order rows from participant 1-n; and within each participant: subclass 1-s; and within each subclass: condition FullAngry/ReduAngry
  # (important because later we assume that each first row is FullAngry, each second row is ReduAngry)
  dfTrial_order <- dfTrial_order[order(dfTrial_order$SUBJECTID, dfTrial_order$electrode, dfTrial_order$subclass), ]
  
  # For each row of subjectEmotionNTable:
  # 1) identify the min number of trials across conditions 
  # 2) select the first min number of trial pairs by assigning TRUE to the first 2*min number of trial
  #    and assign FALSE to all remaining rows (i.e., unpaired rows)
  keepTrialPairs <- unlist(apply(subjectEmotionNTable, 1, 
                                function(x){subjectEmotionN_min = min(x)/3; rep(rep(c(TRUE, FALSE), c(2*subjectEmotionN_min, (sum(x)/3) - 2*subjectEmotionN_min)), times = 3)}))
  
  # Select the rows that have a trial pair and discard all rows that do not
  dfTrial_RP <- dfTrial_order[keepTrialPairs, ]
  
  # Make dataframe wide so we can calculate mean amp difference in next line
  dfTrial_RP <- makeWide_realData(dfTrial_RP) 
  
  # Calculate mean amp difference between trial pairs
  dfTrial_RP$diffWaveMeanAmpNC <- dfTrial_RP$meanAmpNC.FullAngry - dfTrial_RP$meanAmpNC.ReduAngry
  
  # Calculate average presentation number for each trial pair
  dfTrial_RP$presentNumberAvg <- (dfTrial_RP$presentNumber.FullAngry + dfTrial_RP$presentNumber.ReduAngry)/2
  
  # Fit LME model
  LME_RP <- lmer(formulaLME_RP, data=dfTrial_RP, REML=TRUE)

  # If the model does not have problems
  if (check_convergence(LME_RP) && !check_singularity(LME_RP)) {
    # Extract emmeans for ECBQ perceptual predictor
    mLME_RP <- emtrends(LME_RP, ~ECBQ_per, var = "ECBQ_per", 
                        at = list(ECBQ_per = 4, presentNumber = 5.5),
                        mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak',
                        infer = c(TRUE, TRUE))
    
    # Effect is from emtrends, coef is from summary() 
    LME_output_RP <- data.frame(ECBQ_perEffect = unname(summary(mLME_RP))[2],
                                ECBQ_perCoef = unname(fixef(LME_RP)[3]),
                                presentNumCoef = unname(fixef(LME_RP)[2]), 
                                intercept = unname(fixef(LME_RP)[1]),
                                cor_IntPN = cov2cor(vcov(LME_RP))[2,1],
                                cor_IntECBQ = cov2cor(vcov(LME_RP))[3,1],
                                cor_PNECBQ = cov2cor(vcov(LME_RP))[3,2],
                                subjectIntercept_sd = data.frame(VarCorr(LME_RP))[1,5], # this is standard dev!
                                subjectIntercept_var = data.frame(VarCorr(LME_RP))[1,4],
                                resid_sd = data.frame(VarCorr(LME_RP))[2,5], 
                                resid_var = data.frame(VarCorr(LME_RP))[2,4],
                                trialN = 2*nobs(LME_RP),
                                subjectN = nrow(ranef(LME_RP)$SUBJECTID),
                                notConvergeN = 0, singFitN = 0, iterN = i)
  }
  # Check for model problems and update this iteration's output accordingly
  else if (!check_convergence(LME_RP)) {
    LME_output_RP <- data.frame(ECBQ_perEffect = NA,
                                ECBQ_perCoef = NA,
                                presentNumCoef = NA,
                                intercept = NA,
                                cor_IntPN = NA,
                                cor_IntECBQ = NA,
                                cor_PNECBQ = NA,
                                subjectIntercept_sd = data.frame(VarCorr(LME_RP))[1,5], 
                                subjectIntercept_var = data.frame(VarCorr(LME_RP))[1,4],
                                resid_sd = data.frame(VarCorr(LME_RP))[2,5], 
                                resid_var = data.frame(VarCorr(LME_RP))[2,4],
                                trialN = 2*nobs(LME_RP),
                                subjectN = nrow(ranef(LME_RP)$SUBJECTID),
                                notConvergeN = 1, singFitN = 0, iterN = i)
  }
  else if (check_singularity(LME_RP)) {
    LME_output_RP <- data.frame(ECBQ_perEffect = NA,
                                ECBQ_perCoef = NA,
                                presentNumCoef = NA,
                                intercept = NA,
                                cor_IntPN = NA,
                                cor_IntECBQ = NA,
                                cor_PNECBQ = NA,
                                subjectIntercept_sd = data.frame(VarCorr(LME_RP))[1,5], 
                                subjectIntercept_var = data.frame(VarCorr(LME_RP))[1,4],
                                resid_sd = data.frame(VarCorr(LME_RP))[2,5], 
                                resid_var = data.frame(VarCorr(LME_RP))[2,4],
                                trialN = 2*nobs(LME_RP),
                                subjectN = nrow(ranef(LME_RP)$SUBJECTID),
                                notConvergeN = 0, singFitN = 1, iterN = i)
  }
  LME_output_RP
}

# Combine results from list to dataframe 
LME_output_RP_allIterFinal <- do.call(rbind, LME_output_RP_allIter)

# Remove rows where predictor of interest was NA (i.e., model had problems) 
LME_output_RP_allIterFinal <- LME_output_RP_allIterFinal[!is.na(LME_output_RP_allIterFinal$ECBQ_perEffect),]

# Calculate mean and non-parametric bias-corrected and accelerated bootstrap 
# confidence intervals for ECBQ perceptual sensitivity effect
mean(LME_output_RP_allIterFinal$ECBQ_perEffect) # Examine mean of perceptual sensitivity 
bca(LME_output_RP_allIterFinal$ECBQ_perEffect) # Examine if confidence intervals include 0

# Reset functions modified for speeding up random permutation and close cluster
resetNSFunctions(origFuns)
stopCluster(clust)
showConnections()

#-------------------------------------------------------------------------------
# 10. FIT 10-TRIAL REGRESSION MODEL

# Calculate mean-averaged difference waves from condition-level data
dfMeanAvg_10TrialMin <- spread(dfCond_10TrialMin, emotion, meanAmpNC)
dfMeanAvg_10TrialMin$diffWaveMeanAmpNC <- dfMeanAvg_10TrialMin$FullAngry - dfMeanAvg_10TrialMin$ReduAngry

# Specify regression model formula where predictor of interest is ECBQ_per
# (same for both 10-Trial and 15-Trial Regression) 
formulaReg = diffWaveMeanAmpNC ~ ECBQ_per
reg_10TrialMin <- lm(formulaReg, data = dfMeanAvg_10TrialMin)
summary(reg_10TrialMin)

# Check regression assumptions:
# - Linearity
plot(resid(reg_10TrialMin),dfMeanAvg_10TrialMin$diffWaveMeanAmpNC) 
# - Normal distribution of residuals
# We use the shapiro_test function, which requires sample size between 3 and 5000 
# (not applicable for trial-level data)
shapiro_test(residuals(reg_10TrialMin)) 
# - Homoscedasticity of emotion condition not applicable because only 1 difference wave condition 

#-------------------------------------------------------------------------------
# 11. FIT 15-TRIAL REGRESSION MODEL

# Calculate mean-averaged difference waves from condition-level data
dfMeanAvg_15TrialMin <- spread(dfCond_15TrialMin, emotion, meanAmpNC)
dfMeanAvg_15TrialMin$diffWaveMeanAmpNC <- dfMeanAvg_15TrialMin$FullAngry - dfMeanAvg_15TrialMin$ReduAngry

# Fit same regression model formula as for 10-Trial Regression (see Section 10) 
reg_15TrialMin <- lm(formulaReg, data = dfMeanAvg_15TrialMin)
summary(reg_15TrialMin)

# Check regression assumptions:
# - Linearity
plot(resid(reg_15TrialMin),dfMeanAvg_15TrialMin$diffWaveMeanAmpNC) 
# - Normal distribution of residuals
shapiro_test(residuals(reg_15TrialMin))
# - Homoscedasticity of emotion condition not applicable because only 1 difference wave condition 

#-------------------------------------------------------------------------------
# 12. COMPARE ECBQ PERCEPTUAL SENSITIVITY ESTIMATES ACROSS MODELS 

# Use emtrends for all models so that confidence intervals are comparable 
# - We pick ECBQ_per of 4 because this is the median
# - We pick presentNumber of 5.5. because this is the median
mLME_Int <- emtrends(LME_Int, pairwise ~ emotion, var = "ECBQ_per", 
                     at = list(ECBQ_per = 4, presentNumber = 5.5), 
                     mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak',
                     infer=c(TRUE,TRUE))
mLME_EM <- emtrends(LME_EM, ~ECBQ_per, var = "ECBQ_per", 
                    at = list(ECBQ_per = 4, presentNumber = 5.5),
                    mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak',
                    infer = c(TRUE, TRUE))
mLME_RP <- data.frame(estimate = mean(LME_output_RP_allIterFinal$ECBQ_perEffect),
                      lower.CL = bca(LME_output_RP_allIterFinal$ECBQ_perEffect)[1],
                      upper.CL = bca(LME_output_RP_allIterFinal$ECBQ_perEffect)[2])
mLME_NN_N <- emtrends(LME_NN_N, ~ECBQ_per, var = "ECBQ_per", 
                      at = list(ECBQ_per = 4, presentNumberAvg = 5.5),
                      mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak',
                      infer = c(TRUE, TRUE))
mLME_NN_P <- emtrends(LME_NN_P, ~ECBQ_per, var = "ECBQ_per", 
                       at = list(ECBQ_per = 4, presentNumberAvg = 5.5),
                       mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak',
                       infer = c(TRUE, TRUE))
mLME_NN_F <- emtrends(LME_NN_F, ~ECBQ_per, var = "ECBQ_per", 
                       at = list(ECBQ_per = 4, presentNumberAvg = 5.5),
                       mode = "satterthwaite", lmerTest.limit = 240000, adjust='sidak',
                       infer = c(TRUE, TRUE))
mReg_10TrialMin <- emtrends(reg_10TrialMin, ~ECBQ_per, var = "ECBQ_per",
                            at = list(ECBQ_per = 4),
                            adjust='sidak', infer = c(TRUE, TRUE))
mReg_15TrialMin <- emtrends(reg_15TrialMin, ~ECBQ_per, var = "ECBQ_per",
                            at = list(ECBQ_per = 4),
                            adjust='sidak', infer = c(TRUE, TRUE))