# Difference Wave Simulation Helper Functions
  # - induceMissingTrials_v2
  # - allTrials_Int
  # - pairTrials_EM
  # - pairTrials_NN_N
  # - pairTrials_NN_P
  # - pairTrials_NN_F
  # - pairTrials_RP
  # - makeWide
  # - fitANOVA_diffWave
  # - fitLME

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

#-------------------------------------------------------------------------------
# induceMissingTrials_v2: Function to induce missing trials according to whether the
# subject has been selected to have a low trial-count. This function was adapted 
# from the original induceMissingTrials function used in Heise et al. (2022).
# - Format: 
  #   dfMissing_noNA <- induceMissingTrials_v2(dfOriginal, subjectCaseDeletion) 
# - Inputs:
  # - dfOriginal: Dataframe with the simulated "population" dataset before any
  #   induced missing trials. 
  # - subjectCaseDeletion: Array of subject IDs that have been selected to have 
  #   low trial-counts (i.e., less than 10 trials/condition).
# - Output: 
  # - dfMissing_noNA: A copy of the dfOriginal after rows with induced missing 
  #   trials have been removed. 
induceMissingTrials_v2 <- function(dfOriginal, subjectCaseDeletion) {
  # Create copy of the dfOriginal dataframe for inducing missing trials
  dfMissing <- data.frame(dfOriginal)
  
  # Loop through each subject and randomly select a subset of trials from each
  # condition to remove
  for (subject in unique(dfOriginal$SUBJECTID)) {
    if (subject %in% subjectCaseDeletion) { # If the subject has been selected for a low trial-count
      # Generate how many trials are missing so that at least one emotion condition 
      # will have fewer than 10 trials remaining 
      trialMissing <- c(sample(x=(trialMissingThreshold+1):emotionTrialN, size = 1),
                        sample(x=0:emotionTrialN, size = emotionN-1, replace = TRUE))   
      
    } else { # If the subject has NOT been selected for a low trial-count
      # Generate how many missing trials so that all emotion conditions will have
      # at least 10 trials
      trialMissing <- sample(x=0:trialMissingThreshold, size = emotionN, 
                             replace = TRUE)
    }
    
    # Shuffle order of emotion conditions and then loop through each condition
    # (This line is added so that one condition does not consistently have more
    # missing trials due to how the trialMissing variable is defined for subjects
    # with low trial-count.)
    emotionLabelRand <- sample(emotionLabel)
    for (j in 1:length(emotionLabelRand)) {
      # Extract the number of missing trials for this condition
      emotionTrialMissing <- trialMissing[j] 
      
      if (emotionTrialMissing != 0) { # If this emotion condition was not selected to have 0 missing trials
        # Find this subject and condition's rows in the dataframe
        subjectEmotionIndex <- which(dfMissing$SUBJECTID==subject & dfMissing$emotion==emotionLabelRand[j])
        
        # Extract the trial presentation number probability weights for this subject
        # and condition
        subjectEmotionProbWeight <- dfMissing$presentNumberWeight[subjectEmotionIndex]
        
        # Randomly select the missing trials based on the specified probability 
        # weights and the number of missing trials
        subjectEmotionIndexMissing <- sample(subjectEmotionIndex, emotionTrialMissing,
                                             replace=FALSE, prob=subjectEmotionProbWeight)
        
        # For these missing trials only, replace the meanAmpNC value with NA
        dfMissing[subjectEmotionIndexMissing,]$meanAmpNC <- NA
      }
    }
  }
  
  # Return dataframe without NA rows and only necessary columns
  return(dfMissing_noNA = dfMissing[!is.na(dfMissing$meanAmpNC), c("SUBJECTID", "mSens", "age", "emotion", "ACTOR", "presentNumber", "meanAmpNC")]) 
}

#-------------------------------------------------------------------------------
# allTrials_Int: Fit an LME Interaction (Int) model to the input dataframe.
# - Format: 
  #   modelOutput <- allTrials_Int(dfInput_noNA, formulaInput) 
# - Inputs:
  # - dfInput_noNA: Long dataframe with each row corresponding to the mean amplitude
  #   for each emotion condition/actor/trial presentation number/subject. There are
  #   no rows with NA values.
  # - formulaInput: LME model formula (e.g., meanAmpNC ~ emotion + presentNumber + mSens + age + mSens:emotion + age:emotion + presentNumber:emotion + (1|SUBJECTID) + (1|ACTOR))
# - Output: 
  # - modelOutput: Dataframe containing the estimated marginal means for each maternal 
  #   sensitivity group and other model predictors.
allTrials_Int <- function(dfInput_noNA, formulaInput) {
  environment(formulaInput) <- environment() # Added based on https://stackoverflow.com/questions/60510732/lapply-fitting-thousands-of-mixed-models-and-being-able-to-extract-lsmeans
  
  # Convert categorical variables to factor before fitting model
  dfInput_noNA$ACTOR <- as.factor(dfInput_noNA$ACTOR)
  dfInput_noNA$SUBJECTID <- as.factor(dfInput_noNA$SUBJECTID)
  dfInput_noNA$emotion <- as.factor(dfInput_noNA$emotion)
  # Specify effects coding for the age factor 
  dfInput_noNA$age <- as.factor(dfInput_noNA$age)
  contrasts(dfInput_noNA$age) <- contr.Sum(levels(dfInput_noNA$age)) 
  # Specify effects coding for the maternal sensitivity factor
  dfInput_noNA$mSens <- as.factor(dfInput_noNA$mSens)
  contrasts(dfInput_noNA$mSens) <- contr.Sum(levels(dfInput_noNA$mSens)) 
  
  # Fit LME model
  fit.LME <- lmer(formulaInput, data = dfInput_noNA, REML = TRUE)
  
  tryCatch({
    # If the model does not have problems (i.e., non-convergence or singular fit)
    if (check_convergence(fit.LME) && !check_singularity(fit.LME)) { 
      # Calculate marginal means for each maternal sensitivity group. We use 
      # revpairwise so that we extract the condition B-A contrast (comparable to 
      # the B-minus-A mean amplitude difference).
      mLME_mSens <- emmeans(fit.LME, revpairwise ~ emotion|mSens, mode = "satterthwaite", 
                            lmerTest.limit = 240000, at = list(presentNumber = presentNumberRef))
      # Extract maternal sensitivity marginal means
      LME_output_mSens_EM <- data.frame(summary(mLME_mSens)$contrasts)[,1:4]
      
      # Calculate marginal means for each age group
      mLME_age <- emmeans(fit.LME, revpairwise ~ emotion|age, mode = "satterthwaite", 
                          lmerTest.limit = 240000, at = list(presentNumber = presentNumberRef))
      # Extract age marginal means
      LME_output_age_EM <- data.frame(summary(mLME_age)$contrasts)[,1:4]
      
      # Create dataframe with model predictors and other information 
      return(data.frame(lowMSensEMM = LME_output_mSens_EM[2,3], lowMSensSE = LME_output_mSens_EM[2,4],
                        highMSensEMM = LME_output_mSens_EM[1,3], highMSensSE = LME_output_mSens_EM[1,4],
                        mSensEffect = unname(fixef(fit.LME)[6]), # mSens:emotion interaction estimate 
                        mSensPValue = summary(fit.LME)$coefficients[6,5],
                        youngEMM = LME_output_age_EM[2,3], youngSE = LME_output_age_EM[2,4],
                        oldEMM = LME_output_age_EM[1,3], oldSE = LME_output_age_EM[1,4],
                        ageEffect = unname(fixef(fit.LME)[7]), # age:emotion interaction estimate
                        agePValue = summary(fit.LME)$coefficients[7,5],
                        presentNumEffect = unname(fixef(fit.LME)[8]), # presentNumber:emotion interaction estimate
                        presentNumPValue = summary(fit.LME)$coefficients[8,5],
                        subjectIntercept = data.frame(VarCorr(fit.LME))[1,5],
                        actorIntercept = data.frame(VarCorr(fit.LME))[2,5],
                        trialN = nobs(fit.LME),
                        subjectN = nrow(ranef(fit.LME)$SUBJECTID),
                        notConvergeN = 0, singFitN = 0, modelType = "LME"))
    }
    # Check for model problems and return model output accordingly
    else if (!check_convergence(fit.LME)) { # 
      return(data.frame(lowMSensEMM = NA, lowMSensSE = NA,
                        highMSensEMM = NA, highMSensSE = NA,
                        mSensEffect = NA,
                        mSensPValue = NA,
                        youngEMM = NA, youngSE = NA,
                        oldEMM = NA, oldSE = NA,
                        ageEffect = NA,
                        agePValue = NA,
                        presentNumEffect = NA,
                        presentNumPValue = NA,
                        subjectIntercept = data.frame(VarCorr(fit.LME))[1,5],
                        actorIntercept = data.frame(VarCorr(fit.LME))[2,5],
                        trialN = nobs(fit.LME),
                        subjectN = nrow(ranef(fit.LME)$SUBJECTID),
                        notConvergeN = 1, singFitN = 0, modelType = "LME"))
    }
    else if (check_singularity(fit.LME)) {
      return(data.frame(lowMSensEMM = NA, lowMSensSE = NA,
                        highMSensEMM = NA, highMSensSE = NA,
                        mSensEffect = NA,
                        mSensPValue = NA,
                        youngEMM = NA, youngSE = NA,
                        oldEMM = NA, oldSE = NA,
                        ageEffect = NA,
                        agePValue = NA,
                        presentNumEffect = NA,
                        presentNumPValue = NA,
                        subjectIntercept = data.frame(VarCorr(fit.LME))[1,5],
                        actorIntercept = data.frame(VarCorr(fit.LME))[2,5],
                        trialN = nobs(fit.LME),
                        subjectN = nrow(ranef(fit.LME)$SUBJECTID),
                        notConvergeN = 0, singFitN = 1, modelType = "LME"))
    }
  }, 
  error=function(e) {
    if (e$message == "Lapack routine dgesv: system is exactly singular: U[1,1] = 0") {
      return(data.frame(lowMSensEMM = NA, lowMSensSE = NA,
                        highMSensEMM = NA, highMSensSE = NA,
                        mSensEffect = NA,
                        mSensPValue = NA,
                        youngEMM = NA, youngSE = NA,
                        oldEMM = NA, oldSE = NA,
                        ageEffect = NA,
                        agePValue = NA,
                        presentNumEffect = NA,
                        presentNumPValue = NA,
                        subjectIntercept = data.frame(VarCorr(fit.LME))[1,5],
                        actorIntercept = data.frame(VarCorr(fit.LME))[2,5], 
                        trialN = nobs(fit.LME),
                        subjectN = nrow(ranef(fit.LME)$SUBJECTID),
                        notConvergeN = 0, singFitN = 1, modelType = "LME")) 
    }
    else{
      return(data.frame(lowMSensEMM = NA, lowMSensSE = NA,
                        highMSensEMM = NA, highMSensSE = NA,
                        mSensEffect = NA,
                        mSensPValue = NA,
                        youngEMM = NA, youngSE = NA,
                        oldEMM = NA, oldSE = NA,
                        ageEffect = NA,
                        agePValue = NA,
                        presentNumEffect = NA,
                        presentNumPValue = NA,
                        subjectIntercept = data.frame(VarCorr(fit.LME))[1,5],
                        actorIntercept = data.frame(VarCorr(fit.LME))[2,5], 
                        trialN = nobs(fit.LME),
                        subjectN = nrow(ranef(fit.LME)$SUBJECTID),
                        notConvergeN = -99, singFitN = -99, modelType = "LME")) # -99 designates other type of error
    }
  }) 
}

#-------------------------------------------------------------------------------
# pairTrials_EM: Fit an LME Exact Match (EM) model to the input dataframe.
# - Format: 
  #   modelOutput <- pairTrials_EM(dfInput_noNA, formulaInput) 
# - Inputs:
  # - dfInput_noNA: Long dataframe with each row corresponding to the mean amplitude
  #   for each emotion condition/actor/trial presentation number/subject. There are
  #   no rows with NA values.
  # - formulaInput: LME model formula (e.g., meanAmpNC_BMinusA ~ presentNumberAvg + mSens + age + (1|SUBJECTID))
# - Output: 
  # - modelOutput: Dataframe containing the estimated marginal means for each maternal 
  #   sensitivity group and other model predictors.
pairTrials_EM <- function(dfInput_noNA, formulaInput) {
  # Create trial pairs that are matched on all features using reshape function
  dfInput_pairedWide <- stats::reshape(dfInput_noNA, 
                                       v.names = "meanAmpNC",
                                       idvar = c("SUBJECTID", "ACTOR", "presentNumber"),
                                       timevar = "emotion",
                                       direction = "wide")
  
  # Calculate B-minus-A difference wave
  dfInput_pairedWide$meanAmpNC_BMinusA <- dfInput_pairedWide$meanAmpNC.B - dfInput_pairedWide$meanAmpNC.A
  
  # Create presentNumberAvg variable so that we can use the same formulaInput 
  # as the one used for other LME difference wave approaches
  dfInput_pairedWide$presentNumberAvg <- dfInput_pairedWide$presentNumber
  
  # Fit LME model and return model output
  return(fitLME(dfInput_pairedWide, formulaInput))
}

#-------------------------------------------------------------------------------
# pairTrials_NN_N: Fit an LME Nearest Neighbor: No Prioritization (NN:N) model to
# the input dataframe.
# - Format: 
  #   modelOutput <- pairTrials_NN_N(dfInput_noNA, formulaInput) 
# - Inputs:
  # - dfInput_noNA: Long dataframe with each row corresponding to the mean amplitude
  #   for each emotion condition/actor/trial presentation number/subject. There are
  #   no rows with NA values.
  # - formulaInput: LME model formula (e.g., meanAmpNC_BMinusA ~ presentNumberAvg + mSens + age + (1|SUBJECTID))
# - Output: 
  # - modelOutput: Dataframe containing the estimated marginal means for each maternal 
  #   sensitivity group and other model predictors.
pairTrials_NN_N <- function(dfInput_noNA, formulaInput){
  # Create column of 0's and 1's where condition A is set as the 'treatment' group
  # (1) for matchit function
  dfInput_noNA$emotionBinary <- ifelse(dfInput_noNA$emotion == "A", 1, 0)

  # Identify trial pairs that have an exact match based on SUBJECTID and a 
  # greedy match based on trial presentation number and actor 
  matchOutput <- matchit(emotionBinary ~ presentNumber + ACTOR, data = dfInput_noNA,
                         method = "nearest", distance = "mahalanobis", exact = "SUBJECTID")
    
  # Extract data with paired trials only (i.e., exclude unmatched trials)
  dfInput_noNAMatch <- match.data(matchOutput) 
  
  # Order rows from trial pairing ID (i.e., subclass) 1-s; and within each subclass: condition A, then condition B
  # (important for makeWide function which assumes that each trial pair has been sorted by condition) 
  dfInput_noNAMatch <- dfInput_noNAMatch[order(dfInput_noNAMatch$subclass, dfInput_noNAMatch$emotion), ]
  
  # Select the columns of interest
  dfInput_noNAMatch <- dfInput_noNAMatch[, c("SUBJECTID", "mSens", "age", "presentNumber", "meanAmpNC")]
  
  # Make dataframe wide so we can calculate trial-level difference waves
  dfInput_pairedWide <- makeWide(dfInput_noNAMatch) 
  
  # Calculate trial-level difference wave mean amplitudes
  dfInput_pairedWide$meanAmpNC_BMinusA <- dfInput_pairedWide$meanAmpNC.B - dfInput_pairedWide$meanAmpNC.A
  
  # Calculate average trial presentation number for each trial pair
  dfInput_pairedWide$presentNumberAvg <- (dfInput_pairedWide$presentNumber.A + dfInput_pairedWide$presentNumber.B)/2
  
  # Fit LME model and return model output 
  return(fitLME(dfInput_pairedWide, formulaInput))
}

#-------------------------------------------------------------------------------
# pairTrials_NN_P: Fit an LME Nearest Neighbor: Presentation Number Prioritization (NN:P)
# model to the input dataframe.
# - Format: 
  #   modelOutput <- pairTrials_NN_P(dfInput_noNA, formulaInput) 
# - Inputs:
  # - dfInput_noNA: Long dataframe with each row corresponding to the mean amplitude
  #   for each emotion condition/actor/trial presentation number/subject. There are
  #   no rows with NA values.
  # - formulaInput: LME model formula (e.g., meanAmpNC_BMinusA ~ presentNumberAvg + mSens + age + (1|SUBJECTID))
# - Output: 
  # - modelOutput: Dataframe containing the estimated marginal means for each maternal 
  #   sensitivity group and other model predictors.
pairTrials_NN_P <- function(dfInput_noNA, formulaInput){
  # Create column of 0's and 1's where condition A is set as the 'treatment' group
  # (1) for matchit function 
  dfInput_noNA$emotionBinary <- ifelse(dfInput_noNA$emotion == "A", 1, 0)
  
  # Identify trial pairs that have an exact match based on SUBJECTID and trial
  # presentation number and a greedy match based on actor 
  matchOutput <- matchit(emotionBinary ~ ACTOR, data = dfInput_noNA,
                         method = "nearest", distance = "mahalanobis", exact = c("SUBJECTID", "presentNumber"))
  
  # Extract trial-paired data and calculate difference wave mean amplitude using 
  # same process as pairTrials_NN_N (see function above for more detailed comments) 
  dfInput_noNAMatch <- match.data(matchOutput) 
  dfInput_noNAMatch <- dfInput_noNAMatch[order(dfInput_noNAMatch$subclass, dfInput_noNAMatch$emotion), ]
  dfInput_noNAMatch <- dfInput_noNAMatch[, c("SUBJECTID", "mSens", "age", "presentNumber", "meanAmpNC")]
  dfInput_pairedWide <- makeWide(dfInput_noNAMatch)
  dfInput_pairedWide$meanAmpNC_BMinusA <- dfInput_pairedWide$meanAmpNC.B - dfInput_pairedWide$meanAmpNC.A
  dfInput_pairedWide$presentNumberAvg <- (dfInput_pairedWide$presentNumber.A + dfInput_pairedWide$presentNumber.B)/2
  
  return(fitLME(dfInput_pairedWide, formulaInput))
}

#-------------------------------------------------------------------------------
# pairTrials_NN_F: Fit an LME Nearest Neighbor: Stimulus Feature Prioritization (NN:F)
# model to the input dataframe.
# - Format: 
  #   modelOutput <- pairTrials_NN_F(dfInput_noNA, formulaInput) 
# - Inputs:
  # - dfInput_noNA: Long dataframe with each row corresponding to the mean amplitude
  #   for each emotion condition/actor/trial presentation number/subject. There are
  #   no rows with NA values.
  # - formulaInput: LME model formula (e.g., meanAmpNC_BMinusA ~ presentNumberAvg + mSens + age + (1|SUBJECTID))
# - Output: 
  # - modelOutput: Dataframe containing the estimated marginal means for each maternal 
  #   sensitivity group and other model predictors.
pairTrials_NN_F <- function(dfInput_noNA, formulaInput){
  # Create column of 0's and 1's where condition A is set as the 'treatment' group
  # (1) for matchit function 
  dfInput_noNA$emotionBinary <- ifelse(dfInput_noNA$emotion == "A", 1, 0)
  
  # Identify trial pairs that have an exact match based on SUBJECTID and actor
  # (stimulus feature) and a greedy match based on trial presentation number
  matchOutput <- matchit(emotionBinary ~ presentNumber, data = dfInput_noNA,
                         method = "nearest", distance = "mahalanobis", exact = c("SUBJECTID", "ACTOR"))
  
  # Extract trial-paired data and calculate difference wave mean amplitude using 
  # same process as pairTrials_NN_N (see function above for more detailed comments) 
  dfInput_noNAMatch <- match.data(matchOutput) 
  dfInput_noNAMatch <- dfInput_noNAMatch[order(dfInput_noNAMatch$subclass, dfInput_noNAMatch$emotion), ]
  dfInput_noNAMatch <- dfInput_noNAMatch[, c("SUBJECTID", "mSens", "age", "presentNumber", "meanAmpNC")]
  dfInput_pairedWide <- makeWide(dfInput_noNAMatch)
  dfInput_pairedWide$meanAmpNC_BMinusA <- dfInput_pairedWide$meanAmpNC.B - dfInput_pairedWide$meanAmpNC.A
  dfInput_pairedWide$presentNumberAvg <- (dfInput_pairedWide$presentNumber.A + dfInput_pairedWide$presentNumber.B)/2
  
  return(fitLME(dfInput_pairedWide, formulaInput))
}

#-------------------------------------------------------------------------------
# pairTrials_RP: Fit an LME Random Permutation (RP) model to the input dataframe.
# - Format: 
  #   modelOutput <- pairTrials_RP(dfInput_noNA, formulaInput, rpIter) 
# - Inputs:
  # - dfInput_noNA: Long dataframe with each row corresponding to the mean amplitude
  #   for each emotion condition/actor/trial presentation number/subject. There are
  #   no rows with NA values.
  # - formulaInput: LME model formula (e.g., meanAmpNC_BMinusA ~ presentNumberAvg + mSens + age + (1|SUBJECTID))
  # - rpIter: Number of Random Permutation iterations.
# - Output: 
  # - modelOutput: Dataframe containing the estimated marginal means for each maternal 
  #   sensitivity group and other model predictors for each iteration.
pairTrials_RP <- function(dfInput_noNA, formulaInput, rpIter) {
  dfInput_noNAOrder <- dfInput_noNA[order(dfInput_noNA$SUBJECTID, dfInput_noNA$emotion),] # Order rows from subject 1-n; and within each subject: order trials by emotion condition
  subjectEmotionNTable <- table(dfInput_noNA$SUBJECTID, dfInput_noNA$emotion) # Count number of trials for each subject and emotion 
  LME_output_RP_allIter <- vector("list", rpIter)  # Temporary variable for storing output for each random perm iteration
  
  for(i in 1:rpIter) { # Loop through each random perm iteration
    # For each row of subjectEmotionNTable (corresponding to the number of trials per condition for each subject):
    # 1) randomly order integers 1-trialN_A; then randomly order integers 1-trialN_B (no replacement)
    # 2) concatenate these two arrays into one array and then use unlist to save in one array
    # 3) add this array as the subclass column in dfInput_noNAOrder
    # Note: The subclass pairing assumes that dfInput_noNAOrder has already been ordered based on subject and condition! 
    dfInput_noNAOrder$subclass <-unlist(apply(subjectEmotionNTable, 1, 
                                              function(x) c(sample.int(x[1]), sample.int(x[2]))))
    
    # Order rows from subject 1-n; and within each subject: subclass 1-s; and within each subclass: condition A, B
    # (this is important for subsequent processing where we assume that each first row is A, each second row is B)
    dfInput_noNAFinal <- dfInput_noNAOrder[order(dfInput_noNAOrder$SUBJECTID, dfInput_noNAOrder$subclass, dfInput_noNAOrder$emotion), ]
    
    # For each row of subjectEmotionNTable:
    # 1) identify the min number of trials across conditions 
    # 2) select the first min number of trial pairs by assigning TRUE to the first 
    #    2*min number of trial and assign FALSE to all remaining rows (i.e., unpaired rows)
    keepTrialPairs <- unlist(apply(subjectEmotionNTable, 1, 
                                   function(x){subjectEmotionN_min = min(x); rep(c(TRUE, FALSE), c(2*subjectEmotionN_min, sum(x) - 2*subjectEmotionN_min))}))
    
    # Select the rows that have a trial pair and discard all rows that do not;  
    # also select columns of interest
    dfInput_noNAFinal <- dfInput_noNAFinal[keepTrialPairs, c("SUBJECTID", "mSens", "age", "presentNumber", "meanAmpNC")]
    
    # Make dataframe wide
    dfInput_pairedWide <- makeWide(dfInput_noNAFinal) 
    
    # Calculate trial-level difference wave mean amplitudes
    dfInput_pairedWide$meanAmpNC_BMinusA <- dfInput_pairedWide$meanAmpNC.B - dfInput_pairedWide$meanAmpNC.A
    
    # Calculate average trial presentation number for each trial pair
    dfInput_pairedWide$presentNumberAvg <- (dfInput_pairedWide$presentNumber.A + dfInput_pairedWide$presentNumber.B)/2
     
    # Fit LME model and record model output 
    LME_output_RP_allIter[[i]] <- fitLME(dfInput_pairedWide, formulaInput)
    LME_output_RP_allIter[[i]]$sampleRP <- i # Record random iteration 
  }
  # Combine results from list to dataframe 
  LME_output_RP_allIterFinal <- do.call(rbind, LME_output_RP_allIter)
  
  # Return model output for all iterations
  return(LME_output_RP_allIterFinal)
}

#-------------------------------------------------------------------------------
# makeWide: Function to convert long dataframe into wide to facilitate trial-level
# difference wave calculations. 
# - Format:
  #   dfOutput_wide <- makeWide(dfInput_long)
# - Input: 
  # - dfInput_long: Trial-level dataframe in which rows have been sorted based on trial pairing ID
  #   and condition (e.g., row 1 is trial pair 1/condition A, row 2 is trial pair 1/condition B,
  #   row 3 is trial pair 2/condition A, row 4 is trial pair 2/condition B, etc.). 
# - Output: 
  # - dfOutput_wide: Trial-paired dataframe in wide format in which each condition
  #   has columns for mean amplitude and trial presentation number.
makeWide = function(dfInput_long) {
  # Create array of indices for each pair of trials (extract 1st row of each pair)
  i <- seq(1, nrow(dfInput_long) - 1, by = 2)
  
  # For each 2nd row of each trial pair (corresponding to the second condition),
  # we extract mean amplitude and trial presentation number information and save 
  # as new columns in the 1st row of each pair 
  dfInput_long[i, c("meanAmpNC.B", "presentNumber.B")] = dfInput_long[i+1, c("meanAmpNC", "presentNumber")]
  
  # Rename the fourth and fifth column to reflect the first condition (note this 
  # is hard coded!)
  names(dfInput_long)[c(4,5)] = c("presentNumber.A","meanAmpNC.A")
  
  # Return the first row of each trial pair, which now contains the mean amplitude
  # and trial presentation number information for both conditions
  return(dfOutput_wide = dfInput_long[i, ])
}

#-------------------------------------------------------------------------------
# fitANOVA: Fit an ANOVA model to the input dataframe.
# - Format: 
  #   modelOutput <- fitANOVA(dfInput_noNA) 
# - Inputs:
  # - dfInput_noNA: Long dataframe with each row corresponding to the mean amplitude
  #   for each emotion condition/actor/trial presentation number/subject. There are
  #   no rows with NA values. Subjects with low trial-counts have been casewise
  #   deleted.
# - Output: 
  # - modelOutput: Dataframe containing the estimated marginal means for each maternal 
  #   sensitivity group and other model predictors.
fitANOVA = function(dfInput_noNA) {
  # Calculate mean-averaged data for each emotion condition and subject
  dfInputAvg <- aggregate(meanAmpNC ~ SUBJECTID + emotion + age + mSens, dfInput_noNA,
                          mean, na.action = na.omit)
  
  # Convert data from long to wide
  dfInputAvg_pairedWide <- reshape(dfInputAvg, 
                                   v.names = "meanAmpNC",
                                   idvar = "SUBJECTID",
                                   timevar = "emotion",
                                   direction = "wide")
  
  # Calculate mean-averaged difference wave mean amplitudes
  dfInputAvg_pairedWide$meanAmpNC_BMinusA <- dfInputAvg_pairedWide$meanAmpNC.B - dfInputAvg_pairedWide$meanAmpNC.A
  
  # Convert categorical variables to factor before fitting model
  dfInputAvg_pairedWide$SUBJECTID <- as.factor(dfInputAvg_pairedWide$SUBJECTID)
  # Specify effects coding for the age factor 
  dfInputAvg_pairedWide$age <- as.factor(dfInputAvg_pairedWide$age)
  contrasts(dfInputAvg_pairedWide$age) <- contr.Sum(levels(dfInputAvg_pairedWide$age)) 
  # Specify effects coding for the maternal sensitivity factor 
  dfInputAvg_pairedWide$mSens <- as.factor(dfInputAvg_pairedWide$mSens)
  contrasts(dfInputAvg_pairedWide$mSens) <- contr.Sum(levels(dfInputAvg_pairedWide$mSens)) 
  
  # Fit ANOVA model
  fit.ANOVA <- aov_ez(id = "SUBJECTID", dv = "meanAmpNC_BMinusA", data = dfInputAvg_pairedWide,
                      between = c("mSens", "age"))
  
  # Calculate and extract marginal means for each maternal sensitivity group 
  mANOVA_mSens <- data.frame(emmeans(fit.ANOVA, ~ mSens))[,1:3]
  
  # Calculate and extract marginal means for each age group 
  mANOVA_age <- data.frame(emmeans(fit.ANOVA, ~ age))[,1:3]
  
  # Create dataframe with model predictors and other information
  return(data.frame(lowMSensEMM = mANOVA_mSens[2,2], lowMSensSE = mANOVA_mSens[2,3],
                    highMSensEMM = mANOVA_mSens[1,2], highMSensSE = mANOVA_mSens[1,3],
                    mSensEffect = NA, 
                    mSensPValue = summary(fit.ANOVA)[[1,6]], 
                    youngEMM = mANOVA_age[2,2], youngSE = mANOVA_age[2,3], 
                    oldEMM = mANOVA_age[1,2], oldSE = mANOVA_age[1,3], 
                    ageEffect = NA, 
                    agePValue = summary(fit.ANOVA)[[2,6]],
                    presentNumEffect = NA, # Not relevant for ANOVA 
                    presentNumPValue = NA, # Not relevant for ANOVA
                    subjectIntercept = NA, # Not relevant for ANOVA
                    actorIntercept = NA, # Not relevant for ANOVA
                    trialN = nrow(dfInput_noNA), 
                    subjectN = nrow(dfInputAvg)/2, # Return number of subjects NOT casewise deleted
                    notConvergeN = 0, singFitN = 0, modelType = "ANOVA"))
}

#-------------------------------------------------------------------------------
# fitLME: Fit an LME model to the input dataframe with the specified formula. This
# function is used for fitting Exact Match, Nearest Neighbor, and Random Permutation
# models.
# - Format: 
  #   modelOutput <- fitLME(dfInput_pairedWide, formulaInput) 
# - Inputs:
  # - dfInput_pairedWide: Wide dataframe with each row corresponding to the trial-level
  #   difference wave mean amplitude. There are no rows with NA values.
  # - formulaInput: LME model formula (e.g., meanAmpNC_BMinusA ~ presentNumberAvg + mSens + age + (1|SUBJECTID))
# - Output: 
  # - modelOutput: Dataframe containing the estimated marginal means for each maternal 
  #   sensitivity group and other model predictors.
fitLME = function(dfInput_pairedWide, formulaInput) {
  environment(formulaInput) <- environment() # Added based on https://stackoverflow.com/questions/60510732/lapply-fitting-thousands-of-mixed-models-and-being-able-to-extract-lsmeans
  
  # Convert categorical variables to factor before fitting model
  dfInput_pairedWide$SUBJECTID <- as.factor(dfInput_pairedWide$SUBJECTID)
  # Specify effects coding for the age factor
  dfInput_pairedWide$age <- as.factor(dfInput_pairedWide$age)
  contrasts(dfInput_pairedWide$age) <- contr.Sum(levels(dfInput_pairedWide$age)) 
  # Specify effects coding for the maternal sensitivity factor
  dfInput_pairedWide$mSens <- as.factor(dfInput_pairedWide$mSens)
  contrasts(dfInput_pairedWide$mSens) <- contr.Sum(levels(dfInput_pairedWide$mSens)) 
  
  # Fit LME model
  fit.LME <- lmer(formulaInput, data = dfInput_pairedWide, REML = TRUE)
  
  tryCatch({
    # If the model does not have problems (i.e., non-convergence or singular fit)
    if (check_convergence(fit.LME) && !check_singularity(fit.LME)) {
      # Calculate marginal means for each maternal sensitivity group
      mLME_mSens <- emmeans(fit.LME, pairwise ~ mSens, mode = "satterthwaite",
                            lmerTest.limit = 240000, at = list(presentNumberAvg = presentNumberRef))
      # Extract maternal sensitivity marginal means
      LME_output_mSens_EM <- data.frame(summary(mLME_mSens)$emmeans)[,1:3]
      
      # Calculate marginal means for each age group
      mLME_age <- emmeans(fit.LME, pairwise ~ age, mode = "satterthwaite",
                          lmerTest.limit = 240000, at = list(presentNumberAvg = presentNumberRef))
      # Extract age marginal means
      LME_output_age_EM <- data.frame(summary(mLME_age)$emmeans)[,1:3]
      
      # Create dataframe with model predictors and other information 
      return(data.frame(lowMSensEMM = LME_output_mSens_EM[2,2], lowMSensSE = LME_output_mSens_EM[2,3],
                        highMSensEMM = LME_output_mSens_EM[1,2], highMSensSE = LME_output_mSens_EM[1,3],
                        mSensEffect = unname(fixef(fit.LME)[3]),
                        mSensPValue = summary(fit.LME)$coefficients[3,5],
                        youngEMM = LME_output_age_EM[2,2], youngSE = LME_output_age_EM[2,3],
                        oldEMM = LME_output_age_EM[1,2], oldSE = LME_output_age_EM[1,3],
                        ageEffect = unname(fixef(fit.LME)[4]), 
                        agePValue = summary(fit.LME)$coefficients[4,5],
                        presentNumEffect = unname(fixef(fit.LME)[2]), 
                        presentNumPValue = summary(fit.LME)$coefficients[2,5],
                        subjectIntercept = data.frame(VarCorr(fit.LME))[1,5],
                        actorIntercept = NA, # Not relevant for non-Interaction models
                        trialN = 2*nobs(fit.LME),
                        subjectN = nrow(ranef(fit.LME)$SUBJECTID),
                        notConvergeN = 0, singFitN = 0, modelType = "LME"))
    }
    # Check for model problems and return model output accordingly
    else if (!check_convergence(fit.LME)) {
      return(data.frame(lowMSensEMM = NA, lowMSensSE = NA,
                        highMSensEMM = NA, highMSensSE = NA,
                        mSensEffect = NA,
                        mSensPValue = NA,
                        youngEMM = NA, youngSE = NA,
                        oldEMM = NA, oldSE = NA,
                        ageEffect = NA,
                        agePValue = NA,
                        presentNumEffect = NA,
                        presentNumPValue = NA,
                        subjectIntercept = data.frame(VarCorr(fit.LME))[1,5],
                        actorIntercept = NA, # Not relevant for non-Interaction models
                        trialN = 2*nobs(fit.LME),
                        subjectN = nrow(ranef(fit.LME)$SUBJECTID),
                        notConvergeN = 1, singFitN = 0, modelType = "LME"))
    }
    else if (check_singularity(fit.LME)) {
      return(data.frame(lowMSensEMM = NA, lowMSensSE = NA,
                        highMSensEMM = NA, highMSensSE = NA,
                        mSensEffect = NA,
                        mSensPValue = NA,
                        youngEMM = NA, youngSE = NA,
                        oldEMM = NA, oldSE = NA,
                        ageEffect = NA,
                        agePValue = NA,
                        presentNumEffect = NA,
                        presentNumPValue = NA,
                        subjectIntercept = data.frame(VarCorr(fit.LME))[1,5],
                        actorIntercept = NA, # Not relevant for non-Interaction models
                        trialN = 2*nobs(fit.LME),
                        subjectN = nrow(ranef(fit.LME)$SUBJECTID),
                        notConvergeN = 0, singFitN = 1, modelType = "LME"))
    }
  }, 
  error=function(e) {
    if (e$message == "Lapack routine dgesv: system is exactly singular: U[1,1] = 0") {
      return(data.frame(lowMSensEMM = NA, lowMSensSE = NA,
                        highMSensEMM = NA, highMSensSE = NA,
                        mSensEffect = NA,
                        mSensPValue = NA,
                        youngEMM = NA, youngSE = NA,
                        oldEMM = NA, oldSE = NA,
                        ageEffect = NA,
                        agePValue = NA,
                        presentNumEffect = NA,
                        presentNumPValue = NA,
                        subjectIntercept = data.frame(VarCorr(fit.LME))[1,5],
                        actorIntercept = NA, # Not relevant for non-Interaction models
                        trialN = 2*nobs(fit.LME),
                        subjectN = nrow(ranef(fit.LME)$SUBJECTID),
                        notConvergeN = -99, singFitN = 1, modelType = "LME")) # -99 designate this specific singular fit error
    }
    else{
      return(data.frame(lowMSensEMM = NA, lowMSensSE = NA,
                        highMSensEMM = NA, highMSensSE = NA,
                        mSensEffect = NA,
                        mSensPValue = NA,
                        youngEMM = NA, youngSE = NA,
                        oldEMM = NA, oldSE = NA,
                        ageEffect = NA,
                        agePValue = NA,
                        presentNumEffect = NA,
                        presentNumPValue = NA,
                        subjectIntercept = data.frame(VarCorr(fit.LME))[1,5],
                        actorIntercept = NA, # Not relevant for non-Interaction models
                        trialN = 2*nobs(fit.LME),
                        subjectN = nrow(ranef(fit.LME)$SUBJECTID),
                        notConvergeN = -99, singFitN = -99, modelType = "LME")) # -99 designates other type of error
    }
  }) 
}