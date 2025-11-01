# Difference Wave Simulation Script: 5. Compare Model Performance 

# This script calculates the number of usable trials and subjects, statistical
# power to detect the maternal sensitivity effect, and error and bias of the 
# maternal sensitivity estimates for each difference wave approach. 

# ***See Section3_SimulatedData README.md available on the LME_MixedEffectsERPDifferenceWave
# GitHub for pipeline details: https://github.com/basclab/LME_MixedEffectsERPDifferenceWave/tree/main/Section3_SimulatedData

# Requirements: 
  # - Needs R Version 3.6.1 and packages listed below
  # - importFilename: File containing the estimated marginal means for each 
  #   maternal sensitivity group and other model predictors across all models and
  #   simulated datasets (samples) for one type of decay and missingness pattern 
  #   (e.g., different decay/Missingness Pattern #1). See ModelOutput_DataDictionary.xlsx 
  #   file for further information.
  # - See instructions in data environment for specifying other variables

# Script Functions:
  # 1. Load model output file
  # 2. Create final dataframe for analysis
  # 3. Characterize LME model problems
  # 4. Calculate number of usable trials and subjects
  # 5. Calculate power to detect maternal sensitivity effect
  # 6. Calculate error and bias of high maternal sensitivity group estimates
  # 7. Calculate error and bias of low maternal sensitivity group estimates
  # 8. Visualize marginal means for each maternal sensitivity group

# Outputs: 
  # - Tables summarizing number of trials/subjects, statistical power, error, and
  #   bias for each difference wave approach
  # - Visualizations of marginal means for each maternal sensitivity group 

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
library(dplyr) # V.1.0.2; %>% operator
library(ggplot2) # V.3.3.2; graphing functions
library(rstatix) # V.0.6.0; pairwise_t_test function
library(coxed) # V.0.3.3; bootstrapped CI calculation
library(stringr) # V.1.4.9; str_remove function

#-------------------------------------------------------------------------------
# DATA ENVIRONMENT 

# Specify variables needed to import corresponding file 
missType <- 'Miss1'
decayRate <- 'different'
sampleN <- 1000
subjectN <- 48
rpIter <- 10000 
importFolder <- 'C:/Users/basclab/Desktop/Section3_SimulatedData/03_ModelOutput/'
importFilename <- paste0(importFolder, 'ModelOutput_', missType, '_', decayRate, 'Decay_sampleN', sampleN, '_subN', subjectN, '_rpIter', rpIter, '.csv')

# Calculate population value of maternal sensitivity effect of interest at average 
# presentation number simulated in data (this value is specific to the decay rate) 
if (decayRate == 'same') {
  lowMSens <- 1.997
  highMSens <- -1.997
  
} else if (decayRate == 'different') {
  lowMSens <- mean(seq(1.997, 1.997+(1.000*9), by = 1.000))
  highMSens <- mean(seq(-1.997, -1.997+(1.000*9), by = 1.000))
}

#-------------------------------------------------------------------------------
# 1. LOAD MODEL OUTPUT FILE

# Load specified model output file for missingness pattern and decay rate
modelOutput <- fread(importFilename) 

# Create modelMatch column using modelType (e.g., LME) and matchMethod (e.g., Int) columns
modelOutput$modelMatch <- paste0(modelOutput$modelType, '_', modelOutput$matchMethod)
unique(modelOutput$modelMatch)

# Format modelMatch column with labels for graphing 
modelOutput$modelMatch <- factor(modelOutput$modelMatch, levels = c('LME_Int', 'LME_EM', 'LME_RP' , 'LME_NN:N','LME_NN:P','LME_NN:F', 'ANOVA_DW'),
                                 labels = c('Interact', 'Exact M', 'Rand P', 'NN:N', 'NN:P', 'NN:F', 'ANOVA'))

# Format percentage of low trial-count labels for graphing 
modelOutput$caseDeletionPct <- factor(modelOutput$caseDeletionPct, levels = c('Pop.','0%','6%','11%','32%'),
                                      labels = c('Population', '0% Low Trial-Count', '6% Low Trial-Count',
                                                 '11% Low Trial-Count', '32% Low Trial-Count')) 

# Check for occurrence of models with specific problems and then update corresponding variables:
# - notConvergeN = -99 & singFitN = 1: Special singular fit error (Lapack routine dgesv: system is exactly singular: U[1,1] = 0)
spSingFitN_ind <- which(modelOutput$notConvergeN == -99 & modelOutput$singFitN == 1)
modelOutput$notConvergeN[spSingFitN_ind] = 1
# - notConvergeN = -99 & singFitN = -99: Other type of error
#   This did not occur for any simulated sample so no updates were needed 
otherError_ind <- which(modelOutput$notConvergeN == -99 & modelOutput$singFitN == -99)

# Update usable trial and subject count to NA if model had problems 
modelOutput$trialN[modelOutput$notConvergeN == 1 | modelOutput$singFitN == 1] <- NA
modelOutput$subjectN[modelOutput$notConvergeN == 1 | modelOutput$singFitN == 1] <- NA

#-------------------------------------------------------------------------------
# 2. CREATE FINAL DATAFRAME FOR ANALYSIS

# Extract output for all models except for LME Random Permutation and calculate 
# significant effect columns
modelOutput_noRP <- modelOutput %>%
  filter(!(modelMatch %in% c('Rand P'))) %>%
  mutate(mSensPValue_sig = ifelse(mSensPValue <= 0.05, TRUE, FALSE),
         agePValue_sig = ifelse(agePValue <= 0.05, TRUE, FALSE),
         presentNumPValue_sig = ifelse(presentNumPValue <= 0.05, TRUE, FALSE))

# Extract output for LME Random Permutation model and aggregate for each sample 
# (e.g., summarize across 10,000 estimated maternal sensitivity estimates so 
# there is one average maternal sensitivity estimate per simulated sample)
modelOutput_RP <- modelOutput %>% 
  filter(modelMatch %in% c('Rand P') & !is.na(mSensEffect)) %>% # Select Random Permutation model output for iterations in which the model converged (i.e., maternal sensitivity estimate is not NA)
  group_by(modelMatch, caseDeletionPct, sample, modelType, matchMethod) %>% # Aggregate values within each percentage of low trial-count subjects and simulated sample (other variables are included so that they will be present in output)
  summarise(lowMSensEMM_agg = mean(lowMSensEMM, na.rm = T),
            lowMSensSE_agg = sd(lowMSensEMM, na.rm = T),
            highMSensEMM_agg = mean(highMSensEMM, na.rm = T),
            highMSensSE_agg = sd(highMSensEMM, na.rm = T),
            mSensEffect_agg = mean(mSensEffect, na.rm = T),
            mSensPValue = NA, # Set to NA for RP models
            mSensPValue_sig = !(dplyr::between(0, bca(mSensEffect)[1], bca(mSensEffect)[2])), # Examine if confidence intervals include 0 (TRUE if they do not)
            youngEMM_agg = mean(youngEMM, na.rm = T),
            youngSE_agg = sd(youngEMM, na.rm = T),
            oldEMM_agg = mean(oldEMM, na.rm = T),
            oldSE_agg = sd(oldEMM, na.rm = T),
            ageEffect_agg = mean(ageEffect, na.rm = T),
            agePValue = NA, 
            agePValue_sig = !(dplyr::between(0, bca(ageEffect)[1], bca(ageEffect)[2])),
            presentNumEffect_agg = mean(presentNumEffect, na.rm = T),
            presentNumPValue = NA,
            presentNumPValue_sig = !(dplyr::between(0, bca(presentNumEffect)[1], bca(presentNumEffect)[2])),
            subjectIntercept = mean(subjectIntercept, na.rm = T), 
            actorIntercept = NA, # Not relevant for non-Interaction models
            trialN = mean(trialN, na.rm = T), 
            subjectN = mean(subjectN, na.rm = T), 
            notConvergeN = NA, # Non-convergence is summarized separately below 
            singFitN = NA, # Singular fit is summarized separately below 
            sampleRP = max(sampleRP, na.rm = T)) %>%
  ungroup() %>%
  # Rename columns so that we can merge with model output from non-Random Permutation models
  rename(lowMSensEMM = lowMSensEMM_agg,
         lowMSensSE = lowMSensSE_agg,
         highMSensEMM = highMSensEMM_agg,
         highMSensSE = highMSensSE_agg,
         mSensEffect = mSensEffect_agg,
         youngEMM = youngEMM_agg,
         youngSE = youngSE_agg,
         oldEMM = oldEMM_agg,
         oldSE = oldSE_agg,
         ageEffect = ageEffect_agg,
         presentNumEffect = presentNumEffect_agg) 

# Create final dataframe with model output for all difference wave approaches  
modelOutput_final <- rbind(modelOutput_noRP, modelOutput_RP)

#-------------------------------------------------------------------------------
# 3. CHARACTERIZE LME MODEL PROBLEMS 

# Extract model problems for all models except LME Random Permutation
modelOutput_noRP_modPr <- modelOutput_noRP[,c("modelMatch", "caseDeletionPct", "sample",
                                              "modelType", "matchMethod", "notConvergeN",
                                              "singFitN", "sampleRP")]

# Extract model problems for LME Random Permutation model and aggregate for each 
# sample 
modelOutput_RP_modPr <- modelOutput %>% 
  filter(modelMatch %in% c('Rand P')) %>%
  group_by(modelMatch, caseDeletionPct, sample, modelType, matchMethod) %>%
  summarise(notConvergeN = sum(notConvergeN, na.rm = T), # Number of models that did not converge 
            singFitN = sum(singFitN, na.rm = T), # Number of models that had singular fit 
            sampleRP = max(sampleRP, na.rm = T)) %>% # Max random permutation iteration (not used for analysis)
  ungroup() 

# Create final dataframe with model problems for all difference wave approaches  
modelOutput_final_modPr <- rbind(modelOutput_noRP_modPr, modelOutput_RP_modPr)

# Count the number of iterations and samples that had non-convergence or singular fit
modelProblemTable <- modelOutput_final_modPr %>% 
  group_by(caseDeletionPct, modelMatch) %>%
  summarise(notConvergeN_sum = sum(notConvergeN, na.rm = T),
            singFitN_sum = sum(singFitN, na.rm = T)) 

# Calculate the percent of model problems based on the total number of simulated samples
# (for LME Random Permutation, the percent is based on the total number of iterations*samples)
modelProblemTable_final <- modelProblemTable %>%
  mutate(notConvergeN_pct = round(ifelse(modelMatch == "Rand P", 100*notConvergeN_sum/(sampleN*rpIter), 100*notConvergeN_sum/sampleN),3),
         singFitN_pct = round(ifelse(modelMatch == "Rand P", 100*singFitN_sum/(sampleN*rpIter), 100*singFitN_sum/sampleN),3),
         allProb_pct = notConvergeN_pct + singFitN_pct)

#-------------------------------------------------------------------------------
# 4. CALCULATE NUMBER OF USABLE TRIALS AND SUBJECTS

trialSubNTable <- modelOutput_final %>%
  group_by(modelMatch, caseDeletionPct) %>%
  summarise(trialNPerSub_avg = round(mean(trialN/subjectN, na.rm = T), 2),
            trialNPerSub_sd = round(sd(trialN/subjectN, na.rm = T),2),
            subjectN_avg = round(mean(subjectN, na.rm = T),2),
            subjectN_sd = round(sd(subjectN, na.rm = T),2))

#-------------------------------------------------------------------------------
# 5. CALCULATE POWER TO DETECT MATERNAL SENSITIVITY EFFECT 

# Power = number of datasets in which the model detected a significant maternal 
# sensitivity effect divided by the total number of simulated datasets
sigEffectTable <- modelOutput_final %>% 
  group_by(modelMatch, caseDeletionPct) %>%
  summarise(sigMSens = round(sum(mSensPValue_sig, na.rm = T)/sampleN, 2),
            sigAge = round(sum(agePValue_sig, na.rm = T)/sampleN, 2),
            sigPNum = sum(presentNumPValue_sig, na.rm = T)/sampleN)

#-------------------------------------------------------------------------------
# 6. CALCULATE ERROR AND BIAS OF HIGH MATERNAL SENSITIVITY GROUP ESTIMATES

# Calculate average and standard deviation of marginal means
highMSens_marginalMeans <- modelOutput_final %>% 
  group_by(modelMatch, caseDeletionPct) %>%
  summarise(highMSensMean = mean(highMSensEMM, na.rm = TRUE), 
            highMSensStdDev = sd(highMSensEMM, na.rm = TRUE)) 

# Calculate root mean squared error (RMSE) 
highMSens_MSESummary <- modelOutput_final %>%
  mutate(highMSens_SE = (highMSensEMM-highMSens)^2) %>%
  group_by(modelMatch, caseDeletionPct) %>%
  summarise(RMSE = sqrt(mean(highMSens_SE, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(RMSERounded = round(RMSE, 2))

# Calculate average relative bias = average percentage difference of the model's 
# marginal mean from the population value across simulated datasets
highMSens_relBias <- modelOutput_final %>%
  mutate(highMSens_bias = highMSensEMM-highMSens,
         highMSens_relBias = 100*(highMSens_bias/highMSens))
highMSens_relBiasSummary <- highMSens_relBias %>%
  group_by(modelMatch, caseDeletionPct) %>%
  summarise(avgRelBias = mean(highMSens_relBias, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(avgRelBiasRounded = round(avgRelBias, 2))

# Calculate paired t-test for relative bias values for three best LME models 
# (Interaction, Exact Match, Random Permutation) versus ANOVA
highMSens_relBias_forTTest <- highMSens_relBias %>% 
  filter(modelMatch %in% c('Interact', 'Exact M', 'Rand P', 'ANOVA')) %>% 
  select(modelMatch, caseDeletionPct, sample, highMSens_relBias)
highMSens_pairedTTest <- highMSens_relBias_forTTest %>%
  group_by(caseDeletionPct) %>%
  pairwise_t_test(highMSens_relBias ~ modelMatch, paired = TRUE) %>%
  filter(group2 %in% c('ANOVA')) %>%
  mutate(pRounded = str_remove(round(p, 3), "^0+")) # We export p (not p.adj) because we apply a more severe alpha

#-------------------------------------------------------------------------------
# 7. CALCULATE ERROR AND BIAS OF LOW MATERNAL SENSITIVITY GROUP ESTIMATES

# Calculate average and standard deviation of marginal means
lowMSens_marginalMeans <- modelOutput_final %>% 
  group_by(modelMatch, caseDeletionPct) %>%
  summarise(lowMSensMean = mean(lowMSensEMM, na.rm = TRUE), 
            lowMSensStdDev = sd(lowMSensEMM, na.rm = TRUE)) 

# Calculate RMSE 
lowMSens_MSESummary <- modelOutput_final %>%
  mutate(lowMSens_SE = (lowMSensEMM-lowMSens)^2) %>%
  group_by(modelMatch, caseDeletionPct) %>%
  summarise(RMSE = sqrt(mean(lowMSens_SE, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(RMSERounded = round(RMSE, 2))

# Calculate average relative bias 
lowMSens_relBias <- modelOutput_final %>%
  mutate(lowMSens_bias = lowMSensEMM-lowMSens,
         lowMSens_relBias = 100*(lowMSens_bias/lowMSens))
lowMSens_relBiasSummary <- lowMSens_relBias %>%
  group_by(modelMatch, caseDeletionPct) %>%
  summarise(avgRelBias = mean(lowMSens_relBias, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(avgRelBiasRounded = round(avgRelBias, 2))

# Calculate paired t-test for relative bias values
lowMSens_relBias_forTTest <- lowMSens_relBias %>% 
  filter(modelMatch %in% c('Interact', 'Exact M', 'Rand P', 'ANOVA')) %>% 
  select(modelMatch, caseDeletionPct, sample, lowMSens_relBias)
lowMSens_pairedTTest <- lowMSens_relBias_forTTest %>%
  group_by(caseDeletionPct) %>%
  pairwise_t_test(lowMSens_relBias ~ modelMatch, paired = TRUE) %>%
  filter(group2 %in% c('ANOVA')) %>%
  mutate(pRounded = str_remove(round(p, 3), "^0+"))

#-------------------------------------------------------------------------------
# 8. VISUALIZE MARGINAL MEANS FOR EACH MATERNAL SENSITIVITY GROUP

# In this example, we plot the marginal means for the high maternal sensitivity group
plotTitle <- 'A. Missingness Pattern #1: More Missing Data in Later Trials and in Younger Subjects'
highMSens_min <- -6.5 # Y-axis minimum
highMSens_max <- 13 # Y-axis maximum
highMSens_bias_min <- highMSens-(0.1*highMSens) # -10% bias
highMSens_bias_max <- highMSens-(-0.1*highMSens) # +10% bias
modelColors <- c('#e66101','#e66101','#e66101','#e66101','#e66101','#e66101', '#5e3c99')

# Create tiff plotting device and location to save graph
saveFilename <- paste0(importFolder, missType, '_', decayRate, '_highMSens.tiff')
tiff(file = saveFilename,
     width = 13, 
     height = 5, 
     units = 'in',
     res = 600) 

# Create box plot of marginal means for each difference wave approach
ggplot(modelOutput_final) +
  # Create gray shaded area for acceptable 10% bias
  annotate("rect", xmin= -Inf, xmax= Inf,
           ymin=highMSens_bias_min, ymax=highMSens_bias_max,
           fill="gray", alpha = 0.5) +
  # Plot marginal means (we plot for all LME Random Permutation iterations)
  geom_boxplot(mapping = aes(x=modelMatch, y=highMSensEMM, color=modelMatch),
               size = 0.8, outlier.alpha = 0.7, show.legend = FALSE) +
  # Plot dashed line indicating population value
  geom_hline(yintercept = highMSens, linetype = 2, size = 0.8, alpha = 0.3) +
  facet_grid(cols = vars(caseDeletionPct), scales = "free_x") + theme_bw() +
  # Update plot with LME/ANOVA colors
  scale_color_manual("", values = modelColors) +
  theme(text = element_text(size=16),
        axis.text.x= element_text(color = "black", angle=90, vjust=.5, hjust=1),
        axis.text.y= element_text(color = "black"),
        axis.title.x = element_text(vjust=-0.6),
        plot.title.position = "plot") +
  xlab('Model Type') +
  ylab(paste0('Marginal Means of Diff. Wave\n Mean Amp. (\u00B5V)')) +
  labs(title = plotTitle) +
  scale_y_continuous(breaks=seq(-5,20,5), limits=c(highMSens_min, highMSens_max))
dev.off()