
##################################################################
### This script is used to calculate prediction accuracy a PRS (i.e. OR and AUC)
##################################################################

library(data.table)
library(stringr)
library(rcompanion)
library(pROC)

PRS_pheno_file <- '' # This is the file that contains the individual phenotype and PRS information

##################################################
### Writing a function for performing logistic regression
##################################################
evaluate_PRS_prediction_UKB <- function(in_outcome, in_predictor, in_PRS_pheno_file) {
	print(paste0('Evaluating the prediction accuracy for ', in_predictor))

	### Reading the pheno and PRS file
	in_PRS_pheno_df <- fread(in_PRS_pheno_file)

	### Based on the outcome, removing the 'NA' individuals so they are not accounted for in the analysis
	in_PRS_pheno_df <- in_PRS_pheno_df[in_PRS_pheno_df[[in_outcome]] != 'NA', ]

	### Summarising the count of cases and controls
	print(table(in_PRS_pheno_df[[in_outcome]]))

	### Full model with the predictor
	glm_formula <- as.formula(paste(c(in_outcome, ' ~ ', in_predictor, ' + factor(Sex) + recruitment_age + factor(array) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))
	glm_results <- glm(glm_formula, data = in_PRS_pheno_df, family = binomial(link = "logit"))
	print(summary(glm_results))

	### AUC value for full model
	auc_value <- auc(in_PRS_pheno_df[[in_outcome]], glm_results$linear.predictors)
	print(auc_value)
}


#################################################
## Calculating the prediction accurary
#################################################
evaluate_PRS_prediction_UKB('AF_prevalent', 'AF_PRS', PRS_pheno_file)
evaluate_PRS_prediction_UKB('CAD_prevalent', 'CAD_PRS', PRS_pheno_file)
evaluate_PRS_prediction_UKB('HF_prevalent', 'HF_PRS', PRS_pheno_file)
evaluate_PRS_prediction_UKB('MD_prevalent', 'MD_PRS', PRS_pheno_file)
evaluate_PRS_prediction_UKB('SCZ_prevalent', 'SCZ_PRS', PRS_pheno_file)
evaluate_PRS_prediction_UKB('BD_prevalent', 'BD_PRS', PRS_pheno_file)

















	######################## Running the glm ########################
	### NULL model
	print('Evaluating null model, adjusting for only sex, age, genotyping array and 20 genetic PCs')
	null_minimum_formula <- as.formula(paste(c(in_outcome, ' ~ factor(Sex) + recruitment_age + factor(array) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))
	glm_fit_null_minimum <- glm(null_minimum_formula, data = in_PRS_pheno_df, family = binomial(link = "logit"))
	glm_fit_null_minimum_results_file <- paste(c(work_dir, in_outcome, '_UKB_null_minimum_model_glm.txt'), collapse ='')
	sink(glm_fit_null_minimum_results_file) # Start writing to an output file
	print(summary(glm_fit_null_minimum))
	sink() # Stop writing to the file
