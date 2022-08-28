### This script is used to calculate the nagerkerke R2 and the AUC for each PRS

library(data.table)
library(stringr)
library(rcompanion)
library(pROC)

PRS_pheno_file <- '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/5_Assembling_PRS_pheno/UKB_CVD_and_psych_pheno_PRS_covariates.txt' ### This is the file that contains the dataframe with the PRSs and UK Biobank IDs for each CVD ('1'-case, '0'-not a case)
work_dir <- '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/6_PRS_prediction_accuracy/'

##############################
### Writing a function for performing logistic regression
##################################################
running_glm_prediction_accuracy <- function(in_outcome, in_predictor, in_PRS_pheno_file) {
	### Reading the pheno and PRS file
	full_IID_pheno_PRS_df <- fread(in_PRS_pheno_file)

	### Based on the outcome, removing the 'NA' individuals so they are not accounted for in the analysis
	IID_pheno_PRS_df <- full_IID_pheno_PRS_df[full_IID_pheno_PRS_df[[in_outcome]] != 'NA', ]

	### Summarising the outcome (the count of cases and controls)
	summary_df <- as.data.frame(table(IID_pheno_PRS_df[[in_outcome]]))
	colnames(summary_df) <- c('status', 'count')
	summary_file <- paste(c(work_dir, in_outcome, '_UKB_', in_predictor, '_phenotype_summary.txt'), collapse ='')
	write.table(summary_df, summary_file, quote = FALSE, sep = '\t', row.names = FALSE)

	######################## Running the glm ########################
	### NULL model
	print('Evaluating null model')
	null_formula <- as.formula(paste(c(in_outcome, ' ~ factor(Sex) + recruitment_age + factor(array) + mean_BMI + factor(Smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))

	glm_fit_null <- glm(null_formula, data = IID_pheno_PRS_df, family = binomial(link = "logit"))
	R2_results_file <- paste(c(work_dir, in_outcome, '_UKB_null_model_glm.txt'), collapse ='')
	sink(R2_results_file) # Start writing to an output file
	print(summary(glm_fit_null))
	sink() # Stop writing to the file

	### Full model with the predictor
	print('Evaluating the full model')
	full_formula <- as.formula(paste(c(in_outcome, ' ~ ', in_predictor, ' + factor(Sex) + recruitment_age + factor(array) + mean_BMI + factor(Smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))

	glm_fit_full <- glm(full_formula, data = IID_pheno_PRS_df, family = binomial(link = "logit"))
	glm_fit_full_results_file <- paste(c(work_dir, in_outcome, '_UKB_', in_predictor, '_full_model_glm.txt'), collapse ='')
	sink(glm_fit_full_results_file) # Start writing to an output file
	print(summary(glm_fit_full))
	sink() # Stop writing to the file

	### Nagelkerke R2 for full model
	R2_results_file <- paste(c(work_dir, in_outcome, '_UKB_', in_predictor, '_full_model_glm_nagelkerke_R2.txt'), collapse ='')
	sink(R2_results_file) # Start writing to an output file
	print(nagelkerke(glm_fit_full, null = glm_fit_null, restrictNobs = FALSE))
	sink() # Stop writing to the file

	### AUC value for full model
	auc_value = auc(IID_pheno_PRS_df[[in_outcome]], glm_fit_full$linear.predictors)
	AUC_results_file <- paste(c(work_dir, in_outcome, '_UKB_', in_predictor, '_full_model_glm_pROC_AUC.txt'), collapse ='')
	sink(AUC_results_file) # Start writing to an output file
	print(print(auc_value))
	sink() # Stop writing to the file

}

#################################################
## Getting the prediction accuracy (Nagelkelke's R2) for each PRS using the prevalence information for the corresponding trait
#################################################
running_glm_prediction_accuracy('AF_prevalent_withDeath', 'AF_PRS', PRS_pheno_file)
running_glm_prediction_accuracy('CAD_prevalent_withDeath', 'CAD_PRS', PRS_pheno_file)
running_glm_prediction_accuracy('HF_prevalent_withDeath', 'HF_PRS', PRS_pheno_file)
running_glm_prediction_accuracy('MD_prevalent', 'MD_PRS', PRS_pheno_file)
running_glm_prediction_accuracy('SCZ_prevalent', 'SCZ_PRS', PRS_pheno_file)
running_glm_prediction_accuracy('BD_prevalent', 'BD_PRS', PRS_pheno_file)
