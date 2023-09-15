
##################################################################
### This script is used to calculate the nagerkerke R2 for a PRS
##################################################################
library(data.table)
library(stringr)
library(rcompanion)
library(pROC)

PRS_pheno_file <- '' # This is the file that contains the individual phenotype and PRS information

running_glm_prediction_accuracy <- function(in_outcome, in_predictor, in_PRS_pheno_file) {
	### Reading the pheno and PRS file
	full_IID_pheno_PRS_df <- fread(in_PRS_pheno_file)

	### Based on the outcome, removing the 'NA' individuals so they are not accounted for in the analysis
	IID_pheno_PRS_df <- full_IID_pheno_PRS_df[full_IID_pheno_PRS_df[[in_outcome]] != 'NA', ]

	######################## Running the glm ########################
	### NULL model
	print('Evaluating null model')
	null_formula <- as.formula(paste(c(in_outcome, ' ~ factor(Sex) + recruitment_age + factor(array) + baseline_BMI + factor(baseline_smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))
	glm_fit_null <- glm(null_formula, data = IID_pheno_PRS_df, family = binomial(link = "logit"))

	### Full model with the predictor
	print('Evaluating the full model')
	full_formula <- as.formula(paste(c(in_outcome, ' ~ ', in_predictor, ' + factor(Sex) + recruitment_age + factor(array) + baseline_BMI + factor(baseline_smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))
	glm_fit_full <- glm(full_formula, data = IID_pheno_PRS_df, family = binomial(link = "logit"))
	summary(glm_fit_full)

	### Nagelkerke R2 for full model
	nagelkerke(glm_fit_full, null = glm_fit_null, restrictNobs = FALSE)

}
