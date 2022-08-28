### This script is used to perform logistic regression between all three psych and all three CVDs in the entire cohort
### module load R/4.1.0+sf

library(data.table)
library(stringr)

PRS_pheno_file <- '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/5_Assembling_PRS_pheno/UKB_CVD_and_psych_pheno_PRS_covariates.txt' ### This is the file that contains the dataframe with the PRSs and UK Biobank IDs for each CVD ('1'-case, '0'-not a case)
work_dir <- '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/8_CVD_MD_logistic_regression_risk_group/'

##################################################
### Writing a function for performing logistic regression for the top risk groups
##################################################
running_glm_for_risk_group <- function(in_outcome, in_predictor, in_sex, in_cohort, in_PRS_pheno_file, in_low_prob, in_high_prob, in_thresh_label){
	print('Running logistic regression analysis for ')
	print(paste(c(in_outcome, in_predictor, in_sex, in_cohort), collapse = ' | '))

	### The dataframe comtaining the PRS and the incident data, along with the covariates.
	full_IID_pheno_PRS_df <- fread(in_PRS_pheno_file)

	### Stratifying the df by sex if necessary
	print('Stratifying by sex')
	if (in_sex == 'combined') {
	IID_pheno_PRS_df <- full_IID_pheno_PRS_df} else {
		if (in_sex == 'female') {
			IID_pheno_PRS_df <- full_IID_pheno_PRS_df[full_IID_pheno_PRS_df$Sex == '0', ]
			} else {IID_pheno_PRS_df <- full_IID_pheno_PRS_df[full_IID_pheno_PRS_df$Sex == '1', ]}}

	print('Filtering by cohort')
	### Filtering to retain screened control individuals (individuals with no diagnosis of psychiatric disorders or are on any psychiatric medications) if necessary
	if (in_cohort == 'screened_control') {
		IID_pheno_PRS_df <- IID_pheno_PRS_df[IID_pheno_PRS_df$MD_prevalent == '0', ]
	} else {
		IID_pheno_PRS_df <- IID_pheno_PRS_df
	}

	### Based on the outcome, removing the 'NA' individuals so they are not accounted for in the analysis
	IID_pheno_PRS_df <- IID_pheno_PRS_df[IID_pheno_PRS_df[[in_outcome]] != 'NA', ]

	####################################################
	### Finding high-risk and low-risk PRS individuals
	####################################################
	### Getting the top and bottom z-scores
	PRS_top_perc_z <- as.numeric(quantile(IID_pheno_PRS_df[[in_predictor]], probs = c(in_high_prob), names = FALSE))
	PRS_bottom_perc_z <- as.numeric(quantile(IID_pheno_PRS_df[[in_predictor]], probs = c(in_low_prob), names = FALSE))

	### Annotating PRS group.
	risk_group_col <- paste(c(in_predictor, '_risk_group'), collapse = '')
	IID_pheno_PRS_df[[risk_group_col]] <- ifelse(IID_pheno_PRS_df[[in_predictor]] > PRS_top_perc_z, 'Top', 'NA')
	IID_pheno_PRS_df[[risk_group_col]] <- ifelse(IID_pheno_PRS_df[[in_predictor]] < PRS_bottom_perc_z, 'Bottom', IID_pheno_PRS_df[[risk_group_col]])
	IID_pheno_PRS_df <- IID_pheno_PRS_df[IID_pheno_PRS_df[[risk_group_col]] != 'NA', ]
	
	### Summarising the number of ID in each risk group
	PRS_risk_group_summary_file <- paste(c(work_dir, in_outcome, '_', in_predictor, '_', in_sex, '_', in_cohort, '_', in_thresh_label, '_risk_group_phenotype_summary.txt'), collapse = '')
	sink(PRS_risk_group_summary_file) # Start writing to an output file
	print(table(IID_pheno_PRS_df[[risk_group_col]]))
	sink() # Stop writing to the file	

	####################################################
	#### Running the logistic regression analysis
	####################################################
	print('Running the logistic regression analysis')
	### NULL model
	print('Evaluating null model')
	if (in_sex == 'combined') {
		null_formula <- as.formula(paste(c(in_outcome, ' ~ factor(Sex) + recruitment_age + factor(array) + mean_BMI + factor(Smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))} else {
			null_formula <- as.formula(paste(c(in_outcome, ' ~ recruitment_age + factor(array) + mean_BMI + factor(Smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))}
	
	glm_fit_null <- glm(null_formula, data = IID_pheno_PRS_df, family = binomial(link = "logit"))
	glm_fit_null_results_file <- paste(c(work_dir, in_outcome, '_', in_sex, '_', in_cohort, '_', in_thresh_label, '_risk_group_null_model_glm.txt'), collapse ='')
	sink(glm_fit_null_results_file) # Start writing to an output file
	print(summary(glm_fit_null))
	sink() # Stop writing to the file

	### Full model with the predictor
	print('Evaluating the full model')
	if (in_sex == 'combined') {
		full_formula <- as.formula(paste(c(in_outcome, ' ~ factor(', risk_group_col, ') + factor(Sex) + recruitment_age + factor(array) + mean_BMI + factor(Smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))} else {
			full_formula <- as.formula(paste(c(in_outcome, ' ~ factor(', risk_group_col, ') + recruitment_age + factor(array) + mean_BMI + factor(Smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))}
		
	glm_fit_full <- glm(full_formula, data = IID_pheno_PRS_df, family = binomial(link = "logit"))
	Psych_model_results_file <- paste(c(work_dir, in_outcome, '_', in_predictor, '_', in_sex, '_', in_cohort, '_', in_thresh_label, '_risk_group_full_model_glm.txt'), collapse ='')
	sink(Psych_model_results_file) # Start writing to an output file
	print(summary(glm_fit_full))
	sink() # Stop writing to the file

}



##################################################
### Performing logistic regression of CVD incidence on MD PRS top risk groups
##################################################
############################## With death cases included as incidence cases ##############################
### top 25% groups
### 25% threshold for defining risk group
running_glm_for_risk_group('CAD_incident_withDeath', 'MD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 0.25, 0.75, '25perc')
running_glm_for_risk_group('CAD_incident_withDeath', 'MD_PRS', 'male', 'all_cohort', PRS_pheno_file, 0.25, 0.75, '25perc')
running_glm_for_risk_group('CAD_incident_withDeath', 'MD_PRS', 'female', 'all_cohort', PRS_pheno_file, 0.25, 0.75, '25perc')

running_glm_for_risk_group('HF_incident_withDeath', 'MD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 0.25, 0.75, '25perc')
running_glm_for_risk_group('HF_incident_withDeath', 'MD_PRS', 'male', 'all_cohort', PRS_pheno_file, 0.25, 0.75, '25perc')
running_glm_for_risk_group('HF_incident_withDeath', 'MD_PRS', 'female', 'all_cohort', PRS_pheno_file, 0.25, 0.75, '25perc')

running_glm_for_risk_group('AF_incident_withDeath', 'MD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 0.25, 0.75, '25perc')
running_glm_for_risk_group('AF_incident_withDeath', 'MD_PRS', 'male', 'all_cohort', PRS_pheno_file, 0.25, 0.75, '25perc')
running_glm_for_risk_group('AF_incident_withDeath', 'MD_PRS', 'female', 'all_cohort', PRS_pheno_file, 0.25, 0.75, '25perc')

### 10% threshold for defining risk group
running_glm_for_risk_group('CAD_incident_withDeath', 'MD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 0.1, 0.9, '10perc')
running_glm_for_risk_group('CAD_incident_withDeath', 'MD_PRS', 'male', 'all_cohort', PRS_pheno_file, 0.1, 0.9, '10perc')
running_glm_for_risk_group('CAD_incident_withDeath', 'MD_PRS', 'female', 'all_cohort', PRS_pheno_file, 0.1, 0.9, '10perc')

running_glm_for_risk_group('HF_incident_withDeath', 'MD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 0.1, 0.9, '10perc')
running_glm_for_risk_group('HF_incident_withDeath', 'MD_PRS', 'male', 'all_cohort', PRS_pheno_file, 0.1, 0.9, '10perc')
running_glm_for_risk_group('HF_incident_withDeath', 'MD_PRS', 'female', 'all_cohort', PRS_pheno_file, 0.1, 0.9, '10perc')

running_glm_for_risk_group('AF_incident_withDeath', 'MD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 0.1, 0.9, '10perc')
running_glm_for_risk_group('AF_incident_withDeath', 'MD_PRS', 'male', 'all_cohort', PRS_pheno_file, 0.1, 0.9, '10perc')
running_glm_for_risk_group('AF_incident_withDeath', 'MD_PRS', 'female', 'all_cohort', PRS_pheno_file, 0.1, 0.9, '10perc')



