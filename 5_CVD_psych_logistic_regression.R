### This script is used to perform logistic regression between all three psych and all three CVDs in the entire cohort
### module load R/4.1.0+sf

library(data.table)
library(stringr)

PRS_pheno_file <- '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/5_Assembling_PRS_pheno/UKB_CVD_and_psych_pheno_PRS_covariates.txt' ### This is the file that contains the dataframe with the PRSs and UK Biobank IDs for each CVD ('1'-case, '0'-not a case)
work_dir <- '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/7_logistic_regression_analysis/'

##################################################
### Writing a function for performing logistic regression
##################################################
running_glm <- function(in_outcome, in_predictor, in_predictor_str, in_sex, in_cohort, in_PRS_pheno_file, in_subdir) {
	### Reading the pheno and PRS file
	full_IID_pheno_PRS_df <- fread(in_PRS_pheno_file)

	print('running glm for')
	print(in_outcome)
	print(paste(c(in_outcome, in_predictor, in_sex), collapse = ' | '))

	### Stratifying the df by sex if necessary
	print('Stratifying by sex')
	if (in_sex == 'combined') {
	IID_pheno_PRS_df <- full_IID_pheno_PRS_df} else {
		if (in_sex == 'female') {
			IID_pheno_PRS_df <- full_IID_pheno_PRS_df[full_IID_pheno_PRS_df$Sex == '0', ]
			} else {IID_pheno_PRS_df <- full_IID_pheno_PRS_df[full_IID_pheno_PRS_df$Sex == '1', ]}}

	print('Filtering by cohort')
	### Filtering to retain screened control individuals (individuals who have no diagnosis of psychiatric disorders or are on any psychiatric medications) if necessary
	if (in_cohort == 'screened_control') {
		IID_pheno_PRS_df <- IID_pheno_PRS_df[IID_pheno_PRS_df$MD_prevalent == '0', ]
	} else {
		IID_pheno_PRS_df <- IID_pheno_PRS_df
	}

	### Based on the outcome, removing the 'NA' individuals so they are not accounted for in the analysis
	IID_pheno_PRS_df <- IID_pheno_PRS_df[IID_pheno_PRS_df[[in_outcome]] != 'NA', ]

	### Summarising the outcome (the count of cases and controls)
	summary_df <- as.data.frame(table(IID_pheno_PRS_df[[in_outcome]]))
	colnames(summary_df) <- c('status', 'count')
	summary_file <- paste(c(work_dir, in_subdir, in_outcome, '_', in_sex, '_', in_cohort, '_UKB_', in_predictor_str, '_phenotype_summary.txt'), collapse ='')
	write.table(summary_df, summary_file, quote = FALSE, sep = '\t', row.names = FALSE)

	######################## Running the glm ########################
	print('Evaluating the full model')
	if (in_sex == 'combined') {
		full_formula <- as.formula(paste(c(in_outcome, ' ~ ', in_predictor, ' + factor(Sex) + recruitment_age + factor(array) + mean_BMI + factor(Smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))} else {
			full_formula <- as.formula(paste(c(in_outcome, ' ~ ', in_predictor, ' + recruitment_age + factor(array) + mean_BMI + factor(Smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))}
		
	glm_fit_full <- glm(full_formula, data = IID_pheno_PRS_df, family = binomial(link = "logit"))
	glm_fit_full_results_file <- paste(c(work_dir, in_subdir, in_outcome, '_', in_sex, '_', in_cohort, '_UKB_', in_predictor_str, '_full_model_glm.txt'), collapse ='')
	sink(glm_fit_full_results_file) # Start writing to an output file
	print(summary(glm_fit_full))
	sink() # Stop writing to the file

}


#################################################
## Performing logistic regression of CVD incidence on psychiatric disorders PRS in all cohort
#################################################

############################## With death cases included as incidence cases ##############################
### Combined for sex
running_glm('CAD_incident_withDeath', 'MD_PRS', 'MD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_psych_whole_cohort/')
running_glm('HF_incident_withDeath', 'MD_PRS', 'MD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_psych_whole_cohort/')
running_glm('AF_incident_withDeath', 'MD_PRS', 'MD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_psych_whole_cohort/')

running_glm('CAD_incident_withDeath', 'SCZ_PRS', 'SCZ_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_psych_whole_cohort/')
running_glm('HF_incident_withDeath', 'SCZ_PRS', 'SCZ_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_psych_whole_cohort/')
running_glm('AF_incident_withDeath', 'SCZ_PRS', 'SCZ_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_psych_whole_cohort/')

running_glm('CAD_incident_withDeath', 'BD_PRS', 'BD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_psych_whole_cohort/')
running_glm('HF_incident_withDeath', 'BD_PRS', 'BD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_psych_whole_cohort/')
running_glm('AF_incident_withDeath', 'BD_PRS', 'BD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_psych_whole_cohort/')

###  sex-stratified
running_glm('CAD_incident_withDeath', 'MD_PRS', 'MD_PRS', 'female', 'all_cohort', PRS_pheno_file, 'CVD_MD_sex_stratified/')
running_glm('HF_incident_withDeath', 'MD_PRS', 'MD_PRS', 'female', 'all_cohort', PRS_pheno_file, 'CVD_MD_sex_stratified/')
running_glm('AF_incident_withDeath', 'MD_PRS', 'MD_PRS', 'female', 'all_cohort', PRS_pheno_file, 'CVD_MD_sex_stratified/')

running_glm('CAD_incident_withDeath', 'MD_PRS', 'MD_PRS', 'male', 'all_cohort', PRS_pheno_file, 'CVD_MD_sex_stratified/')
running_glm('HF_incident_withDeath', 'MD_PRS', 'MD_PRS', 'male', 'all_cohort', PRS_pheno_file, 'CVD_MD_sex_stratified/')
running_glm('AF_incident_withDeath', 'MD_PRS', 'MD_PRS', 'male', 'all_cohort', PRS_pheno_file, 'CVD_MD_sex_stratified/')

##################################################
### Sensitivity analysis
##################################################
############################## With death cases removed ##############################
### combined
running_glm('CAD_incident_withoutDeath', 'MD_PRS', 'MD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')
running_glm('HF_incident_withoutDeath', 'MD_PRS', 'MD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')
running_glm('AF_incident_withoutDeath', 'MD_PRS', 'MD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')

### female
running_glm('CAD_incident_withoutDeath', 'MD_PRS', 'MD_PRS', 'female', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')
running_glm('HF_incident_withoutDeath', 'MD_PRS', 'MD_PRS', 'female', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')
running_glm('AF_incident_withoutDeath', 'MD_PRS', 'MD_PRS', 'female', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')

### male
running_glm('CAD_incident_withoutDeath', 'MD_PRS', 'MD_PRS', 'male', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')
running_glm('HF_incident_withoutDeath', 'MD_PRS', 'MD_PRS', 'male', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')
running_glm('AF_incident_withoutDeath', 'MD_PRS', 'MD_PRS', 'male', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')

############################## With death cases included, but also including CVD PRS ##############################
### combined
running_glm('CAD_incident_withDeath', 'MD_PRS + CAD_PRS', 'MD_PRS_CAD_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')
running_glm('HF_incident_withDeath', 'MD_PRS + HF_PRS', 'MD_PRS_HF_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')
running_glm('AF_incident_withDeath', 'MD_PRS + AF_PRS', 'MD_PRS_AF_PRS', 'combined', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')

### female
running_glm('CAD_incident_withDeath', 'MD_PRS + CAD_PRS', 'MD_PRS_CAD_PRS', 'female', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')
running_glm('HF_incident_withDeath', 'MD_PRS + HF_PRS', 'MD_PRS_HF_PRS', 'female', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')
running_glm('AF_incident_withDeath', 'MD_PRS + AF_PRS', 'MD_PRS_AF_PRS', 'female', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')

### male
running_glm('CAD_incident_withDeath', 'MD_PRS + CAD_PRS', 'MD_PRS_CAD_PRS', 'male', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')
running_glm('HF_incident_withDeath', 'MD_PRS + HF_PRS', 'MD_PRS_HF_PRS', 'male', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')
running_glm('AF_incident_withDeath', 'MD_PRS + AF_PRS', 'MD_PRS_AF_PRS', 'male', 'all_cohort', PRS_pheno_file, 'CVD_MD_sensitivity/')


##################################################
### Performing logistic regression of CVD incidence on psychiatric disorders PRS in psych controls
##################################################
### combined
running_glm('CAD_incident_withDeath', 'MD_PRS', 'MD_PRS', 'combined', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('HF_incident_withDeath', 'MD_PRS', 'MD_PRS', 'combined', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('AF_incident_withDeath', 'MD_PRS', 'MD_PRS', 'combined', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')

### female
running_glm('CAD_incident_withDeath', 'MD_PRS', 'MD_PRS', 'female', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('HF_incident_withDeath', 'MD_PRS', 'MD_PRS', 'female', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('AF_incident_withDeath', 'MD_PRS', 'MD_PRS', 'female', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')

### male
running_glm('CAD_incident_withDeath', 'MD_PRS', 'MD_PRS', 'male', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('HF_incident_withDeath', 'MD_PRS', 'MD_PRS', 'male', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('AF_incident_withDeath', 'MD_PRS', 'MD_PRS', 'male', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')

##################################################
### Performing logistic regression of CVD prevalence on psychiatric disorders PRS in individuals with no diangosis of psychiatric disorders or are on any psychiatric medications
##################################################
### combined
running_glm('CAD_prevalent_withDeath', 'MD_PRS', 'MD_PRS', 'combined', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('HF_prevalent_withDeath', 'MD_PRS', 'MD_PRS', 'combined', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('AF_prevalent_withDeath', 'MD_PRS', 'MD_PRS', 'combined', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')

### female
running_glm('CAD_prevalent_withDeath', 'MD_PRS', 'MD_PRS', 'female', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('HF_prevalent_withDeath', 'MD_PRS', 'MD_PRS', 'female', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('AF_prevalent_withDeath', 'MD_PRS', 'MD_PRS', 'female', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')

### male
running_glm('CAD_prevalent_withDeath', 'MD_PRS', 'MD_PRS', 'male', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('HF_prevalent_withDeath', 'MD_PRS', 'MD_PRS', 'male', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('AF_prevalent_withDeath', 'MD_PRS', 'MD_PRS', 'male', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')

############################## sensitivity analysis where CVD PRS is included ##############################
### combined
running_glm('CAD_incident_withDeath', 'MD_PRS + CAD_PRS', 'MD_PRS_CAD_PRS', 'combined', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('HF_incident_withDeath', 'MD_PRS + HF_PRS', 'MD_PRS_HF_PRS', 'combined', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('AF_incident_withDeath', 'MD_PRS + AF_PRS', 'MD_PRS_AF_PRS', 'combined', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')

### female
running_glm('CAD_incident_withDeath', 'MD_PRS + CAD_PRS', 'MD_PRS_CAD_PRS', 'female', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('HF_incident_withDeath', 'MD_PRS + HF_PRS', 'MD_PRS_HF_PRS', 'female', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('AF_incident_withDeath', 'MD_PRS + AF_PRS', 'MD_PRS_AF_PRS', 'female', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')

### male
running_glm('CAD_incident_withDeath', 'MD_PRS + CAD_PRS', 'MD_PRS_CAD_PRS', 'male', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('HF_incident_withDeath', 'MD_PRS + HF_PRS', 'MD_PRS_HF_PRS', 'male', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('AF_incident_withDeath', 'MD_PRS + AF_PRS', 'MD_PRS_AF_PRS', 'male', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')

### combined
running_glm('CAD_prevalent_withDeath', 'MD_PRS + CAD_PRS', 'MD_PRS_CAD_PRS', 'combined', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('HF_prevalent_withDeath', 'MD_PRS + HF_PRS', 'MD_PRS_HF_PRS', 'combined', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('AF_prevalent_withDeath', 'MD_PRS + AF_PRS', 'MD_PRS_AF_PRS', 'combined', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')

### female
running_glm('CAD_prevalent_withDeath', 'MD_PRS + CAD_PRS', 'MD_PRS_CAD_PRS', 'female', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('HF_prevalent_withDeath', 'MD_PRS + HF_PRS', 'MD_PRS_HF_PRS', 'female', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('AF_prevalent_withDeath', 'MD_PRS + AF_PRS', 'MD_PRS_AF_PRS', 'female', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')

### male
running_glm('CAD_prevalent_withDeath', 'MD_PRS + CAD_PRS', 'MD_PRS_CAD_PRS', 'male', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('HF_prevalent_withDeath', 'MD_PRS + HF_PRS', 'MD_PRS_HF_PRS', 'male', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
running_glm('AF_prevalent_withDeath', 'MD_PRS + AF_PRS', 'MD_PRS_AF_PRS', 'male', 'screened_control', PRS_pheno_file, 'CVD_MD_screened_control/')
