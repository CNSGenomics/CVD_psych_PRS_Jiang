
##################################################################
### This script is to perform a cox regression analysis 
##################################################################

library(data.table)
library(stringr)
library(survival)
library(knitr)

##################################################
### Importing the pheno, age and PRS dataframe
##################################################
PRS_pheno_file <- '' # This is the file that contains the individual phenotype and PRS information

### Reading the data
PRS_pheno_df <- fread(PRS_pheno_file)

##################################################
### Writing a function for performing sex-stratified cox regression
##################################################
running_coxph <- function(in_outcome, in_outcome_age, in_predictor, in_sex) {
	print(paste0('Running Cox regression for ', in_outcome, '|', in_outcome_age, '|', in_predictor, '|', in_sex))

	### Stratifying the df by sex if necessary
	print('Stratifying by sex')
	if (in_sex == 'female') {
		sex_pheno_PRS_df <- PRS_pheno_df[which(PRS_pheno_df$Sex == '0'), ]
		} else {sex_pheno_PRS_df <- PRS_pheno_df[which(PRS_pheno_df$Sex == '1'), ]}

	### Based on the outcome, removing the 'NA' individuals so they are not accounted for in the analysis
	sex_pheno_PRS_df <- sex_pheno_PRS_df[sex_pheno_PRS_df[[in_outcome]] != 'NA', ]

	### Running coxph
	cox_fit_formula <- as.formula(paste(c('Surv(', in_outcome_age, ', ', in_outcome, ') ~ ', in_predictor, ' + factor(array) + baseline_BMI + factor(baseline_smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))
		
	coxph_fit <- coxph(cox_fit_formula, data = sex_pheno_PRS_df)
	print(summary(coxph_fit))


}

running_coxph('AF_incident', 'AF_time_to_event', 'MD_PRS', 'female')
running_coxph('AF_incident', 'AF_time_to_event', 'MD_PRS', 'male')
running_coxph('CAD_incident', 'CAD_time_to_event', 'MD_PRS', 'female')
running_coxph('CAD_incident', 'CAD_time_to_event', 'MD_PRS', 'male')
running_coxph('HF_incident', 'HF_time_to_event', 'MD_PRS', 'female')
running_coxph('HF_incident', 'HF_time_to_event', 'MD_PRS', 'male')

##################################################
### Writing a function for performing cox regression, stratified by menopause
##################################################
running_coxph_meno <- function(in_outcome, in_outcome_age, in_predictor) {

	### Based on the outcome, removing the 'NA' individuals so they are not accounted for in the analysis
	PRS_pheno_df <- PRS_pheno_df[PRS_pheno_df[[in_outcome]] != 'NA', ]

	######################## Running coxph ########################
	cox_fit_formula <- as.formula(paste(c('Surv(', in_outcome_age, ', ', in_outcome, ') ~ ', in_predictor, ' + factor(array) + baseline_BMI + factor(baseline_smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))

	### For females who have had menopause
	post_curr_pheno_PRS_df <- PRS_pheno_df[((PRS_pheno_df$Sex == '0') & (PRS_pheno_df$baseline_menopause == '1')), ]
	post_coxph_fit <- coxph(cox_fit_formula, data = post_curr_pheno_PRS_df)
	summary(post_coxph_fit)

	### For females who have not had menopause at baseline and were < 50 years of age at baseline
	pre_curr_pheno_PRS_df <- PRS_pheno_df[((PRS_pheno_df$Sex == '0') & (PRS_pheno_df$baseline_menopause == '0')), ]
	pre_coxph_fit <- coxph(cox_fit_formula, data = pre_curr_pheno_PRS_df)
	print(summary(pre_coxph_fit))
}

running_coxph_meno('AF_incident', 'AF_time_to_event', 'MD_PRS')
running_coxph_meno('CAD_incident', 'CAD_time_to_event', 'MD_PRS')
running_coxph_meno('HF_incident', 'HF_time_to_event', 'MD_PRS')

