
##################################################################
### This script is to perform a cox regression analysis 
##################################################################
library(data.table)
library(stringr)
library(survival)
library(knitr)

# work_dir <- '/home/uqjjian3/CVD_Psych_PRS_project/89_checkpoint_24Jul2023/6_Cox_regression/'

##################################################
### Importing the pheno, age and PRS dataframe
##################################################
PRS_pheno_file <- '' # This is the file that contains the individual phenotype and PRS information

full_IID_pheno_PRS_df <- fread(PRS_pheno_file)

##################################################
### Writing a function for performing sex-stratified cox regression
##################################################
running_coxph <- function(in_outcome, in_outcome_age, in_predictor, in_sex) {

	### Stratifying the df by sex if necessary
	print('Stratifying by sex')
	if (in_sex == 'female') {
		IID_pheno_PRS_df <- full_IID_pheno_PRS_df[full_IID_pheno_PRS_df$Sex == '0', ]
		} else {IID_pheno_PRS_df <- full_IID_pheno_PRS_df[full_IID_pheno_PRS_df$Sex == '1', ]}

	### Based on the outcome, removing the 'NA' individuals so they are not accounted for in the analysis
	IID_pheno_PRS_df <- IID_pheno_PRS_df[IID_pheno_PRS_df[[in_outcome]] != 'NA', ]

	######################## Running coxph ########################
	cox_fit_formula <- as.formula(paste(c('Surv(', in_outcome_age, ', ', in_outcome, ') ~ ', in_predictor, ' + factor(array) + baseline_BMI + factor(baseline_smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))
		
	coxph_fit <- coxph(cox_fit_formula, data = IID_pheno_PRS_df)
	coxph_fit_results_file <- paste(c(work_dir, in_subdir, in_outcome, '_', in_outcome_age, '_', in_sex, '_', in_cohort, '_UKB_', in_predictor_str, '_coxph.txt'), collapse ='')
	summary(coxph_fit)


}

##################################################
### Writing a function for performing cox regression, stratified by menopause
##################################################
running_coxph_meno <- function(in_outcome, in_outcome_age, in_predictor) {

	### Based on the outcome, removing the 'NA' individuals so they are not accounted for in the analysis
	curr_pheno_PRS_df <- curr_pheno_PRS_df[curr_pheno_PRS_df[[in_outcome]] != 'NA', ]

	######################## Running coxph ########################
	cox_fit_formula <- as.formula(paste(c('Surv(', in_outcome_age, ', ', in_outcome, ') ~ ', in_predictor, ' + factor(array) + baseline_BMI + factor(baseline_smoked_ever) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))

	### For females who have had menopause
	post_curr_pheno_PRS_df <- curr_pheno_PRS_df[((curr_pheno_PRS_df$Sex == '0') & (curr_pheno_PRS_df$baseline_menopause == '1')), ]
	post_coxph_fit <- coxph(cox_fit_formula, data = post_curr_pheno_PRS_df)
	summary(post_coxph_fit)

	### For females who have not had menopause at baseline and were < 50 at baseline
	pre_curr_pheno_PRS_df <- curr_pheno_PRS_df[((curr_pheno_PRS_df$Sex == '0') & (curr_pheno_PRS_df$baseline_menopause == '0')), ]
	pre_coxph_fit <- coxph(cox_fit_formula, data = pre_curr_pheno_PRS_df)
	summary(pre_coxph_fit)
}
