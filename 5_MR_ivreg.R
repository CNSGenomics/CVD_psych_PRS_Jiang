
##################################################################
### This script is to perform a sex-specific 2sls MR analysis for CVD on MD. 
### Note: the results of a MR performed on binary exposure should be interpreted with EXTRA caution: https://link.springer.com/article/10.1007/s10654-018-0424-6
##################################################################

library(data.table)
library(stringr)
library(broom)
library(AER)

##################################################
### Importing the pheno, age and PRS dataframe
##################################################
PRS_pheno_file <- '' # This is the file that contains the individual phenotype and PRS information

### Reading the data
PRS_pheno_df <- fread(PRS_pheno_file)

##################################################
### Writing a function for performing the ivreg analysis
##################################################
running_ivreg <- function(in_outcome, in_exposure, in_predictor, in_sex) {

	print('Running ivreg for')
	print(paste(c(in_outcome, in_exposure, in_predictor, in_sex, in_cohort), collapse = ' | '))

	### Stratifying the df by sex if necessary
	print('Stratifying by sex')
	if (in_sex == 'female') {
		sex_pheno_PRS_df <- PRS_pheno_df[which(PRS_pheno_df$Sex == '0'), ]
		} else {sex_pheno_PRS_df <- PRS_pheno_df[which(PRS_pheno_df$Sex == '1'), ]}

	### Based on the outcome, removing the 'NA' individuals so they are not accounted for in the analysis
	sex_pheno_PRS_df <- sex_pheno_PRS_df[which((sex_pheno_PRS_df[[in_outcome]] != 'NA') & (sex_pheno_PRS_df[[in_exposure]] != 'NA')), ]

	######################## Running coxph ########################
	ivreg_formula <- as.formula(paste(c(in_outcome, ' ~ ', in_exposure, ' + recruitment_age + array + baseline_BMI + baseline_smoked_ever + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 | ', in_predictor, '  + recruitment_age + array + baseline_BMI + baseline_smoked_ever + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20'), collapse = ''))
	ivreg_result <- summary(ivreg(ivreg_formula, data = sex_pheno_PRS_df))
	print(ivreg_result)
}


for (in_sex in c('female', 'male')){

	running_ivreg('AF_prevalent', 'MD_prevalent', 'MD_PRS', in_sex)
	running_ivreg('CAD_prevalent', 'MD_prevalent', 'MD_PRS', in_sex)
	running_ivreg('HF_prevalent', 'MD_prevalent', 'MD_PRS', in_sex)

}
