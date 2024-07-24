
##################################################################
### This script is used to perform the mediation analysis
##################################################################

library(data.table)
library(stringr)
library(mediation)

##################################################
### Importing the supplied arguments
##################################################
work_dir <- '' # This is the output directory
sim_set <- 10000 # This is the number of simulations
in_predictor <- '' # This is the PRS
in_risk_factor <- '' # This is the risk factor
in_outcome <- '' This is the CVD outcome 
in_seed <- 123 # This is a seeding number
in_sex <- '' # female or male 
in_covariates_c <- c('baseline_BMI', 'recruitment_age', 'array', 'baseline_smoked_ever', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19', 'PC20')

### Removing BMI from the covariates if the mediator investigated is BMI
if (in_risk_factor == 'baseline_BMI'){
	in_covariates_c <- in_covariates_c[in_covariates_c != 'baseline_BMI']
}

### Removing smoking from the covariates if the mediator investigated is smoking
if (in_risk_factor == 'baseline_smoked_ever'){
	in_covariates_c <- in_covariates_c[in_covariates_c != 'baseline_smoked_ever']
}

in_covariates <- paste(in_covariates_c, collapse = ' + ')

##################################################
### Importing the pheno, age and PRS dataframe
##################################################
PRS_pheno_file <- '' # This is the file that contains the individual phenotype and PRS information
PRS_pheno_df <- fread(PRS_pheno_file)

# Removing missing data
curr_PRS_pheno_df <- PRS_pheno_df[which(!is.na(PRS_pheno_df[[in_outcome]]) & !is.na(PRS_pheno_df[[in_risk_factor]])), ]

# Filter for sex
sex_pheno_PRS_df <- curr_PRS_pheno_df[which(curr_PRS_pheno_df$Sex == in_sex), ] # apply the necessary filtering steps for the sex of interest

#############################################################################################################
### Step #1: The total effect
### The total effect describes the total effect that the PRS has on the CVD risk. 
##############################################################################################################
print('Step 1')
PRS_CVD_formula <- as.formula(paste(c(in_outcome,' ~ ', in_predictor, ' + ', in_covariates), collapse = ''))
PRS_CVD_association <- glm(PRS_CVD_formula, data = sex_pheno_PRS_df, family=binomial(link="logit"))		

##############################################################################################################
### Step #2: The effect of the PRS onto the risk factor
##############################################################################################################
print('Step 2')
PRS_RiskFactor_formula <- as.formula(paste(c(in_risk_factor,' ~ ', in_predictor, ' + ', in_covariates), collapse = ''))
if (in_risk_factor == 'baseline_BMI'){
	PRS_RiskFactor_association <- lm(PRS_RiskFactor_formula, data = sex_pheno_PRS_df)
} else {
	PRS_RiskFactor_association <- glm(PRS_RiskFactor_formula, data = sex_pheno_PRS_df, family=binomial(link="logit"))		
}

##############################################################################################################
### Step #3: The effect of the PRS and risk factor on the CVD risks
##############################################################################################################
print('Step 3')
Predictor_outome_formula <- as.formula(paste(c(in_outcome,' ~ ', in_predictor, ' + ', in_risk_factor, ' + ', in_covariates), collapse = ''))
Predictor_outome_association <- glm(Predictor_outome_formula, data = sex_pheno_PRS_df, family=binomial(link="logit"))		

##############################################################################################################
### Step #4: Causal Mediation Analysis
##############################################################################################################
print('Step 4 : mediation analysis')
set.seed(in_seed)

mediation_output <- mediate(PRS_RiskFactor_association, Predictor_outome_association, treat = in_predictor, mediator = in_risk_factor, boot = TRUE, sims = sim_set)
summary(mediation_output)


