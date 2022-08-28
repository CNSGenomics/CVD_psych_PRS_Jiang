### This script is used to compare the log(OR) coefficients from logistic regression betwee female and male models, using a wald test. The equation (b_f-b_m)/(sqrt(se_f^2+se_m^2)) is found on page 1276 on https://www.journals.uchicago.edu/doi/epdf/10.1086/230638 (Statistical Methods for Comparing Regression Coefficients Between Models by Clifford C. Clogg, Eva Petkova, and Adamantios Haritou)
import os
import math
import scipy.stats

### A function to extract the logistic regression estimates and SE values
def getting_estimates(in_log_file):
	with open (in_log_file, 'r') as a:
		line = a.readlines()[12]
		line = line.strip()
		fields = line.split()
		beta = fields[1]
		se = fields[2]
	return([beta, se])

### A function to perform wald test
def wald_log(in_dir):
	### Looping through the directory where the logisic regression results are
	for root, dirs, files in os.walk(in_dir):
		for file in files:
			if ('female' in file) and file.endswith('_full_model_glm.txt'):
				print(file)

				### Conjugating the full path of the file, as well as the paths to the female and male files
				female_path = os.path.join(root, file)
				male_path = female_path.replace('female', 'male')

				### Getting the estimates and se for female and male
				female_beta = float(getting_estimates(female_path)[0])
				female_se = float(getting_estimates(female_path)[1])
				male_beta = float(getting_estimates(male_path)[0])
				male_se = float(getting_estimates(male_path)[1])

				### Calculating z_score
				curr_z = (female_beta - male_beta)/(math.sqrt((female_se**2) + (male_se**2)))
				curr_p = scipy.stats.norm.sf(abs(curr_z))*2

				if curr_p < 0.05:
					print('two-tailed Wald test p < 0.05 :', curr_p)


wald_log('/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/7_logistic_regression_analysis/CVD_MD_sex_stratified/')
wald_log('/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/7_logistic_regression_analysis/CVD_MD_sensitivity/')
wald_log('/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/7_logistic_regression_analysis/CVD_MD_screened_control/')
wald_log('/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/8_CVD_MD_logistic_regression_risk_group/')



      