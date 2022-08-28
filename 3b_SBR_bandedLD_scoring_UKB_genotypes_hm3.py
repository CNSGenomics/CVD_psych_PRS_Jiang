### This script is used to score individuals in the UKB to generate PRSs
import os
import numpy as np

chromosome_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']

def scoring_genotype(in_disease, in_SBR_file):
	print(in_disease)

	### Creating the sub directories
	cmd = 'mkdir /scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/' + in_disease + '_SBR_bandedLD/' + in_disease + '_UKB_genotype_scoring/'
	os.system(cmd)

	### Creating a snplist file containing only the SNPs used to generate the PRS
	in_snplist = in_SBR_file.replace('.snpRes', '.snplist')
	print('Creating a snplist file containing only the SNPs used to generate the PRS')
	cmd = "awk 'NR!=1 {print $2}' " + in_SBR_file + " > " + in_snplist
	os.system(cmd)

	### Scoring individual chromosomes
	print('Scoring individual chromosomes')
	## Using plink1.9 you can just add the individual scores
	for chromosome in chromosome_list:
		print('scoring for chr' + chromosome)
		cmd = 'plink \
			--bfile genotyped_files \
			--remove withdrawn_ID.list \
			--score ' + in_SBR_file + ' 2 5 8 header sum \
			--extract ' + in_snplist + ' \
			--out /scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/' + in_disease + '_SBR_bandedLD/' + in_disease + '_UKB_genotype_scoring/' + in_disease + '_SBR_bandedLD_ukbEURu_imp_plink_sum_chr' + chromosome
		os.system(cmd)

	# print('Adding up the score for individual chromosomes. Using plink1.9 you can just add the individual scores. Not with plink2 though')

	combined_score_list = []
	for chromosome in chromosome_list:
		print(chromosome)

		curr_chr_score_file = '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/' + in_disease + '_SBR_bandedLD/' + in_disease + '_UKB_genotype_scoring/' + in_disease + '_SBR_bandedLD_ukbEURu_imp_plink_sum_chr' + chromosome + '.profile'

		if chromosome == '1':
			FID_list_chr1 = []
			IID_list_chr1 = []
			Sample_score_list = []

			with open (curr_chr_score_file, 'r') as a:
				lines = a.readlines()[1:]
				for line in lines:
					line = line.strip()
					fields = line.split()

					FID_list_chr1.append(fields[0])
					IID_list_chr1.append(fields[1])
					Sample_score_list.append(float(fields[5]))

			print('length check, should be the same: ', len(IID_list_chr1), len(Sample_score_list))
			IID_chr1_string = '|'.join(IID_list_chr1)
			Sample_score_array = np.asarray(Sample_score_list)
			combined_score_list.append(Sample_score_array)

		else: 
			IID_list = []
			Sample_score_list = []

			with open (curr_chr_score_file, 'r') as b:
				lines = b.readlines()[1:]
				for line in lines:
					line = line.strip()
					fields = line.split()

					IID_list.append(fields[1])
					Sample_score_list.append(float(fields[5]))

			print('length check, should be the same: ', len(IID_list), len(Sample_score_list))
			IID_string = '|'.join(IID_list)
			if IID_string == IID_chr1_string:
				Sample_score_array = np.asarray(Sample_score_list)
				combined_score_list.append(Sample_score_array)
			else:
				print('IID not matching')


	print('Number of chromosomes with scores recorded, should be 22: ', len(combined_score_list))

	scores_sum_list = np.sum(combined_score_list, axis = 0)

	out_file = '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/' + in_disease + '_SBR_bandedLD/' + in_disease + '_UKB_genotype_scoring/' + in_disease + '_SBR_bandedLD_ukbEURu_imp_plink_sum_all_sum.profile'
	with open (out_file, 'w') as c:
		c.write('\t'. join(['FID', 'IID', 'Scores']) + '\n')
		i = 0
		while i < len(scores_sum_list):
			curr_line = '\t'. join([FID_list_chr1[i], IID_list_chr1[i], str(scores_sum_list[i])])
			c.write(curr_line + '\n')
			i = i + 1


	print('Z-transforming the sum scores')
	scores_mean = np.mean(scores_sum_list)
	scores_std = np.std(scores_sum_list)
	scores_scale_list = (scores_sum_list - scores_mean)/scores_std

	out_scale_file = '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/' + in_disease + '_SBR_bandedLD/' + in_disease + '_UKB_genotype_scoring/' + in_disease + '_SBR_bandedLD_ukbEURu_imp_plink_sum_all_sum_ztrans.profile'
	with open (out_scale_file, 'w') as c:
		c.write('\t'. join(['FID', 'IID', 'Scores', 'Scores_scale']) + '\n')
		i = 0
		while i < len(scores_sum_list):
			curr_line = '\t'. join([FID_list_chr1[i], IID_list_chr1[i], str(scores_sum_list[i]), str(scores_scale_list[i])])
			c.write(curr_line + '\n')
			i = i + 1


scoring_genotype('AF', '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/AF_SBR_bandedLD/AF_UKB_SNPs_SBR_bandedLD.snpRes')
scoring_genotype('BD', '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/BD_SBR_bandedLD/BD_UKB_SNPs_SBR_bandedLD.snpRes')
scoring_genotype('CAD', '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/CAD_SBR_bandedLD/CAD_UKB_SNPs_SBR_bandedLD.snpRes')
scoring_genotype('HF', '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/HF_SBR_bandedLD/HF_UKB_SNPs_SBR_bandedLD.snpRes')
scoring_genotype('MD', '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/MD_SBR_bandedLD/MD_UKB_SNPs_SBR_bandedLD.snpRes')
scoring_genotype('SCZ', '/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/SCZ_SBR_bandedLD/SCZ_UKB_SNPs_SBR_bandedLD.snpRes')




