#!/bin/bash
#
#PBS -l select=1:ncpus=10:mem=72GB
#PBS -l walltime=48:00:00
#PBS -N SBR_banded
#PBS -e /home/uqjjian3/jobs_sterr_stout/
#PBS -o /home/uqjjian3/jobs_sterr_stout/
#PBS -A UQ-IMB-CNSG
#PBS -J 1-6

# Declare an array of string with type
declare -a TraitArray=("AF" "BD" "CAD" "HF" "MD" "SCZ")
declare -a FileArray=("/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/UKB_SNPs_filtered_sumstats/AF_ukbEURu_imp_v3_HM3_autosome_SNPs_filtered" "/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/UKB_SNPs_filtered_sumstats/BD_ukbEURu_imp_v3_HM3_autosome_SNPs_filtered" "/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/UKB_SNPs_filtered_sumstats/CAD_ukbEURu_imp_v3_HM3_autosome_SNPs_filtered" "/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/UKB_SNPs_filtered_sumstats/HF_ukbEURu_imp_v3_HM3_autosome_SNPs_filtered" "/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/UKB_SNPs_filtered_sumstats/MD_ukbEURu_imp_v3_HM3_autosome_SNPs_filtered" "/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/UKB_SNPs_filtered_sumstats/SCZ_ukbEURu_imp_v3_HM3_autosome_SNPs_filtered")
declare -a OutFileArray=("/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/AF_SBR_bandedLD/AF_UKB_SNPs_SBR_bandedLD" "/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/BD_SBR_bandedLD/BD_UKB_SNPs_SBR_bandedLD" "/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/CAD_SBR_bandedLD/CAD_UKB_SNPs_SBR_bandedLD" "/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/HF_SBR_bandedLD/HF_UKB_SNPs_SBR_bandedLD" "/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/MD_SBR_bandedLD/MD_UKB_SNPs_SBR_bandedLD" "/scratch/user/uqjjian3/CVD_Psych_PRS_project/63_checkpoint_24Aug2022/3_SBR_generation/SCZ_SBR_bandedLD/SCZ_UKB_SNPs_SBR_bandedLD")

i=$((PBS_ARRAY_INDEX-1))

echo ${i}
echo "${TraitArray[i]}"
echo "${FileArray[i]}"
echo "${OutFileArray[i]}"

/home/uqjjian3/utils/GCTB/gctb_2.03beta_Linux/gctb --sbayes R \
--mldm /QRISdata/Q4319/data/GCTB_LD_matrices/band_ukb_10k_hm3/ukb10k.mldm \
--gwas-summary "${FileArray[i]}" \
--pi 0.95,0.02,0.02,0.01 \
--gamma 0.0,0.01,0.1,1 \
--exclude-mhc \
--hsq 0.5 \
--chain-length 50000 \
--burn-in 10000 \
--seed 12345 \
--thread 10 \
--no-mcmc-bin \
--out "${OutFileArray[i]}" 2>&1










