
##################################################################
### This script is used to compute individual PRS using SBayesR
##################################################################

SBR_dir="" # This is the output directory that contains the SBayesR outputs
score_dir="" # This is the output directory that contains the individual scores

geno_dir="" # This is the directory that contains the individual genotype file
LD_matric_dir="" # This is the directory that contains the LD matrices downloaded from the SBayesR website
trait="" # This is the name of the trait
sumstats="" # This is the GWAS summary statistic file that contains only the summary statitistics for SNPs profiled in the prediction cohort
withdrawn_failedQC_ID="" # This is the IDs of individuals who have withdrawn or failed QC

#########################################################
echo "Weight estimation" 

gctb --sbayes R \
--mldm "$LD_matric_dir"band_ukb_10k_hm3/ukb10k.mldm \
--gwas-summary "$sumstats" \
--pi 0.95,0.02,0.02,0.01 \
--gamma 0.0,0.01,0.1,1 \
--exclude-mhc \
--hsq 0.5 \
--chain-length 50000 \
--burn-in 10000 \
--seed 12345 \
--thread 30 \
--no-mcmc-bin \
--out "$SBR_dir""$trait"_UKB_SNPs_SBR_bandedLD 2>&1

#########################################################
echo "Scoring genotypes"

### Getting a list of the snps from the SBR output
awk 'NR!=1 {print $2}' "$SBR_dir""$trait"_UKB_SNPs_SBR_bandedLD.snpRes >  "$SBR_dir""$trait"_UKB_SNPs_SBR_bandedLD.snplist

### Scoring each autosome individually
for chromosome in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22"
do
	plink \
			--bfile "$geno_dir"genotype_file_chr"$chromosome" \
			--remove "$withdrawn_failedQC_ID" \
			--score "$SBR_dir""$trait"_UKB_SNPs_SBR_bandedLD.snpRes 2 5 8 header sum \
			--extract "$SBR_dir""$trait"_UKB_SNPs_SBR_bandedLD.snpRes.snplist \
			--out "$score_dir""$trait"_score_chr"$chromosome"
done


