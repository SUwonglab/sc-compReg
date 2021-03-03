#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 2-00:00:00          # Runtime in D-HH:MM
#SBATCH -p whwong          # Partition to submit to
#SBATCH --mem=125000         # Memory pool for all cores
#SBATCH -o Run-job_2.out     # File to which STDOUT will be written
#SBATCH -e Run-job_2.err     # File to which STDERR will be written


genome=hg19;
####Peak intersect
cat PeakName1.txt|tr '_' '\t' |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}'> PeakName1.bed
cat PeakName2.txt|tr '_' '\t' |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}'> PeakName2.bed
bedtools intersect -a PeakName1.bed -b PeakName2.bed -wa -wb|sort -k1,1 -k2,2n|cut -f 4,8 > PeakName_intersect.txt
cat PeakName_intersect.txt|cut -f 1|tr '_' '\t'|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}'> PeakName_1.bed
bedtools intersect -a PeakName1.bed -b PeakName_1.bed -wa -v >PeakName1_non_overlap.bed
bedtools intersect -a PeakName2.bed -b PeakName_1.bed -wa -v >PeakName2_non_overlap.bed
cat PeakName_1.bed PeakName1_non_overlap.bed PeakName2_non_overlap.bed |sort -k1,1 -k2,2n > PeakName.bed
cat PeakName.bed|cut -f 1-3 > PeakName_3columns.bed
##peak-gene prior
bedtools intersect -a /oak/stanford/groups/whwong/group_scratch/zduren/openness/singleCell/CoupledNMF_v2/common_data/Promoter_100k_${genome}.bed  -b PeakName_3columns.bed -wa -wb -sorted|awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$4,$3-100000-$6}'|sed 's/-//g'|sort -k1,1 -k2,2n >peak_gene_100k_intersect.bed
bedtools intersect -a peak_gene_100k_intersect.bed -b /oak/stanford/groups/whwong/group_scratch/zduren/openness/singleCell/CoupledNMF_v2/common_data/RE_gene_corr_${genome}.bed -wa -wb -sorted|awk 'BEGIN{OFS="\t"}{if ($4==$9) print $1,$2,$3,$4,$5,$10}'|sed 's/\t/\_/1'|sed 's/\t/\_/1'>peak_gene_prior_intersect.bed
####motif scan
findMotifsGenome.pl PeakName.bed ${genome} ./. -size given -find /oak/stanford/groups/whwong/group_scratch/zduren/openness/human/immunotherapy/CAR_Tcell/PECA2_modified/Data/all_motif_rmdup -preparsedDir ./Homer/ > MotifTarget.bed
cat MotifTarget.bed|awk 'NR>1'|cut -f 1,4,6 > MotifTarget.txt
rm MotifTarget.bed
rm motifFindingParameters.txt
