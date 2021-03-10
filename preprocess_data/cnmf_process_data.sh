#/bin/bash


FILE=$1
if [ ! -f "$FILE" ]; then
    echo "Please enter a valid path for the peak names of sample 1."
    echo "Call cnmf_process_data using the following format:"
    echo "bash cnmf_process_data.sh path/to/peak.bed hg19 path/to/prior/data"
	exit 2
fi

if [ -z "$2" ];
  then
    echo "No argument supplied genome version. Please select one of following: {hg19, hg38, mm9, mm10}."
    echo "Call cnmf_process_data using the following format:"
    echo "bash cnmf_process_data.sh path/to/peak.bed hg19 path/to/prior/data"
    exit 2
fi

if [ ! $2 == "hg19" ] && [ ! $2 == "hg38" ] && [ ! $2 == "mm9" ] && [ ! $2 == "mm10" ];
then 
    echo "Please select one of following: {hg19, hg38, mm9, mm10}.'"
    echo "Call cnmf_process_data using the following format:"
    echo "bash cnmf_process_data.sh path/to/peak.bed hg19 path/to/prior/data"
    exit 2
fi 

dir=$3
if [ ! -d "$dir" ]; then
    echo "Please enter a valid path prior data."
    echo "Call cnmf_process_data using the following format:"
    echo "bash cnmf_process_data.sh path/to/peak.bed hg19 path/to/prior/data"
	exit 2
fi

genome=$2

sort -k1,1 -k2,2n $1 > peak_sorted.bed
bedtools intersect -a $3/Promoter_100k_${genome}.bed  -b peak_sorted.bed -wa -wb -sorted|awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$4,$3-100000-$6}'|sed 's/-//g'|sort -k1,1 -k2,2n >peak_gene_100k.bed
bedtools intersect -a peak_gene_100k.bed -b $3/RE_gene_corr_${genome}.bed -wa -wb -sorted|awk 'BEGIN{OFS="\t"}{if ($4==$9) print $1,$2,$3,$4,$5,$10}'|sed 's/\t/\_/1'|sed 's/\t/\_/1'>peak_gene_100k_corr.bed