#!/bin/sh

## This script uses bcftools mpileup to call variants from sorted bam files 
## usage (for testing with just one individual):
## sbatch run_bcf_caller.sh

#SBATCH --account=def-lukens
#SBATCH --time=0-00:30:00 ## days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # number of threads
#SBATCH --mem=16000 # requested memory (in MB)
#SBATCH --mail-type=END


module load bcftools
module load samtools

#assign variable
output_file="bcf_results"

#Run bcftools so that results include depths (DP & AD), and genotype quality GQ, while excluding indels
echo "- Calling variants for all bam files..."
bcftools mpileup -a DP,AD -f burbot_2021.fasta *.sorted.bam | bcftools call -m --variants-only --format-fields GQ --skip-variants indels > $output_file.bcf

#Filter bcf from the previous step to exclude monomorphic SNPs and keep only variants with a quality score >19.
echo "- Further filtering results for a final vcf file..."
bcftools filter -e 'AC==0 || AC==AN' $output_file.bcf | bcftools filter --include 'QUAL > 19' > $output_file.vcf

echo "- Removing intermediate file $output_file.bcf ..."
if [[ -s $output_file.vcf ]]
	then
		rm $output_file.bcf
		echo "Removed intermediate file $output_file.bcf"

else

	echo "$output_file.vcf is empty. Something's fishy..." 
fi

