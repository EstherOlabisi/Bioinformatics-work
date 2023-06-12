#!/bin/sh

## This script uses bcftools mpileup to call variants from sorted bam files 
## usage (for testing with just one individual):
## sbatch run_bcf_caller.sh file.sorted.bam

#SBATCH --account=def-lukens
#SBATCH --time=0-01:00:00 ## days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # number of threads
#SBATCH --mem=16000 # requested memory (in MB)
#SBATCH --mail-type=END

module load samtools
module load picard
module load htslib nixpkgs gcc freebayes

#assign variables
bamfile=$1
rg_bamfile=`echo $bamfile | sed 's/\.sorted/_rg\.sorted/g'`	#new names for files with read groups
output_file="freebayes_result"

echo "- Adding Read Groups to all bam files..."
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$bamfile O=$rg_bamfile RGID=$bamfile RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$bamfile

echo "- Reindexing new bam files containing read groups..."
samtools index $rg_bamfile

#Apply the same filters as in bcftools mpileup
echo "- Calling variants for all bam files..."
freebayes -f burbot_2021.fasta -b *_rg.sorted.bam | vcffilter -f 'TYPE = snp' -f 'QUAL > 19' | vcffilter -f "! ( AC = 0 )" | vcffilter -f "! ( AC = AN )" > $output_file.vcf

if [[ -s $output_file.vcf ]]
	then
		echo "$output_file.vcf is ready."
fi

