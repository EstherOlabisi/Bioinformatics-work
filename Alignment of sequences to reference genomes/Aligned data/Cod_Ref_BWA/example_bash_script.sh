#!/usr/bin/bash

## Usage (within an srun session - could also add SLURM #SBATCH header stuff and run through the queue
## ./example_bash_script.sh

module load samtools

echo "ind raw assembled" > b_bwa_assembled_per_ind.txt

for file in Burbot_Ref_BWA/*sorted.bam
do 
indname=`echo $file | sed 's/.*\///g' | sed 's/\.sorted\.bam//g'`
raw=`samtools stats $file | grep "raw total sequences:" | sed 's/SN\t.*:\t//g' | sed 's/\t# excluding supplementary and secondary reads/ /g'`
assembled=`samtools stats $file | grep "reads mapped:" | sed 's/SN\t.*:\t//g'`
echo "$indname $raw $assembled" >> b_bwa_assembled_per_ind.txt
done
