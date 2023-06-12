#!/bin/sh

## This script converts sam files to bam files and then sorts and indexes the bam files.
## usage (for testing with just one individual):
## sbatch run_sam_to_bam.sh samfile.sam

#SBATCH --account=def-lukens
#SBATCH --time=0-00:5:00 ## days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # number of threads
#SBATCH --mem=16000 # requested memory (in MB)
#SBATCH --mail-type=END

#load module
module load samtools/1.16.1

#assign variables
samfile=$1
basename=`echo $samfile | sed 's/\.gz.sam//'`

#convert from  sam to bam
echo "Converting sam to bam for $basename"
samtools view -b -S -o bam_assem/$basename.bam $samfile

#sort and index bam files
echo "Sorting and indexing bam files for $basename"
samtools sort bam_assem/$basename.bam -o bam_assem/$basename.sorted.bam
samtools index bam_assem/$basename.sorted.bam

#remove intermediate unsorted bam files
if [[ -s bam_assem/$basename.bam ]]
	then
		rm bam_assem/$basename.bam
		echo "removed $basename.bam"
else
	echo "$basename.bam is empty! Something's fishy..."
fi
