#!/bin/sh

## This script uses bowtie2 to map reads (.fastq) to reference genome
## usage (for testing with just one individual):
## sbatch run_bowtie_queuesub.sh EGM16_0001.fastq

#SBATCH --account=def-lukens
#SBATCH --time=0-00:12:00 ## days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # number of threads
#SBATCH --mem=16000 # requested memory (in MB)
#SBATCH --mail-type=END

module load bowtie2/2.4.4
module load samtools/1.16.1

fastq=$1
basename=`echo $fastq | sed 's/\.f.*\.gz.*//'`

echo "Starting alignment of $fastq to reference genome"
bowtie2 -x ../burbot_reference_genome/GCA_900302385.1_ASM90030238v1_genomic.fna -U $fastq -S bowtie2_assem/$basename.sam --very-sensitive-local 
echo "Converting sam to bam for $basename"
samtools view -b -S -o bowtie2_assem/$basename.bam bowtie2_assem/$basename.sam

echo "Sorting and indexing bam files for $basename"
samtools sort bowtie2_assem/$basename.bam -o bowtie2_assem/$basename.sorted.bam
samtools index bowtie2_assem/$basename.sorted.bam

echo "Cleaning up the mess... just a minute!"
if [[ -s bowtie2_assem/$basename.bam ]]
   then
       rm bowtie2_assem/$basename.sam
       echo "removed $basename.sam"

else
    echo "$basename.sam is empty! Something's fishy..."
fi



if [[ -s bowtie2_assem/$basename.bam ]]
   then
       rm bowtie2_assem/$basename.bam
       echo "removed $basename.bam"

else
    echo "$basename.bam is empty! Something's fishy..."
fi
