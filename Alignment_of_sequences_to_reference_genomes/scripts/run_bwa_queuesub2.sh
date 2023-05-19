#!/bin/sh

## This script uses bwa to map reads (.fastq) to reference genome
## usage (for testing with just one individual):
## sbatch run_bwa_queuesub.sh EGM16_0001.fastq

#SBATCH --account=def-lukens
#SBATCH --time=0-00:12:00 ## days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # number of threads
#SBATCH --mem=16000 # requested memory (in MB)
#SBATCH --mail-type=END

module load bwa/0.7.17
module load samtools/1.16.1

fastq=$1
basename=`echo $fastq | sed 's/\.f.*\.gz.*//'`

echo "Starting alignment of $fastq to reference genome"
bwa mem -t 16 ../cod_reference_genome/GCF_902167405.1_gadMor3.0_genomic.fna $fastq >  bwa_assem/$basename.sam
echo "Converting sam to bam for $basename"
samtools view -b -S -o bwa_assem/$basename.bam bwa_assem/$basename.sam

echo "Sorting and indexing bam files for $basename"
samtools sort bwa_assem/$basename.bam -o bwa_assem/$basename.sorted.bam
samtools index bwa_assem/$basename.sorted.bam

echo "Cleaning up the mess... just a minute!"
if [[ -s bwa_assem/$basename.bam ]]
   then
       rm bwa_assem/$basename.sam
       echo "removed $basename.sam"

else
    echo "$basename.sam is empty! Something's fishy..."
fi



if [[ -s bwa_assem/$basename.bam ]]
   then
       rm bwa_assem/$basename.bam
       echo "removed $basename.bam"

else
    echo "$basename.bam is empty! Something's fishy..."
fi
