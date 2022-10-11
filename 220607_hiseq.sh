#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script
#   for TACC Frontera CLX nodes
#
#   *** Serial Job in Small Queue***
# 
# Last revised: 22 June 2021
#
# Notes:
#
#  -- Copy/edit this script as desired.  Launch by executing
#     "sbatch clx.serial.slurm" on a Frontera login node.
#
#  -- Serial codes run on a single node (upper case N = 1).
#       A serial code ignores the value of lower case n,
#       but slurm needs a plausible value to schedule the job.
#
#  -- Use TACC's launcher utility to run multiple serial 
#       executables at the same time, execute "module load launcher" 
#       followed by "module help launcher".
#----------------------------------------------------

#SBATCH -J 220607_hiseq           # Job name
#SBATCH -o /scratch/06800/pcn/data/stdfiles/outputs/220607_hiseq.o%j       # Name of stdout output file
#SBATCH -e /scratch/06800/pcn/data/stdfiles/outputs/220607_hiseq.e%j       # Name of stderr error file
#SBATCH -p small           # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A OTH21143 #Project/Allocation name
#SBATCH --mail-user=pamela.c.nguyen@utexas.edu

# Any other commands must follow all #SBATCH directives...
module list
pwd
date

# Launch serial code...
./myprogram         # Do not use ibrun or any other MPI launcher

GENOME=/work2/06800/pcn/frontera/FASTQ/genome/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa
REFERENCE=/work2/06800/pcn/frontera/FASTQ/reference/WBcel235
HOME=/home1/06800/pcn
SCRATCH=/scratch1/06800/pcn/data
WORK1=/work2/06800/pcn/frontera/FASTQ/July2020_sequencing_run-selected
NUM=56

#WORK2=/work2/06548/esarinay/frontera/FASTQ/SUPRESSOR_NEXTERA


source /home1/06800/pcn/.bashrc

#exec&>> $SCRATCH/logs/ESC-HSL-$NUM.log

module load fastqc #load bio container via tacc

start=`date +%s`

#Creates an index from the celegan reference genome
#NOTE: reference file is .dna_sm NOT .dna or .dna_rm
#NOTE: we only have to run this command once for the rest of the analysis
#bowtie2-build -f\
#       $GENOME \
#       $REFERENCE

# bowtie2-build -f\
#       /home/pcn/utexas/ESC/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa \
#       /home/pcn/utexas/ESC/dna/reference/WBcel235

# Converts from .fastq to .sam given a reference genome
bowtie2 -N 1 -5 5 -3 5 -k 1 -x \
      $REFERENCE -1 \
      $WORK1/ESC-HSL-56_R1.fastq.gz -2 \
      $WORK1/ESC-HSL-56_R2.fastq.gz  -S \
      $SCRATCH/alignment/ESC-HSL-$NUM-aligned.sam

# Converts from .fastq to .sam given a reference genome
# bowtie2 -N 1 -5 5 -3 5 -k 1 -x \
#       /home/pcn/utexas/ESC/dna/reference/WBcel235 -1 \
#       /home/pcn/utexas/ESC/NOV2021_DATA/ESC-HSL-351_R1.fastq.gz -2 \
#       /home/pcn/utexas/ESC/NOV2021_DATA/ESC-HSL-351_R2.fastq.gz  -S \
#       /home/pcn/utexas/ESC/scratch/data/alignment/ESC-HSL-351-aligned.sam


# Converts from .sam to .bam given a sam file
samtools view -q 5 -b \
      $SCRATCH/alignment/ESC-HSL-$NUM-aligned.sam \
      >$SCRATCH/alignment/ESC-HSL-$NUM-unique-aligned.bam

# Sorts the given bam file
samtools sort \
      $SCRATCH/alignment/ESC-HSL-$NUM-unique-aligned.bam -o \
      $SCRATCH/alignment/ESC-HSL-$NUM-unique-aligned-sorted.bam

# Indexes the genome sorted bam file
samtools index \
      $SCRATCH/alignment/ESC-HSL-$NUM-unique-aligned-sorted.bam

#######################################
# Assign Read Groups using GATK
# Refer to https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
# https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-
# for more information
# NOTE: dedup bams were not used
#java -jar $HOME/picard.jar AddOrReplaceReadGroups \
#       I=$SCRATCH/alignment/ESC-HSL-$NUM-unique-aligned-sorted.bam \
#       O=$SCRATCH/alignment/ESC-HSL-$NUM-aligned-sorted-with-reads.bam \
#       RGID=$NUM \
#       RGLB=HSL \
#       RGPL=ILLUMINA \
#       RGPU=unit1 \
#       RGSM=1

#######################################

# java -jar /home/pcn/utexas/ESC/picard.jar AddOrReplaceReadGroups \
#        I=/home/pcn/utexas/ESC/scratch/data/alignment/ESC-HSL-349-unique-aligned-sorted.bam \
#        O=/home/pcn/utexas/ESC/scratch/data/alignment/ESC-HSL-349-aligned-sorted-with-reads.bam \
#        RGID=349 \
#        RGLB=HSL \
#        RGPL=ILLUMINA \
#        RGPU=unit1 \
#        RGSM=1

#######################################
# Indexes read bam files with samtools
#samtools index \
#$SCRATCH/alignment/ESC-HSL-$NUM-aligned-sorted-with-reads.bam
# samtools index \
# /home/pcn/utexas/ESC/scratch/data/alignment/ESC-HSL-349-aligned-sorted-with-reads.bam
#######################################


# # Create reference.dict
# #INPUT: ref.fa
# #OUTPUT: ref.dict
# # https://gatk.broadinstitute.org/hc/en-us/articles/360036729911-CreateSequenceDictionary-Picard-
# #NOTE: we only have to run this command once for the rest of the analysis
# java -jar /home/pcn/utexas/ESC/picard.jar CreateSequenceDictionary \
# R= /home/pcn/utexas/ESC/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa \
# O= /home/pcn/utexas/ESC/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.dict

# #Create the fasta index file:
# #INPUT: ref.fa
# #OUTPUT: ref.fa.fai
# #NOTE: we only have to run this command once for the rest of the analysis
# samtools faidx /home/pcn/utexas/ESC/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa

#######################################
# Runs GATK HaplotypeCaller (Version I)
# refer to https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
# Input: .bam
# Output: .vcf
#      gatk\
#      HaplotypeCaller \
#      -R $GENOME \
#      -I $SCRATCH/alignment/ESC-HSL-$NUM-aligned-sorted-with-reads.bam \
#      -O $SCRATCH/vcf/ESC-HSL-$NUM.g.vcf 

      # /home/pcn/utexas/ESC/gatk-4.2.2.0/./gatk\
      # HaplotypeCaller \
      # -R /home/pcn/utexas/ESC/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa \
      # -I /home/pcn/utexas/ESC/scratch/data/alignment/ESC-HSL-349-aligned-sorted-with-reads.bam \
      # -O /home/pcn/utexas/ESC/scratch/data/vcf/ESC-HSL-349.g.vcf 
#######################################

# Convert from VCF to BED
# Input: .vcf
# Output: .bed
#vcf2bed \
#      < $SCRATCH/vcf/ESC-HSL-$NUM.g.vcf \
#      > $SCRATCH/sorted_bed/ESC-HSL-$NUM-sorted.bed

# vcf2bed \
#       < /home/pcn/utexas/ESC/scratch/data/vcf/ESC-HSL-349.g.vcf \
#       > /home/pcn/utexas/ESC/scratch/data/sorted_bed/ESC-HSL-349-sorted.bed

end=`date +%s`

runtime=$((end-start))
