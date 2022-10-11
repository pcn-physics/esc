#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script
#   for TACC Frontera CLX nodes
#
#   *** MPI Job in Normal Queue ***
# 
# Last revised: 20 May 2019
#
# Notes:
#
#   -- Launch this script by executing
#      "sbatch clx.mpi.slurm" on a Frontera login node.
#
#   -- Use ibrun to launch MPI codes on TACC systems.
#      Do NOT use mpirun or mpiexec.
#
#   -- Max recommended MPI ranks per CLX node: 56
#      (start small, increase gradually).
#
#   -- If you're running out of memory, try running
#      fewer tasks per node to give each task more memory.
#
#----------------------------------------------------

#SBATCH -J Job%j           # Job name
#SBATCH -o /scratch1/06800/pcn/data/stdfiles/outputs/out.o%j       # Name of stdout output file
#SBATCH -e /scratch1/06800/pcn/data/stdfiles/errors/err.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 4               # Total # of nodes 
#SBATCH -n 32              # Total # of mpi tasks
#SBATCH -t 12:30:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=pamela.c.nguyen@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#       # Project/Allocation name (req'd if you have more than 1)


# Any other commands must follow all #SBATCH directives...
module list
pwd
date

# Launch MPI code... 
ibrun ./mycode.exe         # Use ibrun instead of mpirun or mpiexec

#End of submission script

#Header of running script which includes modules to be loaded and defining source and input files
#!/bin/bash

#module load samtools

SCRATCH=/scratch1/06800/pcn
WORK1=/work2/06800/pcn/frontera/FASTQ
WORK2=/work2/06548/esarinay/frontera/FASTQ/SUPRESSOR_NEXTERA
HOME=/home1/06800/pcn

source /home1/06800/pcn/.bashrc

#End of header file... Begin Base script

#Defining output names and creation of respective directories

#output1=$(basename $input1| awk -F "_" '{print($4"_"$5)}')
#output2=$(basename $input2| awk -F "_" '{print($4"_"$5)}')


# VCF to Bed NEED TO FIXXXXX
#/home1/06800/pcn/bedops/bin/vcf2bed < /work/06800/pcn/stampede2/alignedVCF/ESC-HSL-56_unique_aligned.g.vcf > /scratch/06800/pcn/sortedBed/ESC-HSL-55-sorted.bed

#Creates an index from the celegan reference genome
bowtie2-build $WORK1/Caenorhabditis_elegans.WBcel235.dna.toplevel.\
fa $WORK1/WBcel235

# Converts from .fa to .sam given a reference genome
bowtie2 -N 1 -5 5 -3 5 -k 1 -x $WORK1/WBcel235 -1 $WORK2/HSL-214_R1.fastq.gz -2 $WORK2/HSL-214_R2.fastq.gz  -S $SCRATCH/data/HSL-214_aligned.sam

# Converts from .sam to .bam given a sam file
#samtools view -q 5 -b $SCRATCH/data/ESC-HSL-56_aligned.sam >$SCRATCH/data/ESC-HSL-56_unique_aligned.bam

# Sorts the given bam file
#samtools sort $SCRATCH/data/ESC-HSL-56_unique_aligned.bam -o $SCRATCH/data/ESC-HSL-56_unique_aligned_sorted.bam

# Indexes the genome sorted bam file
#samtools index $SCRATCH/data/ESC-HSL-56_unique_aligned_sorted.bam

# Assign Read Groups using GATK
# Refer to https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
# for more information
#java -jar /home1/06800/pcn/picard/picard.jar AddOrReplaceReadGroups \
#       I=/work/06548/esarinay/DATA/200807_ALIGNMENTS/ESC-HSL-69_aligned_sorted_dedup.bam \
#       O=/work/06800/pcn/stampede2/alignments/ESC-HSL-69_aligned_sorted_dedup_with_reads.bam \
#       RGID=69 \
#       RGLB=ESC-HSL \
#       RGPL=ILLUMINA \
#       RGPU=unit1 \
#       RGSM=1

# Indexes dedup read bam files with samtools
#/home1/06800/pcn/samtools/samtools-1.10/samtools index /work/06800/pcn/stampede2/alignments/ESC-HSL-69_aligned_sorted_dedup_with_reads.bam

# Runs GATK HaplotypeCaller (Version I)
# Input: dedup bam
# Output: vcf
#$HOME/gatk/gatk-4.1.8.1/gatk \
#      HaplotypeCaller \
#      -R /work/06800/pcn/stampede2/fasta/caenorhab\
#ditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa \
#      -I /work/06800/pcn/stampede2/alignments/ESC-HSL-69_aligned_sorted_dedup_with_reads\
#.bam \
#      -O /scratch/06800/pcn/alignedVCF/ESC-HSL-69_unique_aligned.g.vcf \
#      -contamination ${default=0 contamination} ${true="-ERC GVCF" false="" make_gvcf}

# Runs GATK HaplotypeCaller (Version II)
# Input: bai
# Output: vcf
# Note: This variation of the HaplotypeCaller was implemented as a requirement due to IGV
# The formatting in the GATK HaplotypeCaller Version I is easier to read and organise
# but the both versions are sufficient outputs the necessary vcf
#/home1/06800/pcn/gatk/gatk-4.1.8.1/gatk HaplotypeCaller -R /work/06800/pcn/stampede2/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa -I /work/06800/pcn/stampede2/esc/alignments/sample7c/sample7c.bai -O sample7c.g.vcf

