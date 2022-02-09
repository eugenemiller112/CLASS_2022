#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=Hepg2_Pol_test
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=john.rinn@colorado.edu
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=10:00:00

pwd; hostname; date
echo "Lets do chipseq"

module load singularity/3.1.1

nextflow run nf-core/chipseq -r 1.2.1 \
-profile singularity \
--single_end \
--input design.csv \
--fasta /Shares/rinn_class/data/genomes/human/gencode/v32/GRCh38.p13.genome.fa \
--gtf /Shares/rinn_class/data/genomes/human/gencode/v32/gencode.v32.annotation.gtf \
--macs_gsize 3.2e9 \
--blacklist /scratch/Shares/rinnclass/data/hg38-blacklist.v2.bed \
--email john.rinn@colorado.edu \
-resume \
-c nextflow.config

date
