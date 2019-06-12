#!/bin/bash
#SBATCH --job-name MapB
#SBATCH --array=1-12
#SBATCH --output temp/MapB.out
#SBATCH --error temp/MapB.err
#SBATCH --mail-type FAIL
#SBATCH --qos=fast
#SBATCH --cpus-per-task 3
#SBATCH --mem=30GB
#SBATCH -p common

### Should be run on SLURM cluster using :
# sbatch RunRIPSeq.sh listRIPSeq.txt
# with text file containing list of dataset names


folder=/pasteur/projets/policy01/BioIT/8001_Pagliuso_RNABinding/
dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)

# prepare sbatch variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# prepare variables
genome=NC_003210
genomeFolder=$folder/Genome/
annotation=$folder/Genome/NC_003210
mappingFolder=$folder/Mapping/bowtie_results/
readFile=$folder/FastQ/${dataName}.fastq
wigfile=$folder/Mapping/${dataName}.bw
bamFileUnfiltered=$folder/Mapping/bowtie_results/${dataName}
bamFileFiltered=$folder/Mapping/${dataName}
output_fcount=${folder}/Expression/FCount_${dataName}

if [ ! -d $mappingFolder ]; then
	mkdir $mappingFolder
fi

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
module add AlienTrimmer/0.4.0
AlienTrimmer -i $scratch/RNASeq_raw/${dataName}.fastq -c $folder/Genome/alienTrimmerPF8contaminants.fasta -o $folder/FastQ/${dataName}.fastq

module add bowtie2
# create genome
#bowtie2-build $genomeFolder/$genome.fna $genomeFolder/$genome
# map datasets
bowtie2 -p 3 --very-sensitive -x $genomeFolder/$genome -U $readFile -S $bamFileUnfiltered.sam 2>$bamFileUnfiltered.log

# # # filter bam files by removing data with score = 0
module add samtools
samtools view -bS ${bamFileUnfiltered}.sam > ${bamFileUnfiltered}.bam
samtools view -b -q 1 ${bamFileUnfiltered}.bam > ${bamFileUnfiltered}_r.bam
samtools sort ${bamFileUnfiltered}_r.bam ${bamFileFiltered}
samtools index ${bamFileFiltered}.bam

# # # # create bigwig
source activate py27
bamCoverage --normalizeUsing RPKM -b ${bamFileFiltered}.bam -o $wigfile

# # FeatureCount for RIPSeq
module add subread
featureCounts -t "exon" -o ${output_fcount}.txt -a ${annotation}_all.gff ${bamFileFiltered}.bam

echo \"Everything is finished\"
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
