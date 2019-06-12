#!/bin/bash
#SBATCH --job-name Map
#SBATCH --array=1-6
#SBATCH --output temp/Map.out
#SBATCH --error temp/Map.err
#SBATCH --mail-type FAIL
#SBATCH --qos=hubbioit
#SBATCH --cpus-per-task 3
#SBATCH --mem=50GB
#SBATCH -p hubbioit


### Should be run on SLURM cluster using :
# sbatch RunRNASeq.sh listRNASeq.txt
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
output_htseq=${folder}/Expression/HTSeq_${dataName}

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
unpigz -p 3 -k $folder/RNASeq_raw/all_gz/$dataName.fastq.gz
mv $folder/RNASeq_raw/all_gz/$dataName.fastq $folder/RNASeq_raw/$dataName.fastq

module add AlienTrimmer/0.4.0
## cutoff at length 45 for RNASeq
AlienTrimmer -l 45 -i $folder/RNASeq_raw/${dataName}.fastq -c $folder/Genome/alienTrimmerPF8contaminants.fasta -o $folder/FastQ/${dataName}.fastq

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

module add HTSeq/0.9.1
# count reads per genes
htseq-count -t exon -i gene_id -a 10 -s no --nonunique all -f 'bam' ${folder}/Mapping/${dataName}.bam ${annotation}.gff > ${output_htseq}.txt

# echo \"Everything is finished\"
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
