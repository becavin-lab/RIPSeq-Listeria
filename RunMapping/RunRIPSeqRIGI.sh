#!/bin/bash
#SBATCH --job-name gsnap
#SBATCH --array=1-12
#SBATCH --output temp/gsnap.out
#SBATCH --error temp/gsnap.err
#SBATCH --mail-type FAIL
#SBATCH --qos=hubbioit
#SBATCH --cpus-per-task 3
#SBATCH --mem=30GB
#SBATCH -p hubbioit

### Should be run on SLURM cluster using :
# sbatch RunRIPSeqRIGI.sh listRIPSeqRIGI.txt
# with text file containing list of dataset names


folder=/pasteur/projets/policy01/BioIT/8001_Pagliuso_RNABinding/MikaelRNASeq/
dataName=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n 1)

# prepare sbatch variables
if [ ! -d $folder/temp/${SLURM_JOB_NAME} ]; then
	mkdir $folder/temp/${SLURM_JOB_NAME}
fi
shFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.sh
logFile=$folder/temp/${SLURM_JOB_NAME}/${dataName}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log

# prepare variables
genome=NC_003210.fna
genomeFolder=$folder/../Genome/
annotation=$genomeFolder/NC_003210
readFile=$folder/FastQ/${dataName}.fastq
wigfile=$folder/Mapping/${dataName}.bw
bamFileUnfiltered=$folder/Mapping/GSNAP_results/${dataName}
bamFileFiltered=$folder/Mapping/${dataName}_SNAP
output=${folder}/Expression/HTSeq_${dataName}

# Create script SH to run in qsub
scriptSH="""
#!/bin/bash -l
unpigz -p 3 -k $folder/RNASeq_raw/all_gz/$dataName.fastq.gz
mv $folder/RNASeq_raw/all_gz/$dataName.fastq $folder/RNASeq_raw/$dataName.fastq

module add AlienTrimmer/0.4.0
AlienTrimmer -i $folder/RNASeq_raw/${dataName}.fastq -c $folder/../Genome/alienTrimmerPF8contaminants.fasta -o $folder/FastQ/${dataName}.fastq

# Run GSNAP
gsnap -d NC_003210 -D $genomeFolder/GMAP/NC_003210/ -m 4 -A sam -o $folder/Mapping/GSNAP_results/$dataName.sam $folder/FastQ/$dataName.fastq
# # # # filter bam files by removing data with score = 0
module add samtools
samtools view -Sb -q 1 ${bamFileUnfiltered}.sam > ${bamFileUnfiltered}.bam
samtools sort ${bamFileUnfiltered}.bam ${bamFileFiltered}
samtools index ${bamFileFiltered}.bam

# # # # create bigwig
source activate py27
bamCoverage --normalizeUsing RPKM -b ${bamFileFiltered}.bam -o $wigfile

module add HTSeq/0.9.1
# # count reads per genes with all elements
htseq-count -t exon -i gene_id -s no -m union --nonunique all -f 'bam' ${folder}/Mapping/${dataName}_SNAP.bam ${annotation}_all.gff > ${output}.txt

# echo \"Everything is finished\"
"""

# save script and run with qsub
echo """$scriptSH""" > ${shFile}
#echo $scriptSH

# Run script
srun -c ${SLURM_CPUS_PER_TASK} -o ${logFile} bash ${shFile}
