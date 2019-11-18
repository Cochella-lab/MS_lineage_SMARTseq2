###################
## bam2fastq script
## CBE optimised
###################
DIR_INPUT="${PWD}/ngs_raw/BAMs"
DIR_OUT="${PWD}/ngs_raw/FASTQs"
echo ${DIR_OUT};

mkdir -p ${DIR_OUT}
mkdir -p $PWD/logs/bam2fastq/

for file in $DIR_INPUT/*.bam;
do
    echo "$file"
    FILENAME="$(basename $file)";
    fname=${FILENAME%.bam};
    echo $fname

    ## creat the script for each sample
    script=$PWD/logs/${fname}_bam2fastq.sh
    cat <<EOF > $script
#!/usr/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=10
#SBATCH --mem=10000
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --output="$PWD/logs/bam2fastq/slurm-%j.out"
#SBATCH --error="$PWD/logs/bam2fastq/slurm-%j.err"
#SBATCH --job-name bam2fq
module load bedtools/2.27.1-foss-2018b;
bamToFastq -i $file -fq ${DIR_OUT}/${fname}.fastq;
EOF

    cat $script;
    sbatch $script
    #break;

done
