###################
## bam2fastq script
###################
DIR_INPUT="${PWD}/ngs_raw/BAMs"
DIR_OUT="${PWD}/ngs_raw/FASTQs"
echo ${DIR_OUT};

mkdir -p ${DIR_OUT}
mkdir -p $PWD/logs

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
#SBATCH --time=60 
#SBATCH --mem=10000 
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $PWD/logs/$fname.out 
#SBATCH -e $PWD/logs/$fname.err 
#SBATCH --job-name bam2fq
module load bedtools/2.25.0-foss-2017a;
bamToFastq -i $file -fq ${DIR_OUT}/${fname}.fastq;
EOF
    
    cat $script;
    sbatch $script
    #break;
    
done