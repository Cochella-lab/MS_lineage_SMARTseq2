###################
## Merging technical replicates
###################
DIR_INPUT="${PWD}/ngs_raw/unmerged_BAMs"
DIR_OUT="${PWD}/ngs_raw/BAMs"

mkdir -p ${DIR_OUT}
mkdir -p $PWD/logs

for file in $DIR_INPUT/CCYTEANXX_4*.bam;
do
    echo "$file"
    FILENAME="$(basename $file)";
    fname=${FILENAME%.bam};
    echo $fname

    INPUT1=$file
    INPUT2=$(echo $file | sed 's|CCYTEANXX_4#80194_|HHGHNBGX9_1#80194_|')
    INPUT3=$(echo $file | sed 's|CCYTEANXX_4#80194_|CD2GTANXX_5#80194_|')

    ## creat the script for each sample
    script=$PWD/logs/${fname}_merger.sh
    cat <<EOF > $script
#!/bin/sh
#SBATCH --cpus-per-task=1
#SBATCH --time=60
#SBATCH --mem=10000
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $PWD/logs/$fname.out
#SBATCH -e $PWD/logs/$fname.err
#SBATCH --job-name samtools_merge
module load samtools/0.1.20-foss-2018b;
samtools merge ${DIR_OUT}/${fname}.bam $INPUT1 $INPUT2 $INPUT3;
EOF

    cat $script;
    sbatch $script
    #break;
done
