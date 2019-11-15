#!/usr/bin/bash

BATCH_SIZE=6 #к-во файлов в батче
BATCH_NAME="${PWD}/velocyto_batches/batch_" #имя папки-батча

DIR_INPUT="${PWD}/results_v2/STAR"
DIR_OUT="${PWD}/LOOMS"
GTF="/groups/cochella/jiwang/annotations/Caenorhabditis_elegans.WBcel235.88.gtf"

mkdir -p $PWD/velocyto_batches
mkdir -p ${DIR_OUT}
mkdir -p $PWD/logs

batchIndex=0
i=0
mkdir "$BATCH_NAME$batchIndex"
for file in $DIR_INPUT/*.bam; # идем по всем .bam файлам в папке
do
    if (( $i == $BATCH_SIZE )) # если превысили лимит батча, то создаем новую папку
    then
        ((i=0))
        ((batchIndex++))
        mkdir "$BATCH_NAME$batchIndex"
    fi

    LINKNAME="$BATCH_NAME$batchIndex/$(basename $file)"
    ln "$file" "$LINKNAME" # копируем файл в батч
    ((i++))
done

for dir in ${BATCH_NAME%batch_}*; # идем по всем батчам
do
    echo "$dir"
    DIRNAME="$(basename $dir)";
    echo $DIRNAME
    echo "$DIR_OUT${PWD##*/}$DIRNAME.loom"
    ## creat the script for each sample
    script=$PWD/logs/${DIRNAME}.sh
    cat <<EOF > $script
#!/usr/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=60
#SBATCH --mem=10000
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $PWD/logs/$fname.out
#SBATCH -e $PWD/logs/$fname.err
#SBATCH --job-name velocyto

module load velocyto/0.17.17-foss-2018b-python-3.6.6;
velocyto run-smartseq2 -e $DIR_OUT/${PWD##*/}_$DIRNAME $dir/*.bam $GTF

EOF

    cat $script;
    sbatch $script # так эти скрипты запускаются

done
