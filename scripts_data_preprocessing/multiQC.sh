###################
## MultiQC script
###################
wdir="${PWD}"
#!/usr/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=60
#SBATCH --mem=128G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $PWD/logs/$fname.out
#SBATCH -e $PWD/logs/$fname.err
#SBATCH --job-name "MultiQc"

module load multiqc/1.3-foss-2017a-python-2.7.13;
multiqc "$wdir";
