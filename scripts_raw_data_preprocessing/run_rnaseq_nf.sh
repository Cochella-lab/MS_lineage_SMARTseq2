##############
# the test script for rnaseq nf 
# https://github.com/lengfei5/nf-core-rnaseq/blob/master/docs/configuration/adding_your_own.md
#
##############
fqs="$PWD/ngs_raw/FASTQs/*.fastq"
output="$PWD/results_v2"
workdir="$PWD/work"

star_index="/groups/cochella/jiwang/Genomes/C_elegans/WBcel235/index_4star"
GTF="/groups/cochella/jiwang/annotations/Caenorhabditis_elegans.WBcel235.88.gtf"

nextflow run /groups/cochella/jiwang/scripts/rnaseq_nf -profile slurm \
--reads "$fqs" \
--singleEnd \
--unstranded \
--star_index "$star_index" \
--gtf "$GTF" \
--max_IntronL 50000 \
--outdir "$output" \
-w "$workdir" \
--skip_preseq \
--saveAlignedIntermediates \
--skip_genebody_coverage \
--skip_dupradar \
--skip_edger \
--sampleLevel \
--skip_multiqc -resume
