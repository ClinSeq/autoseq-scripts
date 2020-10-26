FILES=( alasccaCNA.R \
        bed_to_regions.py \
        cnvgtf2bed.py \
        filter_erroneus_alt.py \
        gatk-klevebring \
        gtf_filter_unused_chrs.py \
        intersect-msi-sites.sh \
        msisensor \
        mutect2 \
        picard_interval_list_to_bed6_converter.py \
        qdnaseq.R \
        qdnaseq2bed.py \
        slopIntervalList.py \
        target_coverage_histogram.py \
        vcf_add_sample.py \
        vcfsorter.pl \
        )

JARS=(GenomeAnalysisTK-3.5.jar \
      GenomeAnalysisTK-Klevebring.jar
      )

mkdir -p $PREFIX/bin

for J in ${JARS[@]}; do
  cp $J $PREFIX/bin/$J
done


for F in ${FILES[@]}; do
  cp $F $PREFIX/bin/$F
done
