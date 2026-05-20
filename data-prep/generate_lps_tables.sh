SCRIPT=hprc-lps/convert_multisample_LPS_table_to_allele_frequency_histograms.py
EXTRACT_SCRIPT=hprc-lps/extract_vcf_interval_metadata.py
META=hprc-lps_2026-05-19/1kGP_metadata.tsv
LPS_TABLE=hprc-lps_2026-05-19/hprc-lps.txt.gz
VCF=hprc-lps_2026-05-19/trgt-hprc.vcf.gz
INTERVAL_TSV=hprc-lps_2026-05-19/trgt-hprc.interval_metadata.tsv.gz

python3 $EXTRACT_SCRIPT --input-vcf $VCF --output-tsv $INTERVAL_TSV

python3 $SCRIPT --sample-metadata-tsv $META --input-table $LPS_TABLE --vcf-interval-tsv $INTERVAL_TSV &

for s in male female; do
    python3 $SCRIPT --sample-metadata-tsv $META --input-table $LPS_TABLE --vcf-interval-tsv $INTERVAL_TSV --sex $s  &
done

for p in AFR  AMR  EAS  EUR  SAS; do
    python3 $SCRIPT --sample-metadata-tsv $META --input-table $LPS_TABLE --vcf-interval-tsv $INTERVAL_TSV --population $p &
    for s in male female; do
	python3 $SCRIPT --sample-metadata-tsv $META --input-table $LPS_TABLE --vcf-interval-tsv $INTERVAL_TSV --population $p --sex $s &
    done
done

wait
