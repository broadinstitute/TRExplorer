SCRIPT=hprc-lps/convert_multisample_LPS_table_to_allele_frequency_histograms.py
META=hprc-lps_2026-05-19/1kGP_metadata.tsv
LPS_TABLE=hprc-lps_2026-05-19/hprc-lps.txt.gz

python3 $SCRIPT --sample-metadata-tsv $META --input-table $LPS_TABLE &

for s in male female; do
    python3 $SCRIPT --sample-metadata-tsv $META --input-table $LPS_TABLE --sex $s  &
done

for p in AFR  AMR  EAS  EUR  SAS; do
    python3 $SCRIPT --sample-metadata-tsv $META --input-table $LPS_TABLE --population $p &
    for s in male female; do
	python3 $SCRIPT --sample-metadata-tsv $META --input-table $LPS_TABLE --population $p --sex $s &
    done
done

wait
