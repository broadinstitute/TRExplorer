SCRIPT=hprc-lps/convert_multisample_LPS_table_to_allele_frequency_histograms.py
META=hprc-lps_2026-05-19/1kGP_metadata.tsv
LPS_TABLE=hprc-lps_2026-05-19/hprc-lps.txt.gz

population=AFR
for n in 10 20 30 40 50 60 70 80 90; do
    python3 $SCRIPT --sample-metadata-tsv $META --input-table $LPS_TABLE --population $population --num-samples $n &
done

wait
