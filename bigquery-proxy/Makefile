.PHONY: deploy logs test undeploy help load

# This script deploys a Google Cloud Function that proxies requests between the client browser and the BigQuery API. 

help:
	@echo "Tandem Repeat Explorer database proxy deployment Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  deploy    - Deploy or update the Google Cloud Function"
	@echo "  logs      - View recent logs from the deployed function"
	@echo "  test      - Run test queries against the deployed function"
	@echo "  undeploy  - Remove the deployed function and its logs"
	@echo "  load      - Load data into BigQuery by running load_bigquery_main_table.py"
	@echo "  help      - Show this help message (default)"

# Deploy the Cloud Function (allowed memory settings: 256Mb, 512M, 1024Mi, 2G, 4Gi, etc.)
# --service-account binds the function to a least-privilege SA (READER on
# cmg-analysis.tandem_repeat_explorer only, jobUser project-wide, objectAdmin on
# gs://tandem-repeat-explorer-export). This is a defence-in-depth mitigation for the SQL-gate
# weakness in main.py: even if the regex is bypassed, the SA cannot read other datasets in the
# project.
#
# For V4 signed export URLs the runtime SA must also have
# roles/iam.serviceAccountTokenCreator on ITSELF (so it can call IAM signBlob):
#   gcloud iam service-accounts add-iam-policy-binding \
#       trexplorer-proxy-sa@cmg-analysis.iam.gserviceaccount.com \
#       --member=serviceAccount:trexplorer-proxy-sa@cmg-analysis.iam.gserviceaccount.com \
#       --role=roles/iam.serviceAccountTokenCreator
deploy:
	gcloud functions deploy query_db \
		--runtime python310 \
		--trigger-http \
		--allow-unauthenticated \
		--service-account=trexplorer-proxy-sa@cmg-analysis.iam.gserviceaccount.com \
		--entry-point query_db \
		--region=us-central1 \
		--memory 256MB \
		--timeout=540s \
		--gen2

# View recent logs
logs:
	gcloud functions logs read query_db --region=us-central1 --limit=100

# Alias for logs target
log: logs

# Run test queries
test:
	@echo "Testing query: SELECT COUNT(*) AS c FROM cmg-analysis.tandem_repeat_explorer.catalog WHERE GencodeGeneName = 'HTT'"
	curl --silent --show-error --fail -X POST "https://us-central1-cmg-analysis.cloudfunctions.net/query_db" \
		-H "Content-Type: application/json" \
		-d '{ "sql": "SELECT COUNT(*) AS c FROM cmg-analysis.tandem_repeat_explorer.catalog WHERE GencodeGeneName = '\''HTT'\''" }' | jq '.'
	@echo "\nTesting query: SELECT chrom, start_0based, end_1based, GencodeGeneName FROM cmg-analysis.tandem_repeat_explorer.catalog WHERE GencodeGeneName = 'HTT' AND chrom = '4' AND start_0based > 3123000 AND end_1based < 3124000 LIMIT 100"
	@curl --silent --show-error --fail -X POST "https://us-central1-cmg-analysis.cloudfunctions.net/query_db" \
		-H "Content-Type: application/json" \
		-d '{ "sql": "SELECT chrom, start_0based, end_1based, GencodeGeneName FROM cmg-analysis.tandem_repeat_explorer.catalog WHERE GencodeGeneName = '\''HTT'\'' AND chrom = '\''4'\'' AND start_0based > 3123000 AND end_1based < 3124000 LIMIT 100" }' | jq '.'
	@echo "\nTesting query: SELECT chrom, start_0based, end_1based, GencodeGeneName FROM cmg-analysis.tandem_repeat_explorer.catalog WHERE GencodeGeneName = 'HTT'"
	curl --silent --show-error --fail -X POST "https://us-central1-cmg-analysis.cloudfunctions.net/query_db" \
		-H "Content-Type: application/json" \
		-d '{ "start_index": 0, "page_size": 3, "sql": "SELECT chrom, start_0based, end_1based, GencodeGeneName FROM cmg-analysis.tandem_repeat_explorer.catalog WHERE GencodeGeneName = '\''HTT'\''" }' | jq '.'

	curl --silent --show-error --fail -X POST "https://us-central1-cmg-analysis.cloudfunctions.net/query_db" \
		-H "Content-Type: application/json" \
		-d '{ "start_index": 1, "page_size": 3, "sql": "SELECT chrom, start_0based, end_1based, GencodeGeneName FROM cmg-analysis.tandem_repeat_explorer.catalog WHERE GencodeGeneName = '\''HTT'\''" }' | jq '.'

# Undeploy the function and clean up logs
undeploy:
	gcloud functions delete query_db --region=us-central1 

# Load the annotated TR catalog into BigQuery
load:
	python3 load_bigquery_main_table.py \
	    --catalog-path ~/code/tandem-repeat-catalogs/results__2026-02-01/release_draft_2026-02-01/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz \
	    --reference-fasta ~/hg38.fa \
	    --tenk10k-tsv ../data-prep/tenk10k_str_mt_rows.reformatted.tsv.gz \
	    --hprc256-tsv ../data-prep/hprc-lps/hprc-lps.per_locus_and_motif.256_samples.tsv.gz \
	    --aou1027-tsv ../data-prep/AoULR_phase1_TRGT_Weisburd_v1_combined.txt.gz \
	    --vamos-tsv ../data-prep/vamos_ori_and_eff_motif_columns.tsv.gz \
	    --hprc256-purity-stratified-tsv ../data-prep/hprc-lps/trgt-hprc.allele_size_purity.stratified.by_population.by_sex.256_samples.tsv.gz \
	    --hprc256-methylation-stratified-tsv ../data-prep/hprc-lps/trgt-hprc.methylation.stratified.by_population.by_sex.256_samples.tsv.gz \
	    --repeat-masker-lookup-json ../data-prep/hg38.RepeatMasker.lookup.json.gz \
	    --non-coding-annotations-bed ../data-prep/ncAnnot.v0.14.jul2024.filtered.bed.gz \
	    --tanudisastro-2024-lookup-json ../data-prep/literature_review/Tanudisastro_2024/Tanudisastro_2024_lookup.json.gz \
	    --manigbas-2024-lookup-json ../data-prep/literature_review/Manigbas_2024/Manigbas_2024_lookup.json.gz \
	    --reason ""
