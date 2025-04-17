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
	@echo "  load      - Load data into BigQuery by running load_data_into_bigquery.py"
	@echo "  help      - Show this help message (default)"

# Deploy the Cloud Function (allowed memory settings: 256Mb, 512M, 1024Mi, 2G, 4Gi, etc.)
deploy:
	gcloud functions deploy query_db \
		--runtime python310 \
		--trigger-http \
		--allow-unauthenticated \
		--entry-point query_db \
		--region=us-central1 \
		--memory 256MB \
		--gen2

# View recent logs
logs:
	gcloud functions logs read query_db --region=us-central1 --limit=10

# Alias for logs target
log: logs

# Run test queries
test:
	@echo "Testing query: SELECT COUNT(*) AS c FROM cmg-analysis.tandem_repeat_explorer.catalog WHERE GencodeGeneName = 'HTT'"
	curl --silent --show-error --fail -X POST "https://us-central1-cmg-analysis.cloudfunctions.net/query_db" \
		-H "Content-Type: application/json" \
		-d '{ "sql": "SELECT COUNT(*) AS c FROM cmg-analysis.tandem_repeat_explorer.catalog WHERE GencodeGeneName = '\''HTT'\''" }' | jq '.'
	@echo "\nTesting query: SELECT chrom, start_0based, end_1based, GencodeGeneName FROM cmg-analysis.tandem_repeat_explorer.catalog WHERE GencodeGeneName = 'HTT' AND chrom = 'chr4' AND start_0based > 3123000 AND end_1based < 3124000 LIMIT 100"
	@curl --silent --show-error --fail -X POST "https://us-central1-cmg-analysis.cloudfunctions.net/query_db" \
		-H "Content-Type: application/json" \
		-d '{ "sql": "SELECT chrom, start_0based, end_1based, GencodeGeneName FROM cmg-analysis.tandem_repeat_explorer.catalog WHERE GencodeGeneName = '\''HTT'\'' AND chrom = '\''chr4'\'' AND start_0based > 3123000 AND end_1based < 3124000 LIMIT 100" }' | jq '.'
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
	python3 load_data_into_bigquery.py
