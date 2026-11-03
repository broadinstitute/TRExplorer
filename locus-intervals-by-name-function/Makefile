.PHONY: deploy logs test unittest help

# This Makefile deploys the Google Cloud Function that the website's gene and transcript name searches use to look up
# where the loci with a name are (see main.py). Its lookup files are written by
# ../data-prep/locus_intervals_by_name/write_locus_intervals_by_name_lookup_file.py.

FUNCTION_NAME = lookup_locus_intervals_by_name
FUNCTION_URL = https://us-central1-cmg-analysis.cloudfunctions.net/$(FUNCTION_NAME)
# The catalog table that the website points to (TABLE_ID in ../website/header_template.html), for the test target
CATALOG_TABLE_ID = $(shell sed -n -E "s/.*const TABLE_ID = '(catalog_[0-9_]+)'.*/\1/p" ../website/header_template.html)

help:
	@echo "Name lookup Cloud Function deployment Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  deploy    - Deploy or update the Google Cloud Function"
	@echo "  logs      - View recent logs from the deployed function"
	@echo "  test      - Send test lookups to the deployed function"
	@echo "  unittest  - Run the unit tests"
	@echo "  help      - Show this help message (default)"

# The function only downloads public files, so it needs no permissions of its own. It runs as the BigQuery proxy's
# least-privilege service account rather than the project's default compute service account, which has broad access.
# 512MB gives an instance a third of a CPU, which reads a lookup file in well under a second.
deploy:
	gcloud functions deploy $(FUNCTION_NAME) \
		--project=cmg-analysis \
		--runtime python313 \
		--trigger-http \
		--allow-unauthenticated \
		--service-account=trexplorer-proxy-sa@cmg-analysis.iam.gserviceaccount.com \
		--entry-point $(FUNCTION_NAME) \
		--region=us-central1 \
		--memory 512MB \
		--timeout=60s \
		--gen2

logs:
	gcloud functions logs read $(FUNCTION_NAME) --project=cmg-analysis --region=us-central1 --limit=100

test:
	curl --silent --show-error --fail -X POST "$(FUNCTION_URL)" \
		-H "Content-Type: application/json" \
		-d '{"names": ["ATXN1", "ENSG00000124788", "nm_000332", "5S_rRNA", "FOOBAR"], "genome_version": "hg38", "catalog_table_id": "$(CATALOG_TABLE_ID)"}' | jq '.'

unittest:
	python3 -m unittest main_tests
