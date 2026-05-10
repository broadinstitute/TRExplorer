This directory defines a Google Cloud Function that proxies requests between the website and the BigQuery API. Even though the underlying BigQuery table is public, the proxy is needed to make auth work. 

Files:

- **load_bigquery_main_table.py** - script that creates and populates the primary BigQuery table which underlies all searches in TR Explorer
- **load_bigquery_HPRC256_stratified_LPS_frequency_data.py** - populates the per-locus HPRC256 stratified LPS allele-frequency table (per-population, per-sex, and Pop_Sex strata)
- **load_bigquery_HPRC256_stratified_allele_purity_and_methylation_data.py** - populates two HPRC256 tables: stratified allele-size-and-purity distributions and stratified allele-size-and-methylation distributions
- **main.py** - the Flask server code that handles SQL query requests (see `query_db(..)` endpoint), and export to file requests (see `export_to_file(..)` endpoint).
- **requirements.txt** - the python dependencies of main.py
- **Makefile** - contains the following targets:
```
  deploy    - Deploy or update the Google Cloud Function
  logs      - View recent logs from the deployed function
  test      - Run test queries against the deployed function
  undeploy  - Remove the deployed function and its logs
  load      - Load data into BigQuery by running load_bigquery_main_table.py
```
