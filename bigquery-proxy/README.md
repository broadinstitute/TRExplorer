This directory defines a Google Cloud Function that proxies requests between the website and the BigQuery API. Even though the underlying BigQuery table is public, the proxy is needed to make auth work. 

Files:

- **load_data_into_bigquery.py** - script that creates and populates the primary BigQuery table which underlies all searches in TR Explorer
- **main.py** - the Flask server code that handles SQL query requests (see `query_db(..)` endpoint), and export to file requests (see `export_to_file(..)` endpoint).
- **requirements.txt** - the python dependencies of main.py
- **Makefile** - contains the following targets:
```
  deploy    - Deploy or update the Google Cloud Function
  logs      - View recent logs from the deployed function
  test      - Run test queries against the deployed function
  undeploy  - Remove the deployed function and its logs
  load      - Load data into BigQuery by running load_data_into_bigquery.py
```
