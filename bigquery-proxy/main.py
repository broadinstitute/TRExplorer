import datetime
from flask import jsonify
import gzip
import json
import os
import re
import uuid

from functions_framework import create_app
from google.cloud import bigquery, storage

BIGQUERY_PROJECT = "cmg-analysis"
BIGQUERY_DATASET = "tandem_repeat_explorer"


def return_query_results(client, result_table, start_index=0, page_size=100):
    """Return query results from BigQuery.
    
    Args:
        client: The BigQuery client
        result_table: The BigQuery table to return results from
        start_index: The index of the first row to return (default: 0) 
        page_size: The number of rows to return (default: 100)
        
    Returns:
        A tuple of (response_body, response_headers)
    """
    
    try:
        # see docs @ https://cloud.google.com/python/docs/reference/bigquery/latest/google.cloud.bigquery.client.Client#google_cloud_bigquery_client_Client_query_and_wait
        # query(..) returns agoogle.cloud.bigquery.job.query.QueryJob
        # https://cloud.google.com/bigquery/docs/paging-results#page_through_query_results

        row_iterator = client.list_rows(result_table, max_results=page_size, start_index=start_index)

        response_dict = {
            "rows": [dict(row) for row in row_iterator],
            "current_page": start_index // page_size + 1,
            "total_pages": result_table.num_rows // page_size + 1,
            "total_results": result_table.num_rows,
            "page_size": page_size,
            "current_page_start_index": start_index,
            "next_page_start_index": start_index + page_size if start_index + page_size < result_table.num_rows else None,
        }

        return jsonify(response_dict), 200

    except Exception as e:
        print(f"ERROR: Unexpected error: {str(e)}")
        response = jsonify({"error": f"Unexpected error: {e}"})
        return response, 500


def export_to_file(client, result_table, export_to_file_format, tool_name=None):
    """Export query results to GCS in the specified format.
    
    Args:
        client: The BigQuery client
        result_table: The BigQuery table to export results from
        export_to_file_format: The format to export to. Must be one of: 'TSV', 'BED', 'JSON'
        
    Returns:
        A tuple of (response_body, response_headers)
    """
    # Generate a unique export path
            
    if export_to_file_format not in ("TSV", "BED", "JSON"):
        response_dict = {"error": f"Invalid export_to_file_format param: {export_to_file_format}. It must be one of: TSV, BED, JSON"}
        print(f"ERROR: {response_dict['error']}")
        return jsonify(response_dict), 400
    
    filename_prefix = uuid.uuid4().hex.upper()
    bucket_name = "tandem-repeat-explorer-export"
    uri = f'gs://{bucket_name}/{filename_prefix}_*.{export_to_file_format.lower()}.gz'
    if export_to_file_format == "JSON":
        destination_format = bigquery.DestinationFormat.NEWLINE_DELIMITED_JSON
        print_header = None
        field_delimiter = None
    else:
        destination_format = bigquery.DestinationFormat.CSV
        print_header = export_to_file_format == "TSV"
        field_delimiter = '\t'

    try:
        # Run the export job
        job_config = bigquery.ExtractJobConfig(
            compression=bigquery.Compression.GZIP,
            print_header=print_header,
            field_delimiter=field_delimiter,
            destination_format=destination_format,
        )

        job = client.extract_table(result_table, destination_uris=uri, job_config=job_config)
        job.result()  # Wait for job to complete

        # Get the list of one or more result files
        storage_client = storage.Client()
        bucket = storage_client.bucket(bucket_name)
        blobs = list(bucket.list_blobs(prefix=filename_prefix))
        if len(blobs) == 0:
            # return an error
            response_dict = {"error": f"Export to file did not complete successfully."}
            print(f"ERROR: {response_dict['error']}")
            return jsonify(response_dict), 400

        public_urls = []
        for i, output_blob in enumerate(blobs):
            print(f"DEBUG: Processing output file {i+1} of {len(blobs)}: gs://{bucket_name}/{output_blob.name}")

            if export_to_file_format == "JSON" and len(blobs) == 1:
                # Download the original file
                temp_file = f"/tmp/{filename_prefix}_original.json.gz"
                output_blob.download_to_filename(temp_file)

                new_filename = f"{filename_prefix}_converted.json.gz"
                new_temp_file = f"/tmp/{new_filename}"

                # Read and convert the newline-delimited JSON to a proper JSON array
                with gzip.open(temp_file, 'rt') as f, gzip.open(new_temp_file, 'wt') as f_out:
                    f_out.write('[')
                    for i, line in enumerate(f):
                        row = json.loads(line)
                        if i > 0:
                            f_out.write(', ')
                        f_out.write(json.dumps(row, indent=4))
                    f_out.write(']\n')

                # Upload the converted file to GCS and clean up temporary files
                new_blob = bucket.blob(new_filename)
                new_blob.upload_from_filename(new_temp_file)

                os.remove(temp_file)
                os.remove(new_temp_file)
                output_blob.delete()

                # Update the blob reference to the new file
                output_blob = new_blob

            # Set the Content-Disposition metadata to specify the download filename
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            filename_label = f"_for_{tool_name}" if tool_name else ""
            if len(blobs) > 1:
                filename_label += f".shard_{i+1:02d}_of_{len(blobs):02d}"
            output_filename = f"TR_catalog{filename_label}.{result_table.num_rows}_loci.{timestamp}.{export_to_file_format.lower()}.gz"
            output_blob.content_type = 'application/gzip'
            output_blob.content_disposition = f'attachment; filename="{output_filename}"'
            output_blob.patch()

            public_urls.append(f"https://storage.cloud.google.com/{bucket_name}/{output_blob.name}")

        return jsonify({
            'status': 'success', 
            'total_results': result_table.num_rows,
            'public_urls': public_urls,
        }), 200
        
    except Exception as e:
        print(f"ERROR: Unexpected error: {str(e)}")
        return jsonify({"error": f"Unexpected error: {e}"}), 500
        


def query_db(request):
    """Cloud Function for forwarding SQL queries to BigQuery.
    
    Args:
        request (flask.Request): The request object.

        Parameters:
            sql (str): The SQL query to execute.
            start_index (int): (optional)The index of the first row to return (default: 0).
            page_size (int): (optional) The number of rows to return (default: 100).
        
    Returns:
        flask.Response: The response object.
    
    """

    response_headers = {
        'Access-Control-Allow-Origin': '*',
        'Access-Control-Allow-Methods': 'POST',
        'Access-Control-Allow-Headers': 'Content-Type',
    }
  
    # set CORS headers for the preflight request
    if request.method == 'OPTIONS':
        return '', 204, response_headers
    
    # validate request body
    try:
        data = request.get_json()
    except Exception as e:
        response_dict = {"error": f"Failed to parse JSON: {request.get_data(as_text=True)}"}
        print(f"ERROR: {response_dict['error']}")
        return jsonify(response_dict), 400, response_headers
    

    sql = data.get("sql")
    if not sql:
        response_dict = {"error": f"Missing SQL query in JSON: {request.get_data(as_text=True)}"}
        print(f"ERROR: {response_dict['error']}")
        return jsonify(response_dict), 400, response_headers
    
    sql = str(sql).strip()

    pattern = rf"^SELECT\s+([\w\*,\s()\"'=;\.\+\-\*\/]+)\s+FROM[\s`]+{BIGQUERY_PROJECT}\.{BIGQUERY_DATASET}"
    if not re.search(pattern, sql, re.IGNORECASE | re.VERBOSE | re.MULTILINE):
        response_dict = {"error": f"Invalid SQL query: {sql}. It must be a SELECT from the allowed project and dataset."}
        print(f"ERROR: {response_dict['error']}")
        return jsonify(response_dict), 400, response_headers

    client = bigquery.Client()
    job = client.query(sql, job_config=bigquery.QueryJobConfig(use_query_cache=True))
    job.result()  # wait for the query to complete
    result_table = client.get_table(job.destination)  # get the temporary destination table object

    export_to_file_format = data.get("export_to_file_format")
    if export_to_file_format:
        # Export to file
        tool_name = data.get("tool_name")  # optional tool name to add to the filename
        if not tool_name or tool_name.lower() not in ("eh", "expansionhunter", "expansion_hunter", "gangstr", "hipstr", "longtr", "trgt"):
            tool_name = None

        response_json, response_status_code = export_to_file(
            client, result_table, export_to_file_format, tool_name=tool_name)
        
        return response_json, response_status_code, response_headers


    else:
        # Return a page of results
        min_page_size = 1
        default_page_size = 100
        max_page_size = 10_000

        try:
            start_index = int(data.get("start_index", 0))
        except Exception as e:
            response_dict = {"error": f"Invalid start_index param: {e}"}
            print(f"ERROR: {response_dict['error']}")
            return jsonify(response_dict), 400, response_headers

        try:
            page_size = max(min_page_size, min(max_page_size, int(data.get("page_size", default_page_size))))
        except Exception as e:
            response_dict = {"error": f"Invalid page_size param: {e}"}
            print(f"ERROR: {response_dict['error']}")
            return jsonify(response_dict), 400, response_headers

        response_json, response_status_code = return_query_results(
            client, result_table, start_index, page_size)

        return response_json, response_status_code, response_headers

# This is required for Cloud Functions Gen2
if __name__ == "__main__":
    port = int(os.getenv("PORT", "8080"))
    print(f"Starting server on port {port}")
    app = create_app("query_db")
    app.run(host="0.0.0.0", port=port)
