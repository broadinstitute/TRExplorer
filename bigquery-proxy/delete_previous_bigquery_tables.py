from google.cloud import bigquery

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"
TABLE_ID = "catalog"

client = bigquery.Client(project=PROJECT_ID)
dataset_ref = client.dataset(DATASET_ID)

# List existing tables in the dataset
print(f"Listing tables in dataset {DATASET_ID}")
tables = client.list_tables(dataset_ref)
table_ids = []
if tables:
    tables = list(tables)
    print(f"Found {len(tables)} existing tables:")
    for t in tables:
        if t.table_id.startswith(TABLE_ID):
            print(f"  {t.table_id}")
            table_ids.append(t.table_id)
else:
    print(f"No existing tables found in dataset {DATASET_ID}")

table_ids.sort()

if len(table_ids) > 1:
    print(f"Will delete the first {len(table_ids) - 1} table(s):")
    for table_id in table_ids[:-1]:
        print(f"  {table_id}")
    print()
    user_input = input("Continue? [Y/n] ")
    if user_input.lower() == "y":
        for table_id in table_ids[:-1]:
            print(f"Deleting table {table_id}...")
            table_ref = dataset_ref.table(table_id)
            client.delete_table(table_ref)
elif len(table_ids) == 1:
    print(f"Only one table found in dataset {DATASET_ID}: {table_ids[0]}. Keeping it.")
else:
    print(f"No tables found in dataset {DATASET_ID}")

print("Exiting...")