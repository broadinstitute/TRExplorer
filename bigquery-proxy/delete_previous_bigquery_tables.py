from google.cloud import bigquery

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"

# Each entry is a prefix shared by a family of timestamped table versions. For each
# family we keep the latest (lexicographically highest, since timestamps are zero-padded
# YYYYMMDD_HHMMSS) and delete the rest.
TABLE_ID_PREFIXES = [
    "catalog",
    "hprc256_stratified_lps",
    "hprc256_allele_purity",
    "hprc256_methylation",
]

client = bigquery.Client(project=PROJECT_ID)
dataset_ref = client.dataset(DATASET_ID)

print(f"Listing tables in dataset {DATASET_ID}")
tables = list(client.list_tables(dataset_ref))
print(f"Found {len(tables)} existing tables in dataset {DATASET_ID}")

prefix_to_table_ids = {prefix: [] for prefix in TABLE_ID_PREFIXES}
for t in tables:
    for prefix in TABLE_ID_PREFIXES:
        if t.table_id.startswith(prefix):
            prefix_to_table_ids[prefix].append(t.table_id)
            break

to_delete = []
for prefix in TABLE_ID_PREFIXES:
    table_ids = sorted(prefix_to_table_ids[prefix])
    print(f"\nFound {len(table_ids)} table(s) with prefix '{prefix}':")
    for table_id in table_ids:
        print(f"  {table_id}")
    if len(table_ids) > 1:
        to_delete.extend(table_ids[:-1])
    elif len(table_ids) == 1:
        print(f"  -> keeping the only version")
    else:
        print(f"  -> nothing to delete")

if to_delete:
    print(f"\nWill delete {len(to_delete)} table(s):")
    for table_id in to_delete:
        print(f"  {table_id}")
    print()
    user_input = input("Continue? [Y/n] ")
    if user_input.lower() in ("y", ""):
        for table_id in to_delete:
            print(f"Deleting table {table_id}...")
            client.delete_table(dataset_ref.table(table_id))
else:
    print("\nNo tables to delete.")

print("Exiting...")
