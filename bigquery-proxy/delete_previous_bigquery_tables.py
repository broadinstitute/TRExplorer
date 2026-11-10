"""Delete the old versions of the TRExplorer BigQuery tables, keeping the ones the website uses.

Each loader writes a new timestamped version of its table (eg. catalog_20261030_170523) and points the website at it
by rewriting a TABLE_ID constant in ../website/header_template.html. For each family of versions, this script keeps:

  - every version that header_template.html names, in the working copy, at HEAD, or at origin/main, so that the table
    the published site queries is kept even when the working copy has already moved to a newer one (or the reverse)
  - for catalog tables, the newest older version that the website has used before (named by a commit of
    header_template.html or index.html), as a fallback the site can be switched back to without a 60-90 minute reload.
    The loaders point the website only at full loads, so partial test loads (load_bigquery_main_table.py -n) never
    become the fallback.
  - non-empty versions newer than those: finished loads the website hasn't switched to yet
  - empty versions created less than a day ago, which may belong to a load that is still running (the main table
    loader creates its table about 20 minutes in and fills it at the end)

It deletes the rest, including empty tables left by loads that died partway, after asking for confirmation. Families
that the website names no version of are left alone.
"""

import datetime
import re
import subprocess
from pathlib import Path

from google.cloud import bigquery

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"

# Each entry is a prefix shared by a family of timestamped table versions (<prefix>_YYYYMMDD_HHMMSS)
TABLE_ID_PREFIXES = [
    "catalog",
    "hprc256_stratified_lps",
    "hprc256_stratified_total_allele_length",
    "hprc256_allele_purity",
    "hprc256_methylation",
    "hprc256_decomposed_alleles",
]

# The families that keep a fallback version older than the one the website uses
TABLE_ID_PREFIXES_THAT_KEEP_A_FALLBACK = {"catalog"}

# An empty table younger than this may belong to a load that is still running
MIN_AGE_OF_EMPTY_TABLES_TO_DELETE = datetime.timedelta(days=1)

REPO_DIR = Path(__file__).resolve().parent.parent
HEADER_TEMPLATE_PATH = "website/header_template.html"

# The git versions of header_template.html, besides the working copy, whose table IDs are kept: the last commit, and
# the last commit pushed to the branch that the website is published from
GIT_REFS_TO_READ_THE_WEBSITE_TABLE_IDS_FROM = ["HEAD", "origin/main"]

# The files whose TABLE_ID constants the loaders rewrite after a full load. Their git history tells which tables the
# website has used.
FILES_THAT_NAME_THE_WEBSITE_TABLES = [HEADER_TEMPLATE_PATH, "index.html"]


def get_table_id_prefix(table_id):
    """Returns the TABLE_ID_PREFIXES entry of a timestamped table version (<prefix>_YYYYMMDD_HHMMSS), or None for any
    other table."""
    for prefix in TABLE_ID_PREFIXES:
        if re.fullmatch(rf"{prefix}_\d{{8}}_\d{{6}}", table_id):
            return prefix
    return None


def find_table_ids_in_text(text):
    """Returns the set of timestamped table versions (<prefix>_YYYYMMDD_HHMMSS) that a text names."""
    prefixes = "|".join(TABLE_ID_PREFIXES)
    return set(re.findall(rf"\b(?:{prefixes})_\d{{8}}_\d{{6}}\b", text))


def read_website_table_ids():
    """Returns the set of table versions that ../website/header_template.html names in the working copy or at the git
    refs in GIT_REFS_TO_READ_THE_WEBSITE_TABLE_IDS_FROM (a ref that doesn't exist is skipped with a warning)."""
    website_table_ids = find_table_ids_in_text((REPO_DIR / HEADER_TEMPLATE_PATH).read_text())
    for git_ref in GIT_REFS_TO_READ_THE_WEBSITE_TABLE_IDS_FROM:
        result = subprocess.run(["git", "-C", str(REPO_DIR), "show", f"{git_ref}:{HEADER_TEMPLATE_PATH}"],
                                capture_output=True, text=True)
        if result.returncode != 0:
            print(f"WARNING: could not read {HEADER_TEMPLATE_PATH} at {git_ref}: {result.stderr.strip()}")
            continue
        website_table_ids |= find_table_ids_in_text(result.stdout)
    return website_table_ids


def read_table_ids_the_website_has_used():
    """Returns the set of table versions named by any commit in the git history of FILES_THAT_NAME_THE_WEBSITE_TABLES.

    Raises:
        SystemExit: if the git history can't be read
    """
    result = subprocess.run(["git", "-C", str(REPO_DIR), "log", "--all", "-p", "--format=", "-G", "TABLE_ID", "--",
                             *FILES_THAT_NAME_THE_WEBSITE_TABLES], capture_output=True, text=True)
    if result.returncode != 0:
        raise SystemExit(f"Could not read the git history of {', '.join(FILES_THAT_NAME_THE_WEBSITE_TABLES)}, which "
                         f"tells which table to keep as the fallback: {result.stderr.strip()}")
    return find_table_ids_in_text(result.stdout)


def decide_which_table_versions_to_keep(table_ids, website_table_ids, num_rows_by_table_id, created_by_table_id, now,
                                        table_ids_the_website_has_used=None):
    """Decides which versions of one family of tables to keep and which to delete.

    Args:
        table_ids (list): the family's table versions (<prefix>_YYYYMMDD_HHMMSS)
        website_table_ids (set): the table versions the website names (read_website_table_ids)
        num_rows_by_table_id (dict): table ID -> number of rows
        created_by_table_id (dict): table ID -> when the table was created (timezone-aware datetime)
        now (datetime.datetime): the current time (timezone-aware)
        table_ids_the_website_has_used (set): for a family that keeps a fallback, the table versions the website has
            used before (read_table_ids_the_website_has_used); the newest of them that is older than the ones the
            website names is kept. None to keep no fallback.

    Returns:
        list: (table ID, True to keep it or False to delete it, the reason) for each version, oldest first. If the
            website names none of the versions, all of them are kept.
    """
    table_ids = sorted(table_ids)  # timestamps are zero-padded, so this is oldest first
    used_by_website = [table_id for table_id in table_ids if table_id in website_table_ids]
    if not used_by_website:
        return [(table_id, True, "the website names no version of this table, so none are deleted")
                for table_id in table_ids]

    fallback_table_id = None
    if table_ids_the_website_has_used is not None:
        older_table_ids_the_website_has_used = [table_id for table_id in table_ids
                                                if table_id < used_by_website[0]
                                                and table_id in table_ids_the_website_has_used]
        if older_table_ids_the_website_has_used:
            fallback_table_id = older_table_ids_the_website_has_used[-1]

    decisions = []
    for table_id in table_ids:
        if table_id in website_table_ids:
            decisions.append((table_id, True, "the website uses it"))
        elif table_id == fallback_table_id:
            decisions.append((table_id, True, "fallback: the newest older version that the website has used"))
        elif num_rows_by_table_id[table_id] == 0:
            if now - created_by_table_id[table_id] < MIN_AGE_OF_EMPTY_TABLES_TO_DELETE:
                decisions.append((table_id, True, "empty, but created less than a day ago, so it may belong to a "
                                                  "load that is still running"))
            else:
                decisions.append((table_id, False, "empty, eg. left by a load that died partway"))
        elif table_id > used_by_website[-1]:
            decisions.append((table_id, True, "newer than the version the website uses: a finished load that the "
                                              "website hasn't switched to yet"))
        else:
            decisions.append((table_id, False, "an older version"))
    return decisions


def main():
    website_table_ids = read_website_table_ids()
    print(f"The website names these tables: {', '.join(sorted(website_table_ids))}")
    table_ids_the_website_has_used = read_table_ids_the_website_has_used()

    client = bigquery.Client(project=PROJECT_ID)
    dataset_ref = client.dataset(DATASET_ID)
    print(f"Listing tables in dataset {DATASET_ID}")
    tables = list(client.list_tables(dataset_ref))
    print(f"Found {len(tables)} existing tables in dataset {DATASET_ID}")

    table_ids_by_prefix = {prefix: [] for prefix in TABLE_ID_PREFIXES}
    for table in tables:
        prefix = get_table_id_prefix(table.table_id)
        if prefix:
            table_ids_by_prefix[prefix].append(table.table_id)

    to_delete = []
    now = datetime.datetime.now(datetime.timezone.utc)
    for prefix in TABLE_ID_PREFIXES:
        table_ids = table_ids_by_prefix[prefix]
        print(f"\nFound {len(table_ids)} table(s) with prefix '{prefix}':")
        details_by_table_id = {table_id: client.get_table(dataset_ref.table(table_id)) for table_id in table_ids}
        for table_id, keep, reason in decide_which_table_versions_to_keep(
                table_ids, website_table_ids,
                {table_id: details.num_rows for table_id, details in details_by_table_id.items()},
                {table_id: details.created for table_id, details in details_by_table_id.items()},
                now, table_ids_the_website_has_used if prefix in TABLE_ID_PREFIXES_THAT_KEEP_A_FALLBACK else None):
            print(f"  {'keep  ' if keep else 'DELETE'} {table_id} ({details_by_table_id[table_id].num_rows:,} rows): "
                  f"{reason}")
            if not keep:
                to_delete.append(table_id)

    if not to_delete:
        print("\nNo tables to delete.")
        return

    print(f"\nWill delete {len(to_delete)} table(s):")
    for table_id in to_delete:
        print(f"  {table_id}")
    print()
    user_input = input("Continue? [Y/n] ")
    if user_input.lower() in ("y", ""):
        for table_id in to_delete:
            print(f"Deleting table {table_id}...")
            client.delete_table(dataset_ref.table(table_id))
    print("Exiting...")


if __name__ == "__main__":
    main()
