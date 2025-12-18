"""Output a TSV where, for each locus in the TRExplorer catalog, it contains values from overlapping RepeatMasker track entries.

The RepeatMasker track table includes the following columns:

$1        #bin : 2
$2     swScore : 536
$3    milliDiv : 349
$4    milliDel : 146
$5    milliIns : 56
$6    genoName : chr1
$7   genoStart : 92274205
$8     genoEnd : 92275925
$9    genoLeft : -156680497
$10     strand : +
$11    repName : L2
$12   repClass : LINE
$13  repFamily : L2
$14   repStart : 406
$15     repEnd : 2306
$16    repLeft : -1113
$17         id : 1
"""

import argparse 
import collections
import gzip
import intervaltree
import ijson
import json
import os
import pandas as pd
import pprint
import tqdm

from str_analysis.utils.misc_utils import parse_interval

parser = argparse.ArgumentParser()
parser.add_argument(
    "--trexplorer-catalog",
    default="~/code/tandem-repeat-catalogs/results__2025-12-07/release_draft_2025-12-07/repeat_catalog_v2.hg38.1_to_1000bp_motifs.EH.json.gz")
parser.add_argument("--repeat-masker-track", default="./hg38.RepeatMasker.tsv.gz")
parser.add_argument("--output-json", default="hg38.RepeatMasker.lookup.json.gz")
parser.add_argument("-n", type=int, default=None)
args = parser.parse_args()

# create interval trees for repeat masker track
print(f"Reading repeat masker table: {os.path.basename(args.repeat_masker_track)}")
repeat_masker_df = pd.read_table(args.repeat_masker_track)
print(f"Parsed {len(repeat_masker_df):,d} rows from repeat masker table: {os.path.basename(args.repeat_masker_track)}")

pprint.pprint(repeat_masker_df.repFamily.value_counts())
pprint.pprint(repeat_masker_df.repClass.value_counts())

excluded_classes = {  # excluded classes
    #"SINE",              #1910631
    #"LINE",              #1614481
    #"LTR",               #770551
    #"Simple_repeat",     #724562
    "DNA",               #512404
    #"Low_complexity",    #106053
    #"Satellite",         #9133
    #"Retroposon",        #5974
    #"LTR?",              #5777
    "Unknown",           #5732
    #"snRNA",             #4686
    "DNA?",              #3364
    #"tRNA",              #2164
    #"rRNA",              #1953
    "RC",                #1820
    #"srpRNA",            #1745
    #"scRNA",             #1484
    "RNA",               #721
    "RC?",               #417
    #"SINE?",             #38
}

before = len(repeat_masker_df)
repeat_masker_df = repeat_masker_df[~repeat_masker_df.repClass.isin(excluded_classes)]

print(f"Filtered out {before - len(repeat_masker_df):,d} rows ({100 * (before - len(repeat_masker_df)) / before:.2f}%) because they were in classes: {', '.join(sorted(excluded_classes))}")

print(f"Building interval trees for repeat masker table")
repeat_masker_interval_trees = collections.defaultdict(intervaltree.IntervalTree)
for row_tuple in tqdm.tqdm(repeat_masker_df.itertuples(), unit=" rows", unit_scale=True, total=len(repeat_masker_df)):
    interval = intervaltree.Interval(int(row_tuple.genoStart), int(row_tuple.genoEnd), data={  # genoStart is 0-based
        "repName": row_tuple.repName,
        "repClass": row_tuple.repClass,
        "repFamily": row_tuple.repFamily,
        "milliDiv": row_tuple.milliDiv,
    })
    repeat_masker_interval_trees[row_tuple.genoName.replace("chr", "")].add(interval)


# iterate through the trexplorer catalog using ijson
max_overlapping_entries = 0
print(f"Parsing TRExplorer catalog: {os.path.basename(args.trexplorer_catalog)}")
fopen = gzip.open if args.output_json.endswith("gz") else open
with gzip.open(os.path.expanduser(args.trexplorer_catalog), "rt") as f, fopen(args.output_json, "wt") as f_out:
    output_counter = 0
    f_out.write("{\n")
    for i, item in tqdm.tqdm(enumerate(ijson.items(f, "item", use_float=True)), unit=" loci", unit_scale=True, total=args.n if args.n is not None else 4_863_041):
        if args.n is not None and i >= args.n:
            break
        locus_id = item["LocusId"]
        reference_region = item["ReferenceRegion"]
        chrom, start, end = parse_interval(reference_region)
        chrom = chrom.replace("chr", "")
        overlapping_intervals = repeat_masker_interval_trees[chrom].overlap(start - 1, end + 1)  # do +/-1 to capture immediately-adjacent intervals
        max_overlapping_entries = max(max_overlapping_entries, len(overlapping_intervals))
        if len(overlapping_intervals) > 0:
            if output_counter > 0:
                f_out.write(",\n")
            output_counter += 1
            f_out.write(f"  \"{locus_id}\": ")
            f_out.write(json.dumps([interval.data for interval in overlapping_intervals]))
            
    f_out.write("}\n")

print(f"Max overlapping entries at a locus: {max_overlapping_entries:,d}")
print(f"Done generating lookup table for {output_counter:,d} loci from TRExplorer catalog")

# write BED file 
output_bed = args.output_json.replace(".json", ".bed").replace(".gz", "").replace(".lookup", "")
print(f"Writing BED file of repeat masker intervals to {output_bed}")
line_counter = 0
with open(output_bed, "wt") as f_out:
    for chrom, interval_tree in tqdm.tqdm(repeat_masker_interval_trees.items(), unit=" chromosomes", unit_scale=True):
        for interval in sorted(interval_tree):
            data = interval.data
            name_field = f"{data['repFamily']}:{data['repName']}"  # f"({data['repClass']}:{data['repFamily']}:{data['repName']})"
            f_out.write("\t".join(map(str, [chrom, interval.begin, interval.end, name_field, ".", "."])) + "\n")
            line_counter += 1

print(f"Wrote {line_counter:,d} lines to BED file: {output_bed}")

os.system(f"bgzip -f {output_bed}")
os.system(f"tabix -f {output_bed}.gz")
print(f"Generated indexed BED file: {output_bed}.gz")
