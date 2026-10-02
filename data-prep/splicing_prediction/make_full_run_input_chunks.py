"""Write the input chunks for the full SpliceAI + SAI-10k run (spliceai_full_run_pipeline.py).

Takes every TRExplorer v2.1 locus that lies inside a GENCODE v50 transcript of the chosen gene set
(the loci SpliceAI would score; optionally limited to --max-distance-bp from the nearest GENCODE v50
splice site) and is polymorphic in HPRC256, builds each locus's population-based contraction and
expansion alleles (see make_tr_expansion_benchmark_alleles.py), and writes them in chunks (loci with no
alleles, because of an N in the anchor base or the tract, are left out and counted in the manifest):

    full_run_inputs/<gene_set>/chunk_<chrom>_<start>_<hash>.json.gz   {"gene_set", "chunk", "loci": [{"locus", "alleles", "target_labels"}, ...]}
    full_run_inputs/<gene_set>/manifest.json         counts, each chunk's name and SHA-256, and the
                                                     fingerprint of the allele design (see
                                                     compute_allele_design_fingerprint)

Loci are in sorted (chrom, start, end, motif) order, and the chunk boundaries are set by the loci
themselves: a new chunk starts at each locus whose ID hashes to a multiple of --mean-loci-per-chunk
(and a chunk past 3 times that many loci ends at the next locus whose ID hashes to a multiple of a
smaller number; see starts_a_new_chunk), and a chunk is named after its first locus. So rerunning
writes identical chunks, and adding or removing a locus (e.g. after an allele-design change) changes
only the chunk that holds it (almost always; see starts_a_new_chunk), not every later chunk as
fixed-size chunks would; the pipeline,
which records each chunk's SHA-256 with its outputs, then scores only the chunks that changed. Chunk
files from an earlier run of this script that are not in the new manifest are removed.

Run locally with a Python that has pyfaidx and numpy:
    python3 make_full_run_input_chunks.py --gene-set basic
"""
import argparse
import gzip
import hashlib
import json
import os

from make_tr_expansion_benchmark_alleles import (
    DISTANCES_TSV, HPRC256_ALLELE_LENGTH_STATS_TSV, MOTIF_RANGE_PERCENTILES_TSV, NO_DISTANCE_LIMIT_BP,
    SPLICEAI_ANNOTATION, load_polymorphic_scored_loci, simulate_alleles_for_locus)
from spliceai_l4_cost_and_disk_benchmark import FULL_RUN_LOCI_PER_CHUNK

HERE = os.path.dirname(os.path.abspath(__file__))


def allele_design_files(gene_set):
    """The files that decide which loci and alleles the inputs contain: the design code and its data."""
    return [os.path.join(HERE, "make_tr_expansion_benchmark_alleles.py"), os.path.join(HERE, "make_full_run_input_chunks.py"),
            HPRC256_ALLELE_LENGTH_STATS_TSV,
            MOTIF_RANGE_PERCENTILES_TSV, SPLICEAI_ANNOTATION[gene_set], DISTANCES_TSV]


def compute_allele_design_fingerprint(gene_set):
    """Returns the SHA-256 of the allele design files' contents, which the manifest records.

    spliceai_full_run_pipeline.py refuses to launch inputs whose recorded fingerprint differs from the
    current one, so inputs built under older rules are never scored. Regenerating the inputs after an
    edit that does not change the alleles (e.g. a comment) writes identical chunks with the same
    SHA-256s, so nothing already scored is redone.
    """
    sha256 = hashlib.sha256()
    for path in allele_design_files(gene_set):
        with open(path, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                sha256.update(block)
    return sha256.hexdigest()


MAX_LOCI_PER_CHUNK_IN_MEAN_CHUNKS = 3


def chunk_name(first_locus_id):
    """Names a chunk after its first locus, e.g. chr1-12345-12400-CAG -> chunk_chr1_000012345_<8 hex digits>."""
    chrom, start, _, _ = first_locus_id.rsplit("-", 3)
    return f"chunk_{chrom}_{int(start):09d}_{hashlib.sha256(first_locus_id.encode()).hexdigest()[:8]}"


def starts_a_new_chunk(locus_id, n_loci_in_current_chunk, mean_loci_per_chunk):
    """True if this locus (next in sorted order) should start a new chunk rather than join the current one.

    A locus starts a chunk when its ID hashes to a multiple of mean_loci_per_chunk (about 1 locus in
    mean_loci_per_chunk, decided by the ID alone). A chunk that already holds
    MAX_LOCI_PER_CHUNK_IN_MEAN_CHUNKS times that many loci also ends, but at the next locus whose hash is
    a multiple of mean_loci_per_chunk // 20 (about 50 loci later), not at a fixed count: a fixed count
    would let one added or removed locus shift every later capped boundary, while a cut point taken
    from the IDs almost always stays put (unless the cap is crossed right at such a locus).
    """
    locus_hash = int(hashlib.sha256(locus_id.encode()).hexdigest()[:15], 16)
    if locus_hash % mean_loci_per_chunk == 0:
        return True
    return (n_loci_in_current_chunk >= MAX_LOCI_PER_CHUNK_IN_MEAN_CHUNKS * mean_loci_per_chunk
            and locus_hash % max(2, mean_loci_per_chunk // 20) == 0)


def write_chunk(output_dir, gene_set, loci_with_alleles):
    """Writes one chunk file and returns its manifest entry: name, SHA-256 of the file, and counts.

    The gzip header's timestamp is fixed so identical contents give an identical file and checksum.
    """
    name = chunk_name(loci_with_alleles[0]["locus"])
    content = json.dumps({"gene_set": gene_set, "chunk": name, "loci": loci_with_alleles}, separators=(",", ":")).encode()
    data = gzip.compress(content, mtime=0)
    with open(os.path.join(output_dir, f"{name}.json.gz"), "wb") as f:
        f.write(data)
    return {
        "chunk": name,
        "sha256": hashlib.sha256(data).hexdigest(),
        "n_loci": len(loci_with_alleles),
        "n_alleles": sum(len(locus["alleles"]) for locus in loci_with_alleles),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gene-set", choices=["basic", "comprehensive"], required=True)
    parser.add_argument("--mean-loci-per-chunk", type=int, default=FULL_RUN_LOCI_PER_CHUNK)
    parser.add_argument("--fasta", default=os.path.expanduser("~/hg38.fa"))
    parser.add_argument("--max-distance-bp", type=int, default=NO_DISTANCE_LIMIT_BP,
                        help="Only loci at most this far from the nearest splice site (default: no limit)")
    args = parser.parse_args()

    import pyfaidx
    fasta = pyfaidx.Fasta(args.fasta)
    counts, loci = load_polymorphic_scored_loci(args.gene_set, max_distance_bp=args.max_distance_bp)

    output_dir = os.path.join(HERE, "full_run_inputs", args.gene_set)
    os.makedirs(output_dir, exist_ok=True)
    chunks, pending, n_loci_without_alleles = [], [], 0
    for i, locus in enumerate(loci):
        labeled_alleles = simulate_alleles_for_locus(fasta, *locus)
        if not labeled_alleles:
            # An N in the anchor base or the tract: the production scorer rejects such alleles
            n_loci_without_alleles += 1
            continue
        locus_id = f"{locus[0]}-{locus[1]}-{locus[2]}-{locus[3]}"
        if pending and starts_a_new_chunk(locus_id, len(pending), args.mean_loci_per_chunk):
            chunks.append(write_chunk(output_dir, args.gene_set, pending))
            pending = []
            if len(chunks) % 50 == 1:
                print(f"wrote {chunks[-1]['chunk']} ({i:,} of {len(loci):,} loci)", flush=True)
        # target_labels[i]: the target lengths (make_tr_expansion_benchmark_alleles.TARGET_LABELS) alleles[i]
        # stands for; the pipeline copies them into each stored record
        pending.append({"locus": locus_id, "alleles": [variant for variant, _ in labeled_alleles],
                        "target_labels": [labels for _, labels in labeled_alleles]})
    if pending:
        chunks.append(write_chunk(output_dir, args.gene_set, pending))
    # Chunk files from an earlier run that are not in this manifest (their first locus changed)
    current_files = {f"{chunk['chunk']}.json.gz" for chunk in chunks}
    for name in os.listdir(output_dir):
        if name.startswith("chunk_") and name.endswith(".json.gz") and name not in current_files:
            os.remove(os.path.join(output_dir, name))

    manifest = {
        "gene_set": args.gene_set,
        "allele_design_fingerprint": compute_allele_design_fingerprint(args.gene_set),
        "mean_loci_per_chunk": args.mean_loci_per_chunk,
        "max_distance_bp": args.max_distance_bp,
        **counts,
        "n_loci_skipped_for_non_acgt_bases": n_loci_without_alleles,
        "n_loci": sum(chunk["n_loci"] for chunk in chunks),
        "n_alleles": sum(chunk["n_alleles"] for chunk in chunks),
        "chunks": chunks,
    }
    with open(os.path.join(output_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"Wrote {len(chunks)} chunks, {manifest['n_loci']:,} loci and {manifest['n_alleles']:,} alleles to {output_dir}")


if __name__ == "__main__":
    main()
