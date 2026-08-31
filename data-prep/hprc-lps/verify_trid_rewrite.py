"""Verify that rewriting variation cluster TRIDs changed nothing except TRID and STRUC.

Walks the original and rewritten VCFs together over a set of chromosomes and asserts,
record by record, that:

  - the two files hold the same records in the same order
  - every column other than INFO is byte-identical
  - within INFO, every key other than TRID and STRUC is byte-identical
  - an isolated repeat record (STRUC=<TR:...>) is completely unchanged
  - a variation cluster record has exactly the intended swap: the new TRID is
    "VC:" plus the span the old STRUC held, and the new STRUC is "<VC:" plus the
    old TRID plus ">"
  - no TRID repeats within the chromosome

Reads through tabix, so it is cheap to run on a few chromosomes rather than the whole
genome.
"""

import argparse
import itertools
import subprocess
import sys

# Distinguishes "this stream ended" from any real VCF line when pairing the two streams.
MISSING = object()


def parse_info(info):
    out = {}
    for kv in info.split(";"):
        key, _, value = kv.partition("=")
        out[key] = value
    return out


def records(vcf_path, chrom):
    proc = subprocess.Popen(["tabix", vcf_path, chrom], stdout=subprocess.PIPE,
                            text=True, bufsize=1 << 20)
    for line in proc.stdout:
        yield line.rstrip("\n")
    proc.stdout.close()
    if proc.wait() != 0:
        raise RuntimeError(f"tabix {vcf_path} {chrom} exited non-zero")


def check_chrom(original_vcf, rewritten_vcf, chrom):
    n_vc = n_tr = 0
    n = 0
    seen_trids = set()
    # Pair the two streams with zip_longest rather than zip: zip would consume a record from the
    # longer stream before noticing the other had ended, so a rewrite that dropped exactly one
    # record would leave both streams exhausted and look clean. The MISSING sentinel makes a length
    # difference visible on the record where it happens. Counting from the loop also avoids
    # re-running tabix afterwards just to compare totals, which would decompress both
    # multi-gigabyte files a second time.
    for n, (before, after) in enumerate(
            itertools.zip_longest(records(original_vcf, chrom), records(rewritten_vcf, chrom),
                                  fillvalue=MISSING), start=1):
        if before is MISSING:
            raise AssertionError(f"{chrom}: the rewrite has more records than the original "
                                 f"(the extra one is record {n})")
        if after is MISSING:
            raise AssertionError(f"{chrom}: the original has more records than the rewrite "
                                 f"(the missing one is record {n})")
        before_fields = before.split("\t")
        after_fields = after.split("\t")
        if len(before_fields) != len(after_fields):
            raise AssertionError(f"{chrom} record {n}: column count changed")
        for i, (b, a) in enumerate(zip(before_fields, after_fields)):
            if i != 7 and b != a:
                raise AssertionError(f"{chrom} record {n}: column {i + 1} changed:\n  {b[:120]}\n  {a[:120]}")

        before_info = parse_info(before_fields[7])
        after_info = parse_info(after_fields[7])
        if list(before_info) != list(after_info):
            raise AssertionError(f"{chrom} record {n}: INFO keys or their order changed")
        for key in before_info:
            if key not in ("TRID", "STRUC") and before_info[key] != after_info[key]:
                raise AssertionError(f"{chrom} record {n}: INFO/{key} changed")

        if before_info["STRUC"].startswith("<TR:"):
            n_tr += 1
            if before_fields[7] != after_fields[7]:
                raise AssertionError(f"{chrom} record {n}: isolated repeat record was modified")
        else:
            n_vc += 1
            span = before_info["STRUC"][4:-1]
            if after_info["TRID"] != f"VC:{span}":
                raise AssertionError(
                    f"{chrom} record {n}: expected TRID=VC:{span}, got {after_info['TRID']}")
            if after_info["STRUC"] != f"<VC:{before_info['TRID']}>":
                raise AssertionError(
                    f"{chrom} record {n}: expected STRUC=<VC:{before_info['TRID']}>, "
                    f"got {after_info['STRUC']}")

        trid = after_info["TRID"]
        if trid in seen_trids:
            raise AssertionError(f"{chrom}: TRID {trid} appears more than once after the rewrite")
        seen_trids.add(trid)

    if n == 0:
        # tabix exits 0 and prints nothing for a sequence name that is not in the index, so without
        # this a misspelled chromosome (say "19" for a VCF whose names are chr-prefixed) would
        # verify nothing and still report success.
        raise AssertionError(f"{chrom}: no records found; is the chromosome name spelled the way the VCF spells it?")
    return n, n_vc, n_tr


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--original-vcf", required=True)
    parser.add_argument("--rewritten-vcf", required=True)
    parser.add_argument("chroms", nargs="+", help="Chromosomes to check, e.g. chr19 chr21 chrY")
    args = parser.parse_args()

    for chrom in args.chroms:
        total, n_vc, n_tr = check_chrom(args.original_vcf, args.rewritten_vcf, chrom)
        print(f"{chrom}: {total:,d} records verified "
              f"({n_vc:,d} variation clusters rewritten, {n_tr:,d} isolated repeats unchanged)")
    print("All checks passed.")


if __name__ == "__main__":
    main()
