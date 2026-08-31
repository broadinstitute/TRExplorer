"""Give variation cluster records their own TRIDs in a TRGT multisample VCF.

TRGT assumes TRIDs are unique within a catalog, but the v2 catalog gave a variation
cluster the TRID of the repeats it contains. A cluster holding a single repeat therefore
ends up with the same TRID as that repeat, and the two become indistinguishable in any
output that carries the TRID but not STRUC -- notably the trgt-lps table. See
https://github.com/PacificBiosciences/trgt-lps/issues/5.

This script swaps the two fields on variation cluster records only:

    before   TRID=19-103817-103838-AC   STRUC=<VC:19:103809-103838>
    after    TRID=VC:19:103809-103838   STRUC=<VC:19-103817-103838-AC>

which is the convention the v2.1 catalog already writes. The swap is lossless: the span
that STRUC held is still readable off the record's own POS and END, and the constituent
repeat IDs that TRID held move into STRUC. Isolated repeat records (STRUC=<TR:...>) are
passed through untouched, as is every other field of every record.

Reads a VCF on stdin and writes one on stdout, so it is meant to be run as

    gzcat in.vcf.gz | python3 rewrite_vc_trids_in_vcf.py | bgzip -@ 12 > out.vcf.gz

Fails rather than guessing if a record's INFO does not look the way it expects, and
verifies as it goes that the rewrite leaves every TRID unique.
"""

import argparse
import sys


def strip_chr(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def rewrite_info(info, chrom, pos):
    """Returns (new_info, kind) for one record's INFO field.

    kind is "VC" if the record was rewritten, "TR" if it was left alone.
    Raises ValueError if INFO is missing a field the rewrite depends on.
    """
    trid = end = struc = None
    for kv in info.split(";"):
        if kv.startswith("TRID="):
            trid = kv[5:]
        elif kv.startswith("END="):
            end = kv[4:]
        elif kv.startswith("STRUC="):
            struc = kv[6:]

    if trid is None or end is None or struc is None:
        raise ValueError(f"INFO is missing TRID, END or STRUC: {info}")

    if struc.startswith("<TR:"):
        return info, "TR"
    if not struc.startswith("<VC:"):
        raise ValueError(f"STRUC is neither <TR:...> nor <VC:...>: {struc}")
    if not struc.endswith(">"):
        raise ValueError(f"STRUC is not closed with '>': {struc}")

    span = struc[4:-1]
    expected_span = f"{strip_chr(chrom)}:{pos}-{end}"
    if span != expected_span:
        raise ValueError(
            f"STRUC span {span!r} does not match the record's own interval "
            f"{expected_span!r}; the rewrite would lose the span")

    new_trid = f"VC:{span}"
    new_struc = f"<VC:{trid}>"
    parts = []
    for kv in info.split(";"):
        if kv.startswith("TRID="):
            parts.append(f"TRID={new_trid}")
        elif kv.startswith("STRUC="):
            parts.append(f"STRUC={new_struc}")
        else:
            parts.append(kv)
    return ";".join(parts), "VC"


EXTRA_HEADER_LINES = [
    '##trexplorerTridRewrite="Variation cluster records were given their own TRID of the '
    'form VC:{chrom}:{start}-{end}, and the repeat IDs that TRID previously held were '
    'moved into STRUC as <VC:id1,id2,...>. Isolated repeat records are unchanged. '
    'See https://github.com/PacificBiosciences/trgt-lps/issues/5"',
]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--no-check-unique", dest="check_unique", action="store_false",
                        help="Skip the duplicate-TRID check, which holds every TRID in memory "
                             "(~1.5 GB for the HPRC256 VCF). Unique TRIDs are the whole point of "
                             "this script, so only skip the check if memory is genuinely short.")
    parser.add_argument("--progress-every", type=int, default=1000000,
                        help="Print a progress line to stderr every N records.")
    args = parser.parse_args()

    seen = set() if args.check_unique else None
    counts = {"VC": 0, "TR": 0}
    n_records = 0
    header_written = False

    out = sys.stdout
    for line in sys.stdin:
        if line.startswith("#"):
            # Emit the provenance lines just before the #CHROM line so they land after
            # every ##-prefixed header the source VCF carries.
            if line.startswith("#CHROM"):
                for extra in EXTRA_HEADER_LINES:
                    out.write(extra + "\n")
                header_written = True
            out.write(line)
            continue

        if not header_written:
            raise ValueError("Reached a record before the #CHROM header line")

        fields = line.split("\t", 8)
        if len(fields) < 9:
            raise ValueError(f"Record has fewer than 9 columns: {line[:200]}")
        new_info, kind = rewrite_info(fields[7], fields[0], fields[1])
        counts[kind] += 1
        n_records += 1
        if seen is not None:
            trid = new_info[new_info.index("TRID=") + 5:].split(";", 1)[0]
            if trid in seen:
                raise ValueError(f"Duplicate TRID after rewrite: {trid}")
            seen.add(trid)
        fields[7] = new_info
        out.write("\t".join(fields))

        if args.progress_every and n_records % args.progress_every == 0:
            print(f"  {n_records:,d} records ({counts['VC']:,d} VC, {counts['TR']:,d} TR)",
                  file=sys.stderr, flush=True)

    print(f"Rewrote {counts['VC']:,d} variation cluster records; "
          f"passed through {counts['TR']:,d} isolated repeat records "
          f"({n_records:,d} total)", file=sys.stderr)


if __name__ == "__main__":
    main()
