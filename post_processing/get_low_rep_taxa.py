#!/usr/bin/env python3

"""
Author: Samuel Aroney
Return phyla with low representation in GTDB
"""

import sys
import argparse
import logging
import polars as pl

def pipeline(gtdb):
    low_rep = (
        gtdb
        .with_columns(
            pl.col("gtdb_taxonomy")
                .str.split_exact(";", 7)
                .struct.rename_fields(["domain", "phylum", "class", "order", "family", "genus", "species"])
                .alias("split_taxonomy")
            )
        .unnest("split_taxonomy")
        .groupby("domain", "phylum")
        .count()
        .filter(pl.col("count") <= 10)
        .sort(["domain", "phylum"])
    )

    return low_rep

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debug", help="output debug information", action="store_true")
    parser.add_argument("--quiet", help="only output errors", action="store_true")

    parser.add_argument("--archaea-metadata", help="Archaeal metadata file from GTDB", required=True)
    parser.add_argument("--bacteria-metadata", help="Bacterial metadata file from GTDB", required=True)
    parser.add_argument("--output", help="Output file", required=True)

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format="%(asctime)s %(levelname)s: %(message)s", datefmt="%Y/%m/%d %I:%M:%S %p")

    metadata_cols = ["accession", "gtdb_representative", "gtdb_taxonomy"]
    gtdb = pl.concat([
        pl.read_csv(args.archaea_metadata, separator="\t", ignore_errors=True).select(metadata_cols),
        pl.read_csv(args.bacteria_metadata, separator="\t", ignore_errors=True).select(metadata_cols)
    ])

    low_rep = pipeline(gtdb)
    low_rep.write_csv(args.output, separator="\t")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
