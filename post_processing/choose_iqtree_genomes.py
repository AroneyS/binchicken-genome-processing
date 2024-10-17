#!/usr/bin/env python3

"""
Author: Samuel Aroney
Choose genomes for IQtree (1 per class)
"""

import os
import sys
import argparse
import logging
import polars as pl
import gzip


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument("--node-names", help="Node names from name_clades.py", required=True)
    parser.add_argument("--genome-taxonomy", help="Genome taxonomy from name_clades.py", required=True)
    parser.add_argument("--gtdb-reps", help="GTDB representatives for known clades", required=True)
    parser.add_argument("--alignment", help="GTDBtk de novo alignment file (fasta.gz)", required=True)
    parser.add_argument("--clade-level", help="Clade level for choosing genomes", required=True)
    parser.add_argument("--outgroup", help="Outgroup taxonomy for tree generation")
    parser.add_argument("--output", help="Output folder", required=True)

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    gtdb_reps = pl.read_csv(args.gtdb_reps, separator="\t")
    clade_letter = args.clade_level[0].lower()

    rep_genomes = pl.concat([
        pl.read_csv(args.node_names, separator="\t")
        .filter(pl.col("clade").str.starts_with(clade_letter))
        .select(genome = "genome_rep"),
        gtdb_reps
        .select(genome = "#Genome ID")
    ])

    os.makedirs(args.output, exist_ok=True)
    logging.info(f"Saving chosen genomes (1 per {args.clade_level})")
    rep_genomes.write_csv(args.output + "/rep_genomes.tsv")

    rep_genomes_dict = {c: False for c in rep_genomes.get_column("genome").to_list()}
    with open(args.output + "/alignment.msa.fasta", 'wb') as of:
        with gzip.open(args.alignment, 'rb') as f:
            print_line = False
            for line in f:
                if line.startswith(b">"):
                    genome = line.decode("utf-8").strip().split(">")[1].split(" ")[0]
                    if genome in rep_genomes_dict:
                        print_line = True
                        rep_genomes_dict[genome] = True
                    else:
                        print_line = False
                if print_line:
                    of.write(line)

    for genome, aligned in rep_genomes_dict.items():
        if not aligned:
            logging.warning(f"Alignment not found for genome: {genome}")

    genome_taxonomy = pl.concat([
        gtdb_reps.select(genome = "#Genome ID", taxonomy = "GTDB taxonomy"),
        pl.read_csv(args.genome_taxonomy, separator="\t"),
    ])

    if args.outgroup:
        logging.info("Recording outgroup genomes")
        outgroup_taxa = (
            rep_genomes
            .join(genome_taxonomy, on="genome")
            .filter(pl.col("taxonomy").str.contains(args.outgroup))
            .get_column("genome")
            .to_list()
        )
        with open(args.output + "/outgroup_taxa.tsv", 'w') as f:
            f.write(",".join(outgroup_taxa))

    logging.info("Done")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
