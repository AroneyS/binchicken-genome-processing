#!/usr/bin/env python3

"""
Author: Samuel Aroney
Process mmseqs2 protein clustering output
"""

import os
import sys
import argparse
import logging
import polars as pl
import numpy as np

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument("--mmseqs-output", help="mmseqs cluster file", required=True)
    parser.add_argument("--genome-magset", help="genome magset tsv", required=True)
    parser.add_argument("--output", help="Output prefix", required=True)

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    genome_magset = (
        pl.scan_csv(args.genome_magset, separator="\t")
        .with_columns(pl.col("genome").str.replace(r"_genomic$", ""))
    )

    logging.info("Loading clusters")
    clusters = (
        pl.scan_csv(args.mmseqs_output, separator="\t", has_header=False, new_columns=["rep", "protein"])
        .with_columns(
            pl.col("rep").hash(),
            pl.col("protein").hash(),
            genome = pl.col("protein").str.extract(r"^([^~]+)~").str.replace(r"^RS_|^GB_", "").str.replace(r"_genomic$", ""),
            )
        .join(genome_magset, on="genome")
        .collect()
    )

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    #####################
    ### General stats ###
    #####################
    # Total proteins: 609,622,874
    # Total protein clusters: 196,305,965

    logging.info("Writing general stats")
    (
        clusters
        .group_by(1)
        .agg(
            total = pl.count(),
            total_clusters = pl.col("rep").unique().count(),
            total_genomes = pl.col("genome").unique().count(),
            )
        .select("total", "total_clusters")
        .write_csv(args.output + "_stats.tsv", separator="\t")
    )

    ###################
    ### Proportions ###
    ###################
    # proportion of global gene clusters that are present in GTDB/GTDB+SPIRE/binchicken

    logging.info("Writing proportions")
    (
        clusters
        .with_columns(
            supermagset =
                pl.when(pl.col("magset") == "binchicken")
                .then(pl.lit("binchicken"))
                .when(pl.col("magset") == "GTDB")
                .then(pl.lit("GTDB"))
                .otherwise(pl.lit("SPIRE+")),
            protein = pl.lit(True),
            )
        .pivot(index="rep", columns="supermagset", values="protein", aggregate_function="first")
        .fill_null(False)
        .group_by("GTDB", "SPIRE+", "binchicken")
        .count()
        .sort("count", descending=True)
        .write_csv(args.output + "_prop.tsv", separator="\t")
    )

    # proportion of global gene clusters (from >=3 species) that are present in GTDB/GTDB+SPIRE/binchicken
    (
        clusters
        .filter(pl.count().over("rep") >= 3)
        .with_columns(
            supermagset =
                pl.when(pl.col("magset") == "binchicken")
                .then(pl.lit("binchicken"))
                .when(pl.col("magset") == "GTDB")
                .then(pl.lit("GTDB"))
                .otherwise(pl.lit("SPIRE+")),
            protein = pl.lit(True),
            )
        .pivot(index="rep", columns="supermagset", values="protein", aggregate_function="first")
        .fill_null(False)
        .group_by("GTDB", "SPIRE+", "binchicken")
        .count()
        .sort("count", descending=True)
        .write_csv(args.output + "_prop_gt3.tsv", separator="\t")
    )

    ###################
    ### Rarefaction ###
    ###################
    # Number of new protein clusters per new species (accumulation/rarefaction curves)
    def rarefy(proteins, genomes, column="rep", steps=100, iterations=5):
        max_size = genomes.height
        step_size = max_size // steps
        results = []

        for i in range(1, steps + 1):
            sample_size = step_size * i
            unique_counts = []
            for _ in range(iterations):
                sample_genomes = genomes.sample(n=sample_size)
                unique_counts.append(
                    proteins
                    .join(sample_genomes, on="genome")
                    .select(pl.col(column).unique().count())
                    .to_numpy()[0,0]
                )
            avg_unique = np.mean(unique_counts)
            results.append((sample_size, avg_unique))

        return pl.DataFrame(results, schema=["step", "unique"])

    logging.info("Writing novel binchicken rarefaction (>=3 species)")
    (
        rarefy(
            clusters.filter(pl.count().over("rep") >= 3).filter((pl.col("magset") == "binchicken").all().over("rep")),
            clusters.filter(pl.col("magset") == "binchicken").select("genome").unique(),
            )
        .write_csv(args.output + "_rarefy_bc_gt3.tsv", separator="\t")
    )

    logging.info("Writing novel binchicken rarefaction")
    (
        rarefy(
            clusters.filter((pl.col("magset") == "binchicken").all().over("rep")),
            clusters.filter(pl.col("magset") == "binchicken").select("genome").unique(),
            )
        .write_csv(args.output + "_rarefy_bc.tsv", separator="\t")
    )

    logging.info("Writing non-binchicken rarefaction (>=3 species)")
    (
        rarefy(
            clusters.filter(pl.count().over("rep") >= 3).filter(pl.col("magset") != "binchicken"),
            clusters.filter(pl.col("magset") != "binchicken").select("genome").unique(),
            steps=10,
            )
        .write_csv(args.output + "_rarefy_nonbc_gt3.tsv", separator="\t")
    )

    logging.info("Writing all rarefaction (>=3 species)")
    (
        rarefy(
            clusters.filter(pl.count().over("rep") >= 3),
            clusters.select("genome").unique(),
            steps=10,
            )
        .write_csv(args.output + "_rarefy_all_gt3.tsv", separator="\t")
    )

    logging.info("Writing non-binchicken rarefaction")
    (
        rarefy(
            clusters.filter(pl.col("magset") != "binchicken"),
            clusters.filter(pl.col("magset") != "binchicken").select("genome").unique(),
            steps=10,
            )
        .write_csv(args.output + "_rarefy_nonbc.tsv", separator="\t")
    )

    logging.info("Writing all rarefaction")
    (
        rarefy(
            clusters,
            clusters.select("genome").unique(),
            steps=10,
            )
        .write_csv(args.output + "_rarefy_all.tsv", separator="\t")
    )

    logging.info("Done")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
