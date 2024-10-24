#!/usr/bin/env python3

"""
Author: Samuel Aroney
Calculate rarefaction curves for genomes
"""

import os
import sys
import argparse
import logging
import polars as pl
import numpy as np


def rarefy(taxonomy, genomes, column="clade", join="name", steps=100, iterations=10):
    max_size = genomes.height
    step_size = max_size // steps
    integer_steps = max_size // step_size
    results = []

    for i in range(1, integer_steps + 1):
        sample_size = step_size * i
        unique_counts = []
        for _ in range(iterations):
            sample_genomes = genomes.sample(n=sample_size, with_replacement=False)
            unique_counts.append(
                taxonomy
                .join(sample_genomes, on=join)
                .select(pl.col(column).unique().count())
                .to_numpy()[0,0]
            )
        avg_unique = np.mean(unique_counts)
        results.append((sample_size, avg_unique))

    return pl.DataFrame(results, schema=["step", "unique"])

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument("--genome-metadata", help="Genome metadata file", required=True)
    parser.add_argument("--coassembly-metadata", help="Coassembly metadata file", required=True)
    parser.add_argument("--archaeal-named", help="Archaeal named taxonomy", required=True)
    parser.add_argument("--bacterial-named", help="Bacterial named taxonomy", required=True)
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

    named_taxonomy = pl.concat([
        pl.read_csv(args.archaeal_named, separator="\t"),
        pl.read_csv(args.bacterial_named, separator="\t"),
    ])

    coassembly_metadata = (
        pl.read_csv(args.coassembly_metadata, separator="\t")
        .select("samples", coassembly = "name")
        # Remove single-sample genomes from rarefaction
        .filter(pl.col("samples").str.contains(","))
    )

    genome_taxonomy = (
        pl.read_csv(args.genome_metadata, separator="\t")
        .select("coassembly", "name", "taxonomy")
        .join(named_taxonomy, left_on="name", right_on="genome", how="left")
        .join(coassembly_metadata, on="coassembly")
        .with_columns(
            taxonomy = 
                pl.when(pl.col("taxonomy_right").is_null())
                .then(pl.col("taxonomy"))
                .otherwise(pl.col("taxonomy_right"))
            )
        .filter(pl.col("taxonomy").is_not_null())
        .with_columns(
            clade = pl.col("taxonomy").str.split(";"),
            domain = pl.col("taxonomy").str.split(";").list.first(),
            )
        .explode("clade")
        .with_columns(taxon = pl.col("clade").str.slice(0, 1))
        .select("domain", "samples", "name", "taxon", "clade")
    )

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    ####################
    ### Coassemblies ###
    ####################
    # Phyla
    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").filter(pl.col("taxon") == "p").filter(pl.col("clade").str.contains("binchicken")),
            coassembly_metadata.select("samples").unique(),
            iterations=100,
            join="samples",
            )
        .write_csv(args.output + "_coassembly_bacteria_phyla.tsv", separator="\t")
    )

    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").filter(pl.col("taxon") == "p").filter(pl.col("clade").str.contains("binchicken")),
            coassembly_metadata.select("samples").unique(),
            iterations=100,
            join="samples",
            )
        .write_csv(args.output + "_coassembly_archaea_phyla.tsv", separator="\t")
    )

    # Classes
    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").filter(pl.col("taxon") == "c").filter(pl.col("clade").str.contains("binchicken")),
            coassembly_metadata.select("samples").unique(),
            iterations=100,
            join="samples",
            )
        .write_csv(args.output + "_coassembly_bacteria_classes.tsv", separator="\t")
    )

    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").filter(pl.col("taxon") == "c").filter(pl.col("clade").str.contains("binchicken")),
            coassembly_metadata.select("samples").unique(),
            iterations=100,
            join="samples",
            )
        .write_csv(args.output + "_coassembly_archaea_classes.tsv", separator="\t")
    )

    # Orders
    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").filter(pl.col("taxon") == "o").filter(pl.col("clade").str.contains("binchicken")),
            coassembly_metadata.select("samples").unique(),
            iterations=100,
            join="samples",
            )
        .write_csv(args.output + "_coassembly_bacteria_orders.tsv", separator="\t")
    )

    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").filter(pl.col("taxon") == "o").filter(pl.col("clade").str.contains("binchicken")),
            coassembly_metadata.select("samples").unique(),
            iterations=100,
            join="samples",
            )
        .write_csv(args.output + "_coassembly_archaea_orders.tsv", separator="\t")
    )

    # Families
    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").filter(pl.col("taxon") == "f").filter(pl.col("clade").str.contains("binchicken")),
            coassembly_metadata.select("samples").unique(),
            join="samples",
            )
        .write_csv(args.output + "_coassembly_bacteria_families.tsv", separator="\t")
    )

    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").filter(pl.col("taxon") == "f").filter(pl.col("clade").str.contains("binchicken")),
            coassembly_metadata.select("samples").unique(),
            join="samples",
            )
        .write_csv(args.output + "_coassembly_archaea_families.tsv", separator="\t")
    )

    # Genera
    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").filter(pl.col("taxon") == "g").filter(pl.col("clade").str.contains("binchicken")),
            coassembly_metadata.select("samples").unique(),
            join="samples",
            )
        .write_csv(args.output + "_coassembly_bacteria_genera.tsv", separator="\t")
    )

    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").filter(pl.col("taxon") == "g").filter(pl.col("clade").str.contains("binchicken")),
            coassembly_metadata.select("samples").unique(),
            join="samples",
            )
        .write_csv(args.output + "_coassembly_archaea_genera.tsv", separator="\t")
    )

    # Species
    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").filter(pl.col("taxon") == "s").filter(pl.col("clade").str.contains(pl.col("name"))),
            coassembly_metadata.select("samples").unique(),
            join="samples",
            )
        .write_csv(args.output + "_coassembly_bacteria_species.tsv", separator="\t")
    )

    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").filter(pl.col("taxon") == "s").filter(pl.col("clade").str.contains(pl.col("name"))),
            coassembly_metadata.select("samples").unique(),
            join="samples",
            )
        .write_csv(args.output + "_coassembly_archaea_species.tsv", separator="\t")
    )

    ###############
    ### Species ###
    ###############
    # Phyla
    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").filter(pl.col("taxon") == "p").filter(pl.col("clade").str.contains("binchicken")),
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").select("name").unique(),
            iterations=100
            )
        .write_csv(args.output + "_bacteria_phyla.tsv", separator="\t")
    )

    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").filter(pl.col("taxon") == "p").filter(pl.col("clade").str.contains("binchicken")),
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").select("name").unique(),
            iterations=100
            )
        .write_csv(args.output + "_archaea_phyla.tsv", separator="\t")
    )

    # Classes
    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").filter(pl.col("taxon") == "c").filter(pl.col("clade").str.contains("binchicken")),
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").select("name").unique(),
            iterations=100
            )
        .write_csv(args.output + "_bacteria_classes.tsv", separator="\t")
    )

    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").filter(pl.col("taxon") == "c").filter(pl.col("clade").str.contains("binchicken")),
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").select("name").unique(),
            iterations=100
            )
        .write_csv(args.output + "_archaea_classes.tsv", separator="\t")
    )

    # Orders
    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").filter(pl.col("taxon") == "o").filter(pl.col("clade").str.contains("binchicken")),
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").select("name").unique(),
            iterations=100
            )
        .write_csv(args.output + "_bacteria_orders.tsv", separator="\t")
    )

    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").filter(pl.col("taxon") == "o").filter(pl.col("clade").str.contains("binchicken")),
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").select("name").unique(),
            iterations=100
            )
        .write_csv(args.output + "_archaea_orders.tsv", separator="\t")
    )

    # Families
    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").filter(pl.col("taxon") == "f").filter(pl.col("clade").str.contains("binchicken")),
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").select("name").unique()
            )
        .write_csv(args.output + "_bacteria_families.tsv", separator="\t")
    )

    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").filter(pl.col("taxon") == "f").filter(pl.col("clade").str.contains("binchicken")),
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").select("name").unique()
            )
        .write_csv(args.output + "_archaea_families.tsv", separator="\t")
    )

    # Genera
    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").filter(pl.col("taxon") == "g").filter(pl.col("clade").str.contains("binchicken")),
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").select("name").unique()
            )
        .write_csv(args.output + "_bacteria_genera.tsv", separator="\t")
    )

    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").filter(pl.col("taxon") == "g").filter(pl.col("clade").str.contains("binchicken")),
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").select("name").unique()
            )
        .write_csv(args.output + "_archaea_genera.tsv", separator="\t")
    )

    # Species
    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").filter(pl.col("taxon") == "s").filter(pl.col("clade").str.contains(pl.col("name"))),
            genome_taxonomy.filter(pl.col("domain") == "d__Bacteria").select("name").unique()
            )
        .write_csv(args.output + "_bacteria_species.tsv", separator="\t")
    )

    (
        rarefy(
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").filter(pl.col("taxon") == "s").filter(pl.col("clade").str.contains(pl.col("name"))),
            genome_taxonomy.filter(pl.col("domain") == "d__Archaea").select("name").unique()
            )
        .write_csv(args.output + "_archaea_species.tsv", separator="\t")
    )

    logging.info("Done")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
