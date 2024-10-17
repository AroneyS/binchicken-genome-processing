#!/usr/bin/env python3

"""
Author: Samuel Aroney
Assign TEAs to genomes from eggnog filtered outputs
"""

import sys
import argparse
import logging
import polars as pl

pathway_key_genes = pl.DataFrame([
    ["Acetate_fermentation_1", "K01512", 1],
    ["Acetate_fermentation_2", "K00925", 1],
    ["Ethanol_fermentation_1", "K13953", 1],
    ["Ethanol_fermentation_2", "K04072", 1],
    ["Lactate_fermentation_1", "K00016", 1],
    ["Lactate_fermentation_2", "K00101", 1],
    ["Lactate_fermentation_3", "K00102", 1],
    ["Lactate_fermentation_4", "K03777", 1],
    ["Lactate_fermentation_5", "K03778", 1],
    ["Propionate_fermentation_1", "K11264", 2],
    ["Propionate_fermentation_1", "K01026", 2],
    ["Propionate_fermentation_2", "K03416", 2],
    ["Propionate_fermentation_2", "K01026", 2],
    ["Butanoate_fermentation_1", "K00634", 2],
    ["Butanoate_fermentation_1", "K00929", 2],
    ["Nitrate_reduction_1", "K00370", 3],
    ["Nitrate_reduction_1", "K00371", 3],
    ["Nitrate_reduction_1", "K00374", 3],
    ["Nitrate_reduction_2", "K02567", 2],
    ["Nitrate_reduction_2", "K02568", 2],
    ["Nitrite_reduction_1", "K00362", 2],
    ["Nitrite_reduction_1", "K00363", 2],
    ["Nitrite_reduction_2", "K03385", 2],
    ["Nitrite_reduction_2", "K15876", 2],
    ["Nitrite_reduction_3", "K00368", 1],
    ["Nitrite_reduction_4", "K15864", 1],
    ["Anammox_1", "K20932", 4],
    ["Anammox_1", "K20933", 4],
    ["Anammox_1", "K20934", 4],
    ["Anammox_1", "K20935", 4],
    ["Sulfate_reduction_1", "K00958", 1],
    ["Sulfate_reduction_2", "K00955", 1],
    ["Sulfate_reduction_3", "K00956", 2],
    ["Sulfate_reduction_3", "K00957", 2],
    ["Sulfite_reduction_1", "K11180", 2],
    ["Sulfite_reduction_1", "K11181", 2],
    ["Methanogenesis_1", "K00399", 3],
    ["Methanogenesis_1", "K00401", 3],
    ["Methanogenesis_1", "K00402", 3],
    ["Fe_reduction_1", "K13255", 1],
    ["Se_reduction_1", "K07309", 4],
    ["Se_reduction_1", "K07310", 4],
    ["Se_reduction_1", "K07311", 4],
    ["Se_reduction_1", "K07312", 4],
    ["Se_reduction_2", "K12527", 3],
    ["Se_reduction_2", "K12529", 3],
    ["Se_reduction_2", "K12528", 3],
    ["Se_reduction_3", "K17050", 3],
    ["Se_reduction_3", "K17051", 3],
    ["Se_reduction_3", "K17052", 3],
], schema={"pathway": str, "KEGG_ko": str, "key_gene_count": int})

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument("--comb-results", help="Output from filter_eggnog", required=True)
    parser.add_argument("--mhc-results", help="Output from filter_eggnog", required=True)
    parser.add_argument("--aer-results", help="Output from aerobicity", required=True)
    parser.add_argument("--output", help="Output file")

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.info("Loading files")
    comb = pl.concat([
        pl.read_csv(args.comb_results, separator="\t")
            .with_columns(pl.col("KEGG_ko").str.split(","))
            .explode("KEGG_ko")
            .with_columns(pl.col("KEGG_ko").str.replace(r"^ko:", ""))
            .select("genome", "query", "KEGG_ko")
            .join(pathway_key_genes, on="KEGG_ko")
            .unique()
            .group_by("genome", "pathway", "key_gene_count")
            .count()
            .select(
                "genome", "pathway",
                evidence =
                    pl.when(pl.col("key_gene_count") == 1)
                    .then(pl.lit("single_gene"))
                    .when(pl.col("count") >= pl.col("key_gene_count"))
                    .then(pl.lit("full"))
                    .otherwise(pl.lit("partial"))
                ),
        pl.read_csv(args.mhc_results, separator="\t")
            .group_by("genome")
            .count()
            .select(
                "genome",
                pathway = pl.lit("MHC"),
                evidence = pl.when(pl.col("count") == 1).then(pl.lit("single")).otherwise(pl.lit("multiple"))
                ),
        pl.read_csv(args.aer_results, separator="\t")
            .filter(pl.col("prediction") == 1)
            .select(
                genome = pl.col("node").str.extract(r".*/([^/]+)(_transcript|_protein)"),
                pathway = pl.lit("aerobic"),
                evidence = pl.lit("full")
                ),
    ])

    comb.write_csv(args.output, separator="\t")

    logging.info("Done")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
