#!/usr/bin/env python3

"""
Author: Samuel Aroney
Filter Eggnog mapper results
"""

import os
import sys
import argparse
import logging
import polars as pl

EGGNOG_COLUMNS={
    "query": str,
    "seed_ortholog": str,
    "evalue": float,
    "score": float,
    "eggNOG_OGs": str,
    "max_annot_lvl": str,
    "COG_category": str,
    "Description": str,
    "Preferred_name": str,
    "GOs": str,
    "EC": str,
    "KEGG_ko": str,
    "KEGG_Pathway": str,
    "KEGG_Module": str,
    "KEGG_Reaction": str,
    "KEGG_rclass": str,
    "BRITE": str,
    "KEGG_TC": str,
    "CAZy": str,
    "BiGG_Reaction": str,
    "PFAMs": str,
}

kos_of_interest = [
    "K11180", # (r)sulfite reductase dsrA
    "K17224", # Sulfur oxidation soxB
    "K00370", # Nitrate (oxido)reductase narG/nxrA
    "K02567", # Nitrate reductase napA
    "K00368", # Nitrite reductase nirK
    "K15864", # Nitrite reductase nirS
    "K10944", # Ammonia/methane monooxygenase amoA/pmoA
    "K16157", # Methane monooxygenase mmoX
    "K00399", # Methanogenesis mcrA
    "K02556", # Flagella motA
    "K02410", # Flagella fliG
    "K03301", # Parasitism ATP/ADP translocase
    "K02588", # Nitrogen fixation nifH
    "K22899", # Nitrogen fixation vnfH
    "K01601", # Calvin-Benson RubisCO
    "K15230", # Reductive TCA aclA
    "K14138", # Wood-Ljungdahl acsB
    "K14468", # Fuchs-Holo mcr
    "K02437", # Reductive glycine gcvH
    "K02689", # Oxygenic photosynthesis psaA
    "K02703", # Oxygenic photosynthesis psbA
    "K08928", # Anoxygenic photosynthesis (purple bacteria) pufL
    "K08929", # Anoxygenic photosynthesis (purple bacteria) pufM
    "K08944", # Anoxygenic photosynthesis (green sulfur) fmoA
]

kos_of_key_genes = [
    "K01512", # Acetate fermentation
    "K00925", # Acetate fermentation
    "K13953", # Ethanol fermentation
    "K04072", # Ethanol fermentation
    "K00016", # Lactate fermentation
    "K00101", # Lactate fermentation
    "K00102", # Lactate fermentation
    "K03777", # Lactate fermentation
    "K03778", # Lactate fermentation
    "K11264", # Propionate fermentation
    "K03416", # Propionate fermentation
    "K01026", # Propionate fermentation
    "K00634", # Butanoate fermentation
    "K00929", # Butanoate fermentation
    "K00370", # Nitrate reduction
    "K00371", # Nitrate reduction
    "K00374", # Nitrate reduction
    "K02567", # Nitrate reduction
    "K02568", # Nitrate reduction
    "K00362", # Nitrite reduction
    "K00363", # Nitrite reduction
    "K03385", # Nitrite reduction
    "K15876", # Nitrite reduction
    "K00368", # Nitrite reduction
    "K15864", # Nitrite reduction
    "K20932", # Anammox
    "K20933", # Anammox
    "K20934", # Anammox
    "K20935", # Anammox
    "K00958", # Sulfate reduction
    "K00955", # Sulfate reduction
    "K00956", # Sulfate reduction
    "K00957", # Sulfate reduction
    "K11180", # Sulfite reduction
    "K11181", # Sulfite reduction
    "K00399", # Methanogenesis/ANME
    "K00401", # Methanogenesis/ANME
    "K00402", # Methanogenesis/ANME
    "K13255", # Fe reduction
    "K07309", # Se reduction
    "K07310", # Se reduction
    "K07311", # Se reduction
    "K07312", # Se reduction
    "K12527", # Se reduction
    "K12529", # Se reduction
    "K12528", # Se reduction
    "K17050", # Se reduction
    "K17051", # Se reduction
    "K17052", # Se reduction
]

kos_of_hydrogenases = [
    "K00436", # Hydrogenase
    "K00437", # Hydrogenase
    "K18008", # Hydrogenase
    "K00440", # Hydrogenase
    "K00441", # Hydrogenase
    "K00442", # Hydrogenase
    "K00443", # Hydrogenase
    "K05586", # Hydrogenase
    "K05587", # Hydrogenase
    "K05588", # Hydrogenase
    "K05922", # Hydrogenase
    "K05927", # Hydrogenase
    "K06281", # Hydrogenase
    "K06282", # Hydrogenase
    "K13942", # Hydrogenase
    "K14068", # Hydrogenase
    "K14069", # Hydrogenase
    "K14070", # Hydrogenase
    "K14086", # Hydrogenase
    "K14087", # Hydrogenase
    "K14088", # Hydrogenase
    "K14089", # Hydrogenase
    "K14090", # Hydrogenase
    "K14091", # Hydrogenase
    "K14092", # Hydrogenase
    "K14093", # Hydrogenase
    "K14094", # Hydrogenase
    "K14095", # Hydrogenase
    "K14096", # Hydrogenase
    "K14097", # Hydrogenase
    "K14098", # Hydrogenase
    "K14099", # Hydrogenase
    "K14100", # Hydrogenase
    "K14101", # Hydrogenase
    "K14102", # Hydrogenase
    "K14103", # Hydrogenase
    "K14104", # Hydrogenase
    "K14105", # Hydrogenase
    "K14106", # Hydrogenase
    "K14107", # Hydrogenase
    "K14108", # Hydrogenase
    "K14109", # Hydrogenase
    "K14110", # Hydrogenase
    "K14111", # Hydrogenase
    "K14112", # Hydrogenase
    "K14113", # Hydrogenase
    "K14114", # Hydrogenase
    "K14115", # Hydrogenase
    "K14116", # Hydrogenase
    "K14117", # Hydrogenase
    "K14118", # Hydrogenase
    "K14119", # Hydrogenase
    "K14120", # Hydrogenase
    "K14121", # Hydrogenase
    "K14122", # Hydrogenase
    "K14123", # Hydrogenase
    "K14124", # Hydrogenase
    "K14125", # Hydrogenase
    "K06862", # Hydrogenase
    "K17992", # Hydrogenase
    "K17993", # Hydrogenase
    "K17994", # Hydrogenase
    "K17995", # Hydrogenase
    "K17996", # Hydrogenase
    "K17997", # Hydrogenase
    "K17998", # Hydrogenase
    "K17999", # Hydrogenase
    "K18005", # Hydrogenase
    "K18006", # Hydrogenase
    "K18007", # Hydrogenase
    "K18016", # Hydrogenase
    "K18017", # Hydrogenase
    "K18023", # Hydrogenase
    "K18330", # Hydrogenase
    "K18331", # Hydrogenase
    "K18332", # Hydrogenase
    "K23548", # Hydrogenase
    "K23549", # Hydrogenase
]

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument("--annotations-list", help="List of eggnog mapper annotations files, newline separated", required=True)
    parser.add_argument("--cxxch-results", help="CXXCH.py output", required=True)
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

    with open(args.annotations_list) as f:
        annotations_files = [l.strip() for l in f]

    logging.info(f"Loading {len(annotations_files)} annotations files")
    queries = []
    ko_regex = "|".join(kos_of_interest + kos_of_key_genes + kos_of_hydrogenases)
    pfam_regex = "|".join(["Cytochrome_C", "Cytochrome_c", "Multi-haem_cyto", "Cytochrom_C", "Cytochrom_c"])
    for file in annotations_files:
        queries.append(
            pl.scan_csv(
                file,
                separator="\t",
                has_header=False,
                new_columns=list(EGGNOG_COLUMNS.keys()),
                dtypes=EGGNOG_COLUMNS,
                comment_char="#",
                )
            .with_columns(genome = pl.lit(file).str.extract(r"([^/]+).emapper.annotations").str.replace("_protein.faa", ""))
            .filter((pl.col("KEGG_ko").str.contains(ko_regex)) | (pl.col("PFAMs").str.contains(pfam_regex)))
        )

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    pl.concat(queries).sink_csv(args.output, separator="\t")

    logging.info("Filtering CXXCH genes")
    cxxch = (
        pl.read_csv(args.cxxch_results, separator="\t")
        .select(
            "genome",
            query = "contig",
            CXXCH = pl.col("CXXCH") + pl.col("CXXXCH")
            )
        .lazy()
    )

    (
        pl.scan_csv(args.output, separator="\t")
        .join(cxxch, on=["genome", "query"])
        .filter(pl.col("PFAMs").str.contains(pfam_regex))
        .select(
            "genome", "query", "KEGG_ko", "PFAMs", "CXXCH",
            MHC = pl.lit("MHC")
        )
        .sink_csv(args.output.replace(".tsv", "_MHCs.tsv"), separator="\t")
    )

    logging.info("Filtering hydrogenases")
    hydrogenases_regex = "|".join(kos_of_hydrogenases)
    (
        pl.scan_csv(args.output, separator="\t")
        .filter(pl.col("KEGG_ko").str.contains(hydrogenases_regex))
        .sink_csv(args.output.replace(".tsv", "_hydrogenases.tsv"), separator="\t")
    )

    logging.info("Done")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
