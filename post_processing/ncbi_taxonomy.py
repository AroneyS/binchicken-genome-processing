#!/usr/bin/env python3

"""
Author: Samuel Aroney
Process NCBI taxonomy into table

python scripts/ncbi_taxonomy.py \
    --sra-taxonomy-dir /home/aroneys/m/big_data_microbiome/14_sra_metadata_20240520/sra_taxonomy_table_20240520 \
    --output data/ncbi_taxonomy.tsv
"""

import os
import sys
import argparse
import logging
import polars as pl
import json

class Taxonomy:
    def __str__(self):
        return '{} {} {} {}'.format(self.tax_id, self.xlevel, self.parent_id, self.sci_name)
    def __repr__(self):
        return str(self)

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument('--sra-taxonomy-dir', help='sra_taxonomy_table.json files folder', required=True)
    parser.add_argument("--output", help="Output file", required=True)

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.info("Loading taxonomy JSON ..")
    taxonomies = {}
    filenames = ["000000000000", "000000000002", "000000000001"]
    for filename in filenames:
        with open(os.path.join(args.sra_taxonomy_dir, filename)) as f:
            for line in f:
                j = json.loads(line)
                # taxonomy = Taxonomy()
                # taxonomy.sci_name = j['sci_name']
                # taxonomy.parent_id = int(j['parent_id'])
                # taxonomy.xlevel = int(j['xlevel'])
                # taxonomy.tax_id = int(j['tax_id'])
                taxonomy = (j['sci_name'], int(j['parent_id']), int(j['xlevel']), int(j['tax_id']))
                try:
                    taxonomies[filename].append(taxonomy)
                except KeyError:
                    taxonomies[filename] = [taxonomy]

    taxonomy_df = (
        pl.DataFrame(
            taxonomies[filenames[0]],
            schema={"sci_name": str, "parent_id": int, "xlevel": int, "tax_id": int}
            )
        .join(
            pl.DataFrame(
                taxonomies[filenames[1]],
                schema={"sci_name": str, "parent_id": int, "xlevel": int, "tax_id": int}
                )
            .select(
                parent_id = pl.col("tax_id"),
                parent_name = pl.col("sci_name"),
                superparent_id = pl.col("parent_id")
                ),
            on = "parent_id",
            )
        .join(
            pl.DataFrame(
                taxonomies[filenames[2]],
                schema={"sci_name": str, "parent_id": int, "xlevel": int, "tax_id": int}
                )
            .select(
                superparent_id = pl.col("tax_id"),
                superparent_name = pl.col("sci_name"),
                ),
            on = "superparent_id",
            )
    )

    (
        taxonomy_df
        .select("sci_name", "parent_name", "superparent_name")
        .write_csv(args.output, separator="\t")
    )
    logging.info("Done")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
