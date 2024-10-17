#!/usr/bin/env python3

import unittest
import polars as pl
from polars.testing import assert_frame_equal
from post_processing import get_coassemblies, name_coassemblies, get_coassembly_info, get_genomes, name_genomes, get_genome_info, get_genome_novelty, get_genome_taxa_recovery, update_metadata

GET_COASSEMBLIES_INPUT_COLUMNS={
    "binchicken_dir": str,
    "style": str,
    "scope": str,
    "target": str,
    "appraise_level": str,
    "aviary_version": str,
    }
GET_COASSEMBLIES_OUTPUT_COLUMNS=GET_COASSEMBLIES_INPUT_COLUMNS | {
    "assemble_dir": str,
    "recover_dir": str,
    }

NAME_COASSEMBLIES_OLD_COLUMNS={
    "name": str,
    }
NAME_COASSEMBLIES_OUTPUT_COLUMNS={
    "name": str,
    } | GET_COASSEMBLIES_OUTPUT_COLUMNS

GET_COASSEMBLY_INFO_OUTPUT_COLUMNS = NAME_COASSEMBLIES_OUTPUT_COLUMNS | {
    "binchicken_version": str,
    "original_name": str,
    "samples": str,
    "recover_samples": str,
    "unmapping_identity": float,
    "unmapping_alignment": float,
    "assembler": str,
    "assembly_size": int,
    "contigs_total": int,
    "contig_N50": int,
    "gc": float,
}

GET_GENOMES_COLUMNS = {
    "name": str,
    "recover_dir": str,
}

GET_GENOMES_OUTPUT_COLUMNS = {
    "coassembly": str,
    "recover_dir": str,
    "genome_path": str,
}

NAME_GENOMES_OUTPUT_COLUMNS={
    "name": str,
    } | GET_GENOMES_OUTPUT_COLUMNS

GET_GENOME_INFO_OUTPUT_COLUMNS = {
    "name": str,
    "original_name": str,
    } | GET_GENOMES_OUTPUT_COLUMNS | {
    "binner": str,
    "genome_size": int,
    "contigs_total": int,
    "contig_N50": int,
    "contig_max": int,
    "gc": float,
    "coding_density": float,
    "translation_table": int,
    "completeness": float,
    "contamination": float,
    "CheckM_model": str,
    "CheckM_version": str,
    "taxonomy": str,
    "red_value": float,
    "GTDBtk_method": str,
    "GTDBtk_version": str,
    "GTDBtk_note": str,
    "GTDBtk_warnings": str,
}

GET_GENOME_NOVELTY_INPUT_COLUMNS = {
    "name": str,
    "completeness": float,
    "contamination": float,
    "95ANI_representative": str,
    "taxonomy": str,
}

GET_GENOME_NOVELTY_OUTPUT_COLUMNS = {
    "domain": str,
    "novelty": str,
    "distinct": int,
    "total": int,
}

GET_GENOME_TAXA_RECOVERY_INPUT_COLUMNS = {
    "name": str,
    "target": str,
}

GET_GENOME_TAXA_RECOVERY_OUTPUT_COLUMNS = {
    "target": str,
    "domain": str,
    "phyla": int,
    "class": int,
    "order": int,
    "family": int,
    "genus": int,
    "species": int,
    "strain": int,
}

CLUSTER_COLUMNS = {
    "name": str,
    "95ANI_representative": str,
}

GTDBTK_COLUMNS = {
    "user_genome": str,
    "classification": str,
    "red_value": float,
    "classification_method": str,
    "note": str,
    "warnings": str,
}

GENE_COUNTS_COLUMNS = {
    "name": str,
    "gene_count": int,
}

GUNC_COLUMNS = {
    "genome": str,
    "taxonomic_level": str,
    "contamination_portion": float,
    "reference_representation_score": float,
    "pass.GUNC": bool,
}

RRNA_COLUMNS = {
    "name": str,
    "seqname": str,
    "source": str,
    "feature": str,
    "start": str,
    "end": str,
    "score": str,
    "strand": str,
    "frame": str,
    "attribute": str,
}

TRNA_COLUMNS = {
    "name": str,
    "seqname": str,
    "tRNA_no": str,
    "begin": str,
    "end": str,
    "type": str,
    "codon": str,
    "intron_begin": str,
    "intron_end": str,
    "score": str,
    "note": str,
}

UPDATE_METADATA_OUTPUT_COLUMNS = {
    "name": str,
    "original_name": str,
    "coassembly": str,
    "genome_path": str,
    "binner": str,
    "genome_size": int,
    "contigs_total": int,
    "contig_N50": int,
    "contig_max": int,
    "gc": float,
    "coding_density": float,
    "translation_table": int,
    "completeness": float,
    "contamination": float,
    "CheckM_model": str,
    "CheckM_version": str,
    "95ANI_representative": str,
    "taxonomy": str,
    "red_value": float,
    "GTDBtk_method": str,
    "GTDBtk_version": str,
    "GTDBtk_note": str,
    "GTDBtk_warnings": str,
    "gene_count": int,
    "GUNC_tax_level": str,
    "GUNC_contamination": float,
    "GUNC_representation": float,
    "GUNC_pass": bool,
    "rRNA_5S": int,
    "rRNA_16S": int,
    "rRNA_23S": int,
    "tRNAs": int,
}

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtype=False, check_row_order=False)

    @unittest.skip("Requires actual output files to test")
    def test_get_coassemblies(self):
        df = pl.DataFrame([
            ["results/aviary/target/20230619/p__Abyssobacteria", "taxa-targeting", "global", "p__Abyssobacteria", "0.86", "v0.6.0"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0"],
        ], schema=GET_COASSEMBLIES_INPUT_COLUMNS)

        expected_output = pl.DataFrame([
            ["results/aviary/target/20230619/p__Abyssobacteria", "taxa-targeting", "global", "p__Abyssobacteria", "0.86", "v0.6.0", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/assemble", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/recover"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/recover"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover"],
        ], schema=GET_COASSEMBLIES_OUTPUT_COLUMNS)

        observed_output = get_coassemblies(df)
        self.assertDataFrameEqual(expected_output, observed_output)

    def test_name_coassemblies(self):
        old = pl.DataFrame([
            ["SRP110325_co1"],
            ["SRP110325_co2"],
            ["SRP110325_co3"],
            ["SRP110325_co4"],
            ["SRP110325_co5"],
            ["SRP110325_co6"],
        ], schema = NAME_COASSEMBLIES_OLD_COLUMNS)

        new = pl.DataFrame([
            ["results/aviary/target/20230619/p__Abyssobacteria", "taxa-targeting", "global", "p__Abyssobacteria", "0.86", "v0.6.0", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/assemble", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/recover"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/recover"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover"],
        ], schema=GET_COASSEMBLIES_OUTPUT_COLUMNS)


        expected_output = pl.DataFrame([
            ["SRP110325_co7", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/recover"],
            ["SRP110325_co8", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover"],
            ["SRP110325_co9", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria", "taxa-targeting", "global", "p__Abyssobacteria", "0.86", "v0.6.0", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/assemble", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover"],
            ["SRP110325_co10", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/recover"],
        ], schema=NAME_COASSEMBLIES_OUTPUT_COLUMNS)

        observed_output = name_coassemblies(new, old)
        self.assertDataFrameEqual(expected_output, observed_output)

    def test_name_coassemblies_no_old(self):
        old = pl.DataFrame(schema = NAME_COASSEMBLIES_OLD_COLUMNS)

        new = pl.DataFrame([
            ["results/aviary/target/20230619/p__Abyssobacteria", "taxa-targeting", "global", "p__Abyssobacteria", "0.86", "v0.6.0", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/assemble", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/recover"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/recover"],
            ["results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover"],
        ], schema=GET_COASSEMBLIES_OUTPUT_COLUMNS)


        expected_output = pl.DataFrame([
            ["SRP110325_co1", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/recover"],
            ["SRP110325_co2", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover"],
            ["SRP110325_co3", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover"],
            ["SRP110325_co4", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/recover"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria", "taxa-targeting", "global", "p__Abyssobacteria", "0.86", "v0.6.0", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/assemble", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover"],
        ], schema=NAME_COASSEMBLIES_OUTPUT_COLUMNS)

        observed_output = name_coassemblies(new, old)
        self.assertDataFrameEqual(expected_output, observed_output)

    @unittest.skip("Requires actual output files to test")
    def test_get_coassembly_info(self):
        coassemblies = pl.DataFrame([
            ["SRP110325_co1", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/recover"],
            ["SRP110325_co2", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover"],
            ["SRP110325_co3", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/recover"],
            ["SRP110325_co4", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria", "taxa-targeting", "global", "p__Abyssobacteria", "0.86", "v0.6.0", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/assemble", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover"],
        ], schema=NAME_COASSEMBLIES_OUTPUT_COLUMNS)


        expected_output = pl.DataFrame([
            [
                "SRP110325_co1", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0",
                "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/recover",
                "unknown", "coassembly_1", "SRR5753876,SRR5753878", "SRR5753868,SRR5753875,SRR5753876,SRR5753877,SRR5753878,SRR5753880",
                95, 90, "metaspades", 8653744625, 20021866, 5107389, 0.3542,
            ],
            [
                "SRP110325_co2", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0",
                "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover",
                "unknown", "coassembly_1", "SRR5753876,SRR5753878", "SRR5753868,SRR5753875,SRR5753876,SRR5753877,SRR5753878,SRR5753880",
                95, 90, "metaspades", 8653744625, 20021866, 5107389, 0.3542,
            ],
            [
                "SRP110325_co3", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0",
                "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/recover",
                "unknown", "coassembly_0", "SRR5753868,SRR5753874", "SRR5753867,SRR5753868,SRR5753872,SRR5753873,SRR5753874,SRR5753876,SRR5753878",
                95, 90, "metaspades", 7251820768, 15780132, 3560191, 0.3558,
            ],
            [
                "SRP110325_co4", "results/aviary/phylum_2sample/20230718/SRP110325", "project-specific", "SRP110325", "", "0.64", "v0.6.0",
                "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/assemble", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover",
                "unknown", "coassembly_0", "SRR5753868,SRR5753874", "SRR5753867,SRR5753868,SRR5753872,SRR5753873,SRR5753874,SRR5753876,SRR5753878",
                95, 90, "metaspades", 7251820768, 15780132, 3560191, 0.3558,
            ],
            [
                "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria", "taxa-targeting", "global", "p__Abyssobacteria", "0.86", "v0.6.0",
                "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/assemble", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover",
                "unknown", "coassembly_1", "SRR7048218,SRR7048226", "SRR11241200,SRR14699966,SRR14699972,SRR14699978,SRR14699980,SRR14699982,SRR14699983,SRR14699986,SRR6963352,SRR6963480,SRR6963587,SRR6963591,SRR7048218,SRR7048219,SRR7048220,SRR7048223,SRR7048224,SRR7048226,SRR7048227,SRR7615342",
                95, 90, "metaspades", 1943949021, 3202261, 272183, 0.5362,
            ],
        ], schema=GET_COASSEMBLY_INFO_OUTPUT_COLUMNS)

        observed_output = get_coassembly_info(coassemblies)
        self.assertDataFrameEqual(expected_output, observed_output)

    @unittest.skip("Requires actual output files to test")
    def test_get_genomes(self):
        df = pl.DataFrame([
            ["SRP110325_co1", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1/recover"],
            ["SRP110325_co2", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover"],
            ["SRP110325_co3", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0/recover"],
            ["SRP110325_co4", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover"],
        ], schema=GET_GENOMES_COLUMNS)


        expected_output = pl.DataFrame([
            ["SRP110325_co2", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover/bins/final_bins/bin.1891.fna"],
            ["SRP110325_co4", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover/bins/final_bins/bin.18.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/concoct_bins.tsv.183_sub.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/metabat2_refined_bins.tsv.180_sub.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/maxbin_bins.tsv.187_sub.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/semibin_refined_bins.tsv.186_sub.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/maxbin_bins.tsv.188_sub.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/concoct_bins.tsv.184_sub.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/metabat2_refined_bins.tsv.187.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/vamb_bins.tsv.182_sub.fna"],
        ], schema=GET_GENOMES_OUTPUT_COLUMNS)

        observed_output = get_genomes(df).filter(pl.col("genome_path").str.contains(".18", literal=True))
        self.assertDataFrameEqual(expected_output, observed_output)

    def test_name_genomes(self):
        df = pl.DataFrame([
            ["SRP110325_co2", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover/bins/final_bins/bin.1891.fna"],
            ["SRP110325_co4", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover/bins/final_bins/bin.18.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/concoct_bins.tsv.183_sub.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/metabat2_refined_bins.tsv.180_sub.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/maxbin_bins.tsv.187_sub.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/semibin_refined_bins.tsv.186_sub.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/maxbin_bins.tsv.188_sub.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/concoct_bins.tsv.184_sub.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/metabat2_refined_bins.tsv.187.fna"],
            ["binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/vamb_bins.tsv.182_sub.fna"],
        ], schema=GET_GENOMES_OUTPUT_COLUMNS)


        expected_output = pl.DataFrame([
            ["SRP110325_co2_1", "SRP110325_co2", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover/bins/final_bins/bin.1891.fna"],
            ["SRP110325_co4_1", "SRP110325_co4", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover/bins/final_bins/bin.18.fna"],
            ["binchicken_co1_1", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/concoct_bins.tsv.183_sub.fna"],
            ["binchicken_co1_2", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/concoct_bins.tsv.184_sub.fna"],
            ["binchicken_co1_3", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/maxbin_bins.tsv.187_sub.fna"],
            ["binchicken_co1_4", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/maxbin_bins.tsv.188_sub.fna"],
            ["binchicken_co1_5", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/metabat2_refined_bins.tsv.180_sub.fna"],
            ["binchicken_co1_6", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/metabat2_refined_bins.tsv.187.fna"],
            ["binchicken_co1_7", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/semibin_refined_bins.tsv.186_sub.fna"],
            ["binchicken_co1_8", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/vamb_bins.tsv.182_sub.fna"],
        ], schema=NAME_GENOMES_OUTPUT_COLUMNS)

        observed_output = name_genomes(df)
        self.assertDataFrameEqual(expected_output, observed_output)

    @unittest.skip("Requires actual output files to test")
    def test_get_genome_info(self):
        genomes = pl.DataFrame([
            ["SRP110325_co2_1", "SRP110325_co2", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover/bins/final_bins/bin.1891.fna"],
            ["SRP110325_co4_1", "SRP110325_co4", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover", "results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover/bins/final_bins/bin.18.fna"],
            ["binchicken_co1_1", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/concoct_bins.tsv.183_sub.fna"],
            ["binchicken_co1_2", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/concoct_bins.tsv.184_sub.fna"],
            ["binchicken_co1_3", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/maxbin_bins.tsv.187_sub.fna"],
            ["binchicken_co1_4", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/maxbin_bins.tsv.188_sub.fna"],
            ["binchicken_co1_5", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/metabat2_refined_bins.tsv.180_sub.fna"],
            ["binchicken_co1_6", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/metabat2_refined_bins.tsv.187.fna"],
            ["binchicken_co1_7", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/semibin_refined_bins.tsv.186_sub.fna"],
            ["binchicken_co1_8", "binchicken_co1", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover", "results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/vamb_bins.tsv.182_sub.fna"],
        ], schema=NAME_GENOMES_OUTPUT_COLUMNS)


        expected_output = pl.DataFrame([
            [
                'SRP110325_co2_1', 'bin.1891', 'SRP110325_co2', 'results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover',
                'results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover/bins/final_bins/bin.1891.fna', None,
                729355, 140, 5592, 15740, 0.35, 0.205, 11, 4.52, 0.0,
                'Neural Network (Specific Model)', 'CheckM2v1.0.2',
                'Unclassified', None, None, 'gtdbtk-2.3.0', None, 'No bacterial or archaeal marker'
            ],
            [
                'SRP110325_co4_1', 'bin.18', 'SRP110325_co4', 'results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover',
                'results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_0_semibin2/recover/bins/final_bins/bin.18.fna', None,
                4830055, 1231, 3985, 9776, 0.35, 0.327, 11, 63.79, 9.77,
                'Gradient Boost (General Model)', 'CheckM2v1.0.2',
                'Unclassified', None, None, 'gtdbtk-2.3.0', None, 'No bacterial or archaeal marker'
            ],
            [
                'binchicken_co1_1', 'concoct_bins.tsv.183_sub', 'binchicken_co1', 'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover',
                'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/concoct_bins.tsv.183_sub.fna', 'concoct',
                1685615, 661, 2610, 7076, 0.59, 0.83, 11, 58.74, 0.55,
                'Gradient Boost (General Model)', 'CheckM2v1.0.2',
                'd__Bacteria;p__Desulfobacterota;c__Desulfobaccia;o__Desulfobaccales;f__0-14-0-80-60-11;g__DTHB01;s__', 0.86156, 'taxonomic classification defined by topology and ANI', 'gtdbtk-2.3.0', 'classification based on placement in class-level tree', None
            ],
            [
                'binchicken_co1_2', 'concoct_bins.tsv.184_sub', 'binchicken_co1', 'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover',
                'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/concoct_bins.tsv.184_sub.fna', 'concoct',
                406793, 174, 2232, 6450, 0.5, 0.896, 11, 18.78, 0.82,
                'Neural Network (Specific Model)', 'CheckM2v1.0.2',
                'd__Bacteria;p__Patescibacteria;c__Gracilibacteria;o__UM-FILTER-43-11;f__UM-FILTER-43-11;g__;s__', 0.90369, 'taxonomic novelty determined using RED', 'gtdbtk-2.3.0', 'classification based on placement in class-level tree', None
            ],
            [
                'binchicken_co1_3', 'maxbin_bins.tsv.187_sub', 'binchicken_co1', 'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover',
                'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/maxbin_bins.tsv.187_sub.fna', 'maxbin',
                2806436, 1163, 2183, 55868, 0.43, 0.883, 11, 41.63, 1.26,
                'Neural Network (Specific Model)', 'CheckM2v1.0.2',
                'd__Bacteria;p__OLB16;c__OLB16;o__SURF-12;f__;g__;s__', 0.62201, 'taxonomic novelty determined using RED', 'gtdbtk-2.3.0', 'classification based on placement in class-level tree', None
            ],
            [
                'binchicken_co1_4', 'maxbin_bins.tsv.188_sub', 'binchicken_co1', 'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover',
                'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/maxbin_bins.tsv.188_sub.fna', 'maxbin',
                2441711, 777, 3445, 14427, 0.4, 0.912, 11, 67.64, 7.88,
                'Gradient Boost (General Model)', 'CheckM2v1.0.2',
                'd__Bacteria;p__Acidobacteriota;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__WJMT01;s__', 0.92929, 'taxonomic classification defined by topology and ANI', 'gtdbtk-2.3.0', 'classification based on placement in class-level tree', None
            ],
            [
                'binchicken_co1_5', 'metabat2_refined_bins.tsv.180_sub', 'binchicken_co1', 'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover',
                'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/metabat2_refined_bins.tsv.180_sub.fna', 'metabat2',
                3563088, 1074, 3625, 15877, 0.49, 0.914, 11, 69.66, 3.46,
                'Gradient Boost (General Model)', 'CheckM2v1.0.2',
                'd__Bacteria;p__Planctomycetota;c__GCA-2746535;o__;f__;g__;s__', 0.49065, 'taxonomic classification defined by topology and ANI', 'gtdbtk-2.3.0', 'classification based on placement in class-level tree', None
            ],
            [
                'binchicken_co1_6', 'metabat2_refined_bins.tsv.187', 'binchicken_co1', 'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover',
                'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/metabat2_refined_bins.tsv.187.fna', 'metabat2',
                3154926, 155, 39961, 156589, 0.52, 0.904, 11, 96.22, 4.06,
                'Gradient Boost (General Model)', 'CheckM2v1.0.2',
                'd__Archaea;p__Hadarchaeota;c__Hadarchaeia;o__;f__;g__;s__', 0.28352, 'taxonomic novelty determined using RED', 'gtdbtk-2.3.0', None, None
            ],
            [
                'binchicken_co1_7', 'semibin_refined_bins.tsv.186_sub', 'binchicken_co1', 'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover',
                'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/semibin_refined_bins.tsv.186_sub.fna', 'semibin',
                859692, 243, 3292, 13164, 0.42, 0.907, 11, 22.0, 0.72,
                'Neural Network (Specific Model)', 'CheckM2v1.0.2',
                'd__Bacteria;p__Chlamydiota;c__Chlamydiia;o__Chlamydiales;f__SM23-39;g__;s__', 0.78077, 'taxonomic novelty determined using RED', 'gtdbtk-2.3.0', 'classification based on placement in class-level tree', None
            ],
            [
                'binchicken_co1_8', 'vamb_bins.tsv.182_sub', 'binchicken_co1', 'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover',
                'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/vamb_bins.tsv.182_sub.fna', 'vamb',
                606959, 237, 2583, 19336, 0.46, 0.926, 11, 32.93, 0.23,
                'Neural Network (Specific Model)', 'CheckM2v1.0.2',
                'd__Archaea;p__Micrarchaeota;c__Micrarchaeia;o__Anstonellales;f__;g__;s__', 0.55037, 'taxonomic classification fully defined by topology', 'gtdbtk-2.3.0', None, None
            ],
        ], schema=GET_GENOME_INFO_OUTPUT_COLUMNS)

        observed_output = get_genome_info(genomes)
        self.assertDataFrameEqual(expected_output, observed_output)

    def test_get_genome_novelty(self):
        genome_info = pl.DataFrame([
            # original genomes
            ["binchicken_co1_1", 63.79, 1.77, "binchicken_co1_1", 'Unclassified'],
            ["binchicken_co1_2", 58.74, 0.55, "binchicken_co1_2", 'd__Bacteria;p__Desulfobacterota;c__Desulfobaccia;o__Desulfobaccales;f__0-14-0-80-60-11;g__DTHB01;s__'],
            ["binchicken_co1_3", 67.64, 2.88, "binchicken_co1_3", 'd__Bacteria;p__Acidobacteriota;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__WJMT01;s__'],
            ["binchicken_co1_4", 69.66, 2.46, "binchicken_co1_4", 'd__Bacteria;p__Planctomycetota;c__GCA-2746535;o__;f__;g__;s__'],
            ["binchicken_co1_5", 96.22, 4.06, "binchicken_co1_5", 'd__Archaea;p__Hadarchaeota;c__Hadarchaeia;o__;f__;g__;s__'],
            # mock genomes
            ["binchicken_co1_6", 90.11, 5.22, "binchicken_co1_6", 'Unclassified'],
            ["binchicken_co1_7", 90.11, 5.22, "binchicken_co1_7", 'd__Bacteria;p__;c__;o__;f__;g__;s__'],
            ["binchicken_co1_8", 90.11, 5.22, "binchicken_co1_8", 'd__Bacteria;p__;c__;o__;f__;g__;s__'],
            ["binchicken_co1_9", 90.11, 5.22, "binchicken_co1_9", 'd__Bacteria;p__Acidobacteriota;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__WJMT01;s__'],
            ["binchicken_co1_10", 90.11, 5.22, "binchicken_co1_10", 'd__Bacteria;p__Planctomycetota;c__GCA-2746535;o__;f__;g__;s__'],
            ["binchicken_co1_11", 90.11, 5.22, "binchicken_co1_11", 'd__Archaea;p__Hadarchaeota;c__Hadarchaeia2;o__;f__;g__;s__'],
            ["binchicken_co1_12", 90.11, 5.22, "binchicken_co1_12", 'd__Bacteria;p__Acidobacteriota;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__WJMT01;s__mock'],
            # Non-representatives
            ["binchicken_co1_18", 90.11, 5.22, "binchicken_co1_10", 'd__Bacteria;p__Planctomycetota;c__GCA-2746535;o__;f__;g__;s__'],
            ["binchicken_co1_19", 90.11, 5.22, "binchicken_co1_11", 'd__Archaea;p__Hadarchaeota;c__Hadarchaeia2;o__;f__;g__;s__'],
            ["binchicken_co1_20", 90.11, 5.22, "binchicken_co1_12", 'd__Bacteria;p__Acidobacteriota;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__WJMT01;s__mock'],
            # All fail quality
            ["binchicken_co1_13", 4.52, 0.0, "binchicken_co1_13", 'Unclassified'],
            ["binchicken_co1_14", 18.78, 0.82, "binchicken_co1_14", 'd__Bacteria;p__Patescibacteria;c__Gracilibacteria;o__UM-FILTER-43-11;f__UM-FILTER-43-11;g__;s__'],
            ["binchicken_co1_15", 41.63, 1.26, "binchicken_co1_15", 'd__Bacteria;p__OLB16;c__OLB16;o__SURF-12;f__;g__;s__'],
            ["binchicken_co1_16", 22.0, 0.72, "binchicken_co1_16", 'd__Bacteria;p__Chlamydiota;c__Chlamydiia;o__Chlamydiales;f__SM23-39;g__;s__'],
            ["binchicken_co1_17", 32.93, 0.23, "binchicken_co1_17", 'd__Archaea;p__Micrarchaeota;c__Micrarchaeia;o__Anstonellales;f__;g__;s__'],
        ], schema=GET_GENOME_NOVELTY_INPUT_COLUMNS)

        expected_output = pl.DataFrame([
            ["", "unclassified", 1, 2],
            ["d__Archaea", "phylum", 0, 0],
            ["d__Archaea", "class", 0, 0],
            ["d__Archaea", "order", 2, 2],
            ["d__Archaea", "family", 0, 0],
            ["d__Archaea", "genus", 0, 0],
            ["d__Archaea", "species", 0, 0],
            ["d__Archaea", "strain", 0, 0],
            ["d__Bacteria", "phylum", 1, 2],
            ["d__Bacteria", "class", 0, 0],
            ["d__Bacteria", "order", 1, 2],
            ["d__Bacteria", "family", 0, 0],
            ["d__Bacteria", "genus", 0, 0],
            ["d__Bacteria", "species", 2, 3],
            ["d__Bacteria", "strain", 1, 1],
        ], schema=GET_GENOME_NOVELTY_OUTPUT_COLUMNS)

        observed_output = get_genome_novelty(genome_info)
        self.assertDataFrameEqual(expected_output, observed_output)

    def test_get_genome_taxa_recovery(self):
        genome_info = pl.DataFrame([
            # original genomes
            ["binchicken_co1_1", 63.79, 9.77, "binchicken_co1_1", 'Unclassified'],
            ["binchicken_co1_2", 58.74, 0.55, "binchicken_co1_2", 'd__Bacteria;p__Desulfobacterota;c__Desulfobaccia;o__Desulfobaccales;f__0-14-0-80-60-11;g__DTHB01;s__'],
            ["binchicken_co1_3", 67.64, 7.88, "binchicken_co1_3", 'd__Bacteria;p__Acidobacteriota;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__WJMT01;s__'],
            ["binchicken_co1_4", 69.66, 3.46, "binchicken_co1_4", 'd__Bacteria;p__Planctomycetota;c__GCA-2746535;o__;f__;g__;s__'],
            ["binchicken_co1_5", 96.22, 4.06, "binchicken_co1_5", 'd__Archaea;p__Hadarchaeota;c__Hadarchaeia;o__;f__;g__;s__'],
            # mock genomes
            ["binchicken_co1_6", 90.11, 5.22, "binchicken_co1_6", 'Unclassified'],
            ["binchicken_co1_7", 90.11, 5.22, "binchicken_co1_7", 'd__Bacteria;p__Thing1;c__;o__;f__;g__;s__'],
            ["binchicken_co1_8", 90.11, 5.22, "binchicken_co1_8", 'd__Bacteria;p__Thing1;c__;o__;f__;g__;s__'],
            ["binchicken_co1_9", 90.11, 5.22, "binchicken_co1_9", 'd__Bacteria;p__Thing-2;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__WJMT01;s__'],
            ["binchicken_co1_10", 90.11, 5.22, "binchicken_co1_10", 'd__Bacteria;p__Thing-2;c__GCA-2746535;o__;f__;g__;s__'],
            ["binchicken_co1_11", 90.11, 5.22, "binchicken_co1_11", 'd__Archaea;p__Thing-2;c__Hadarchaeia2;o__;f__;g__;s__'],
            ["binchicken_co1_12", 90.11, 5.22, "binchicken_co1_12", 'd__Bacteria;p__Thing-2;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__WJMT01;s__mock'],
            # Non-representatives
            ["binchicken_co1_18", 90.11, 5.22, "binchicken_co1_10", 'd__Bacteria;p__Thing-2;c__GCA-2746535;o__;f__;g__;s__'],
            ["binchicken_co1_19", 90.11, 5.22, "binchicken_co1_11", 'd__Archaea;p__Thing-2;c__Hadarchaeia2;o__;f__;g__;s__'],
            ["binchicken_co1_20", 90.11, 5.22, "binchicken_co1_12", 'd__Bacteria;p__Thing-2;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__WJMT01;s__mock'],
            # All fail quality
            ["binchicken_co1_13", 4.52, 0.0, "binchicken_co1_13", 'Unclassified'],
            ["binchicken_co1_14", 18.78, 0.82, "binchicken_co1_14", 'd__Bacteria;p__Patescibacteria;c__Gracilibacteria;o__UM-FILTER-43-11;f__UM-FILTER-43-11;g__;s__'],
            ["binchicken_co1_15", 41.63, 1.26, "binchicken_co1_15", 'd__Bacteria;p__OLB16;c__OLB16;o__SURF-12;f__;g__;s__'],
            ["binchicken_co1_16", 22.0, 0.72, "binchicken_co1_16", 'd__Bacteria;p__Chlamydiota;c__Chlamydiia;o__Chlamydiales;f__SM23-39;g__;s__'],
            ["binchicken_co1_17", 32.93, 0.23, "binchicken_co1_17", 'd__Archaea;p__Micrarchaeota;c__Micrarchaeia;o__Anstonellales;f__;g__;s__'],
        ], schema=GET_GENOME_NOVELTY_INPUT_COLUMNS)

        coassemble_info = pl.DataFrame([
            ["Thing1_co1", "Thing1"],
            ["Thing2_co1", "Thing-2"],
            ["Thing3_co1", "Thing3"],
            ["OtherThing_co1", None],
        ], schema=GET_GENOME_TAXA_RECOVERY_INPUT_COLUMNS, orient="row")

        expected_output = pl.DataFrame([
            ["Thing1", "d__Bacteria",  0, 2, 0, 0, 0, 0, 0],
            ["Thing-2", "d__Archaea",  0, 0, 1, 0, 0, 0, 0],
            ["Thing-2", "d__Bacteria",  0, 0, 1, 0, 0, 1, 1],
            ["Thing3", "d__Archaea",  0, 0, 0, 0, 0, 0, 0],
            ["Thing3", "d__Bacteria",  0, 0, 0, 0, 0, 0, 0],
        ], schema=GET_GENOME_TAXA_RECOVERY_OUTPUT_COLUMNS)

        observed_output = get_genome_taxa_recovery(genome_info, coassemble_info)
        self.assertDataFrameEqual(expected_output, observed_output)

    def test_get_genome_taxa_recovery_no_targets(self):
        genome_info = pl.DataFrame([
            # original genomes
            ["binchicken_co1_1", 63.79, 9.77, "binchicken_co1_1", 'Unclassified'],
            ["binchicken_co1_2", 58.74, 0.55, "binchicken_co1_2", 'd__Bacteria;p__Desulfobacterota;c__Desulfobaccia;o__Desulfobaccales;f__0-14-0-80-60-11;g__DTHB01;s__'],
            ["binchicken_co1_3", 67.64, 7.88, "binchicken_co1_3", 'd__Bacteria;p__Acidobacteriota;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__WJMT01;s__'],
            ["binchicken_co1_4", 69.66, 3.46, "binchicken_co1_4", 'd__Bacteria;p__Planctomycetota;c__GCA-2746535;o__;f__;g__;s__'],
            ["binchicken_co1_5", 96.22, 4.06, "binchicken_co1_5", 'd__Archaea;p__Hadarchaeota;c__Hadarchaeia;o__;f__;g__;s__'],
            # mock genomes
            ["binchicken_co1_6", 90.11, 5.22, "binchicken_co1_6", 'Unclassified'],
            ["binchicken_co1_7", 90.11, 5.22, "binchicken_co1_7", 'd__Bacteria;p__Thing1;c__;o__;f__;g__;s__'],
            ["binchicken_co1_8", 90.11, 5.22, "binchicken_co1_8", 'd__Bacteria;p__Thing1;c__;o__;f__;g__;s__'],
            ["binchicken_co1_9", 90.11, 5.22, "binchicken_co1_9", 'd__Bacteria;p__Thing2;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__WJMT01;s__'],
            ["binchicken_co1_10", 90.11, 5.22, "binchicken_co1_10", 'd__Bacteria;p__Thing2;c__GCA-2746535;o__;f__;g__;s__'],
            ["binchicken_co1_11", 90.11, 5.22, "binchicken_co1_11", 'd__Archaea;p__Thing2;c__Hadarchaeia2;o__;f__;g__;s__'],
            ["binchicken_co1_12", 90.11, 5.22, "binchicken_co1_12", 'd__Bacteria;p__Thing2;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__WJMT01;s__mock'],
            # All fail quality
            ["binchicken_co1_13", 4.52, 0.0, "binchicken_co1_13", 'Unclassified'],
            ["binchicken_co1_14", 18.78, 0.82, "binchicken_co1_14", 'd__Bacteria;p__Patescibacteria;c__Gracilibacteria;o__UM-FILTER-43-11;f__UM-FILTER-43-11;g__;s__'],
            ["binchicken_co1_15", 41.63, 1.26, "binchicken_co1_15", 'd__Bacteria;p__OLB16;c__OLB16;o__SURF-12;f__;g__;s__'],
            ["binchicken_co1_16", 22.0, 0.72, "binchicken_co1_16", 'd__Bacteria;p__Chlamydiota;c__Chlamydiia;o__Chlamydiales;f__SM23-39;g__;s__'],
            ["binchicken_co1_17", 32.93, 0.23, "binchicken_co1_17", 'd__Archaea;p__Micrarchaeota;c__Micrarchaeia;o__Anstonellales;f__;g__;s__'],
        ], schema=GET_GENOME_NOVELTY_INPUT_COLUMNS)

        coassemble_info = pl.DataFrame([
            ["Thing1_co1", None],
            ["Thing2_co1", None],
            ["Thing3_co1", None],
            ["OtherThing_co1", None],
        ], schema=GET_GENOME_TAXA_RECOVERY_INPUT_COLUMNS, orient="row")

        expected_output = pl.DataFrame([
        ], schema=GET_GENOME_TAXA_RECOVERY_OUTPUT_COLUMNS)

        observed_output = get_genome_taxa_recovery(genome_info, coassemble_info)
        self.assertDataFrameEqual(expected_output, observed_output)

    def test_update_metadata(self):
        genome_metadata = pl.DataFrame([
            [
                'SRP110325_co2_1', 'bin.1891', 'SRP110325_co2', 'results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover',
                'results/aviary/phylum_2sample/20230718/SRP110325/coassemble/coassemble/coassembly_1_semibin2/recover/bins/final_bins/bin.1891.fna', None,
                729355, 140, 5592, 15740, 0.35, 0.205, 11, 4.52, 0.0,
                'Neural Network (Specific Model)', 'CheckM2v1.0.2',
                'Unclassified', None, None, 'gtdbtk-2.3.0', None, 'No bacterial or archaeal marker'
            ],
            [
                'binchicken_co1_1', 'concoct_bins.tsv.183_sub', 'binchicken_co1', 'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover',
                'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/concoct_bins.tsv.183_sub.fna', 'concoct',
                1685615, 661, 2610, 7076, 0.59, 0.83, 11, 58.74, 0.55,
                'Gradient Boost (General Model)', 'CheckM2v1.0.2',
                'd__Bacteria;p__Desulfobacterota;c__Desulfobaccia;o__Desulfobaccales;f__0-14-0-80-60-11;g__DTHB01;s__', 0.86156, 'taxonomic classification defined by topology and ANI', 'gtdbtk-2.3.0', 'classification based on placement in class-level tree', None
            ],
            [
                'binchicken_co1_8', 'vamb_bins.tsv.182_sub', 'binchicken_co1', 'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover',
                'results/aviary/target/20230619/p__Abyssobacteria/coassemble/coassemble/coassembly_1/recover/bins/final_bins/vamb_bins.tsv.182_sub.fna', 'vamb',
                606959, 237, 2583, 19336, 0.46, 0.926, 11, 32.93, 0.23,
                'Neural Network (Specific Model)', 'CheckM2v1.0.2',
                'd__Archaea;p__Micrarchaeota;c__Micrarchaeia;o__Anstonellales;f__;g__;s__', 0.55037, 'taxonomic classification fully defined by topology', 'gtdbtk-2.3.0', None, None
            ],
        ], schema=GET_GENOME_INFO_OUTPUT_COLUMNS)

        cluster = pl.DataFrame([
            ["/mnt/hpccs01/work/microbiome/msingle/sam/projects/23-SRA-coassembly-r214/results/genomes/20230913/SRP110325_co2/genomes/SRP110325_co2_1.fna", "/mnt/hpccs01/work/microbiome/msingle/sam/projects/23-SRA-coassembly-r214/results/genomes/20230913/SRP110325_co2/genomes/SRP110325_co2_1.fna"],
            ["/mnt/hpccs01/work/microbiome/msingle/sam/projects/23-SRA-coassembly-r214/results/genomes/20230913/binchicken_co1/genomes/binchicken_co1_1.fna", "/mnt/hpccs01/work/microbiome/msingle/sam/projects/23-SRA-coassembly-r214/results/genomes/20230913/binchicken_co1/genomes/binchicken_co1_1.fna"],
            ["/mnt/hpccs01/work/microbiome/msingle/sam/projects/23-SRA-coassembly-r214/results/genomes/20230913/binchicken_co1/genomes/binchicken_co1_8.fna", "/mnt/hpccs01/work/microbiome/msingle/sam/projects/23-SRA-coassembly-r214/results/genomes/20230913/SRP110325_co2/genomes/SRP110325_co2_1.fna"],
        ], schema=CLUSTER_COLUMNS)

        gtdbtk = pl.DataFrame([
            ["SRP110325_co2_1", "d__Archaea", 0.31 ,"topology", "note", "warnings"],
            ["binchicken_co1_1", "d__Archaea;p__Iainarchaeota", 0.32 ,"topology", "note", "warnings"],
            ["binchicken_co1_8", "d__Archaea;p__Micrarchaeota", 0.33 ,"topology", "note", "warnings"],
        ], schema=GTDBTK_COLUMNS)

        gene_counts = pl.DataFrame([
            ["SRP110325_co2_1", 3000],
            ["binchicken_co1_1", 2000],
            ["binchicken_co1_8", 1000],
        ], schema=GENE_COUNTS_COLUMNS)

        gunc = pl.DataFrame([
            ["SRP110325_co2_1", "kingdom", 0.07, 0.4, True],
            ["binchicken_co1_1", "genus", 0.05, 0.5, True],
            ["binchicken_co1_8", "species", 0.30, 0.6, False],
        ], schema=GUNC_COLUMNS)

        rrna = pl.DataFrame([
            ["SRP110325_co2_1", "NODE_1250_length_62074_cov_6.650381", "barrnap:0.9", "rRNA", "100", "200", "0", "-", ".", "Name=23S_rRNA;product=23S ribosomal RNA"],
            ["SRP110325_co2_1", "NODE_1250_length_62074_cov_6.650381", "barrnap:0.9", "rRNA", "100", "200", "0", "-", ".", "Name=23S_rRNA;product=23S ribosomal RNA"],
            ["SRP110325_co2_1", "NODE_1250_length_62074_cov_6.650381", "barrnap:0.9", "rRNA", "100", "200", "0", "-", ".", "Name=5_8S_rRNA;product=5.8S ribosomal RNA (partial);note=aligned only 62 percent of the 5.8S ribosomal RNA"],
            ["SRP110325_co2_1", "NODE_1250_length_62074_cov_6.650381", "barrnap:0.9", "rRNA", "100", "200", "0", "-", ".", "Name=16S_rRNA;product=16S ribosomal RNA"],
            ["binchicken_co1_1", "NODE_1250_length_62074_cov_6.650381", "barrnap:0.9", "rRNA", "100", "200", "0", "-", ".", "Name=23S_rRNA;product=23S ribosomal RNA"],
            ["binchicken_co1_1", "NODE_1250_length_62074_cov_6.650381", "barrnap:0.9", "rRNA", "100", "200", "0", "-", ".", "Name=16S_rRNA;product=16S ribosomal RNA"],
            ["binchicken_co1_8", "NODE_1250_length_62074_cov_6.650381", "barrnap:0.9", "rRNA", "100", "200", "0", "-", ".", "Name=5_8S_rRNA;product=5.8S ribosomal RNA (partial);note=aligned only 62 percent of the 5.8S ribosomal RNA"],
            ["binchicken_co1_8", "NODE_1250_length_62074_cov_6.650381", "barrnap:0.9", "rRNA", "100", "200", "0", "-", ".", "Name=16S_rRNA;product=16S ribosomal RNA"],
        ], schema=RRNA_COLUMNS)

        trna = pl.DataFrame([
            ["SRP110325_co2_1", "NODE_2_length_697855_cov_12.452031", "1", "begin", "end", "Phe", "GAA", "0", "0", "score", None],
            ["SRP110325_co2_1", "NODE_2_length_697855_cov_12.452031", "2", "begin", "end", "Leu", "TAG", "0", "0", "score", None],
            ["SRP110325_co2_1", "NODE_2_length_697855_cov_12.452031", "3", "begin", "end", "Asn", "GTT", "0", "0", "score", None],
            ["binchicken_co1_1", "NODE_2_length_697855_cov_12.452031", "1", "begin", "end", "Phe", "GAA", "0", "0", "score", None],
            ["binchicken_co1_1", "NODE_2_length_697855_cov_12.452031", "2", "begin", "end", "Leu", "TAG", "0", "0", "score", None],
            ["binchicken_co1_1", "NODE_2_length_697855_cov_12.452031", "2", "begin", "end", "Leu", "TAG", "0", "0", "score", None],
            ["binchicken_co1_1", "NODE_2_length_697855_cov_12.452031", "2", "begin", "end", "Undet", "TAG", "0", "0", "score", None],
            ["binchicken_co1_8", "NODE_2_length_697855_cov_12.452031", "1", "begin", "end", "Phe", "GAA", "0", "0", "score", None],
        ], schema=TRNA_COLUMNS)


        expected_output = pl.DataFrame([
            [
                'SRP110325_co2_1', 'bin.1891', 'SRP110325_co2', "SRP110325_co2/genomes/SRP110325_co2_1.fna", None,
                729355, 140, 5592, 15740, 0.35, 0.205, 11, 4.52, 0.0,
                'Neural Network (Specific Model)', 'CheckM2v1.0.2',
                "SRP110325_co2_1",
                "d__Archaea", 0.31 ,"topology", "v2.3.0", "note", "warnings",
                3000,
                "kingdom", 0.07, 0.4, True,
                0, 1, 2, 3,
            ],
            [
                'binchicken_co1_1', 'concoct_bins.tsv.183_sub', 'binchicken_co1', "binchicken_co1/genomes/binchicken_co1_1.fna", 'concoct',
                1685615, 661, 2610, 7076, 0.59, 0.83, 11, 58.74, 0.55,
                'Gradient Boost (General Model)', 'CheckM2v1.0.2',
                "binchicken_co1_1",
                "d__Archaea;p__Iainarchaeota", 0.32 ,"topology", "v2.3.0", "note", "warnings",
                2000,
                "genus", 0.05, 0.5, True,
                0, 1, 1, 2,
            ],
            [
                'binchicken_co1_8', 'vamb_bins.tsv.182_sub', 'binchicken_co1', "binchicken_co1/genomes/binchicken_co1_8.fna", 'vamb',
                606959, 237, 2583, 19336, 0.46, 0.926, 11, 32.93, 0.23,
                'Neural Network (Specific Model)', 'CheckM2v1.0.2',
                "SRP110325_co2_1",
                "d__Archaea;p__Micrarchaeota", 0.33 ,"topology", "v2.3.0", "note", "warnings",
                1000,
                "species", 0.30, 0.6, False,
                0, 1, 0, 1,
            ],
        ], schema=UPDATE_METADATA_OUTPUT_COLUMNS)

        observed_output = update_metadata(genome_metadata, cluster, gtdbtk, gene_counts, gunc, rrna, trna)
        self.assertDataFrameEqual(expected_output, observed_output)


if __name__ == '__main__':
    unittest.main()
