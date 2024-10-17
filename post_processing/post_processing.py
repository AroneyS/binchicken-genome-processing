#!/usr/bin/env python3

"""
Author: Samuel Aroney
Post_processing module
"""

import os
import polars as pl
import yaml
import re

TAXA_TARGETING = "taxa-targeting"
PROJECT_SPECIFIC = "project-specific"

def find_coassemblies(folder):
    return [f.path for f in os.scandir(os.path.join(folder, "coassemble", "coassemble")) if f.is_dir() and not f.path.endswith("logs")]

def get_coassemblies(df):
    output = (
        df
        .with_columns(coassemble_dir = pl.col("binchicken_dir").map_elements(find_coassemblies))
        .explode("coassemble_dir")
        .with_columns(
            assemble_dir = pl.concat_str("coassemble_dir", pl.lit("/assemble")).str.replace("_semibin2", ""),
            recover_dir = pl.concat_str("coassemble_dir", pl.lit("/recover")),
            )
        .drop("coassemble_dir")
    )

    return output

def name_coassemblies(new, old):
    old = (
        old
        .select(
            basename = pl.col("name").str.extract(r"(.*[^\d])"),
            number = pl.col("name").str.extract(r"(\d+$)").cast(int),
            )
        .group_by("basename")
        .agg(pl.max("number"))
    )

    output = (
        new
        .with_columns(
            new_name =
                pl.when(pl.col("style") == TAXA_TARGETING)
                .then(pl.lit("binchicken_co"))
                .otherwise(pl.concat_str("scope", pl.lit("_co"))),
            recover_hash = pl.col("recover_dir").hash(),
            )
        .join(old, left_on="new_name", right_on="basename", how="left")
        .with_columns(
            pl.col("number").fill_null(0),
            ones = pl.lit(1),
            )
        .sort("recover_hash")
        .with_columns(
            row_count = pl.col("ones").cumsum().over("new_name").flatten()
            )
        .select(
            (pl.col("new_name") + (pl.col("number") + pl.col("row_count")).cast(str)).alias("name"),
            "binchicken_dir", "style", "scope", "target", "appraise_level", "aviary_version", "assemble_dir", "recover_dir",
            )
    )

    return output

def read_binchicken_config(path):
    with open(path) as f:
        config = yaml.safe_load(f)

    try:
        binchicken_version = config["Ibis_version"]
    except KeyError:
        binchicken_version = "unknown"

    unmapping_identity = config["unmapping_max_identity"]
    unmapping_alignment = config["unmapping_max_alignment"]

    return {
        "binchicken_version": binchicken_version,
        "unmapping_identity": unmapping_identity,
        "unmapping_alignment": unmapping_alignment,
        }

def get_samples_lists(struct):
    coassembly = re.search(r"coassembly_\d+", struct["recover_dir"])[0]
    elusive_clusters = (
        pl.read_csv(struct["binchicken_dir"], separator="\t")
        .filter(pl.col("coassembly") == coassembly)
        )

    samples = elusive_clusters.get_column("samples").to_list()[0]
    recover_samples = elusive_clusters.get_column("recover_samples").to_list()[0]

    return {
        "original_name": coassembly,
        "samples": samples,
        "recover_samples": recover_samples,
    }

def read_assemble_config(path):
    with open(path) as f:
        config = yaml.safe_load(f)

    assembler = "megahit" if config["use_megahit"] else "metaspades"

    return {
        "assembler": assembler,
        }


def read_assembly_stats(path):
    with open(path) as f:
        for line in f:
            if line.startswith("A\tC\tG\tT\tN\tIUPAC\tOther\tGC"):
                gc = next(f).strip().split("\t")[7]

            if line.startswith("Main genome contig N/L50"):
                contig_N50 = line.strip().split("\t")[1].split("/")[0]

            if line.startswith("Length"):
                next(f)
                all_line = next(f).strip().split("\t")
                contigs_total = all_line[2].strip().replace(",", "")
                assembly_size = all_line[4].strip().replace(",", "")

    return {
        "assembly_size": assembly_size,
        "contigs_total": contigs_total,
        "contig_N50": contig_N50,
        "gc": gc,
        }

def get_coassembly_info(df):
    output = (
        df
        .with_columns(
            binchicken_config = pl.concat_str("binchicken_dir", pl.lit("/config.yaml")).map_elements(read_binchicken_config),
            elusive_clusters = pl.struct(
                pl.concat_str("binchicken_dir", pl.lit("/coassemble/target/elusive_clusters.tsv")),
                pl.col("recover_dir")
                ).map_elements(get_samples_lists),
            assemble_config = pl.concat_str("assemble_dir", pl.lit("/config.yaml")).map_elements(read_assemble_config),
            assembly_stats = pl.concat_str("assemble_dir", pl.lit("/www/assembly_stats.txt")).map_elements(read_assembly_stats),
        )
        .unnest("binchicken_config")
        .unnest("elusive_clusters")
        .unnest("assemble_config")
        .unnest("assembly_stats")
        .select(
            "name", "binchicken_dir", "style", "scope", "target", "appraise_level", "aviary_version", "assemble_dir", "recover_dir",
            "binchicken_version", "original_name", "samples", "recover_samples", "unmapping_identity", "unmapping_alignment", "assembler",
            pl.col("assembly_size").cast(int), pl.col("contigs_total").cast(int), pl.col("contig_N50").cast(int), pl.col("gc").cast(float),
        )
    )

    return output

def find_genomes(folder):
    try:
        output = [os.path.join(folder, "bins", "final_bins", f) for f in os.listdir(os.path.join(folder, "bins", "final_bins")) if f.endswith(".fna")]
    except FileNotFoundError:
        output = []
    return output

def get_genomes(df):
    output = (
        df
        .with_columns(genome_path = pl.col("recover_dir").map_elements(find_genomes))
        .explode("genome_path")
        .filter(pl.col("genome_path").is_not_null())
        .rename({"name": "coassembly"})
    )

    return output

def name_genomes(df):
    output = (
        df
        .with_columns(ones = pl.lit(1))
        .sort("coassembly", "genome_path")
        .with_columns(
            row_count = pl.col("ones").cumsum().over("coassembly").flatten()
            )
        .select(
            (pl.col("coassembly") + pl.lit("_") + pl.col("row_count").cast(str)).alias("name"),
            "coassembly", "recover_dir", "genome_path",
            )
    )

    return output

def get_gtdbtk_version(path):
    with open(path) as f:
        config = yaml.safe_load(f)

    gtdbtk_path = config["gtdbtk_folder"]

    return re.search(r"gtdb(tk-|_release)(\d+\.)?(\d+\.)?\d+", gtdbtk_path)[0]

def read_bin_info(path):
    bin_info = (
        pl.read_csv(path + "/bins/bin_info.tsv", separator="\t")
        .select(
            recover_dir = pl.lit(path),
            original_name = pl.col("Bin Id"),
            genome_size = pl.col("Genome_Size").cast(int),
            contigs_total = pl.col("Total_Contigs").cast(int),
            contig_N50 = pl.col("Contig_N50").cast(int),
            contig_max = pl.col("Max_Contig_Length").cast(int),
            gc = pl.col("GC_Content").cast(float),
            coding_density = pl.col("Coding_Density").cast(float),
            translation_table = pl.col("Translation_Table_Used").cast(int),
            completeness = pl.col("Completeness (CheckM2)").cast(float),
            contamination = pl.col("Contamination (CheckM2)").cast(float),
            CheckM_model = pl.col("Completeness_Model_Used"),
            CheckM_version = pl.lit("CheckM2v1.0.2"),
            taxonomy = pl.col("classification"),
            red_value = pl.col("red_value").cast(float),
            GTDBtk_method = pl.col("classification_method"),
            GTDBtk_note = pl.col("note"),
            GTDBtk_warnings = pl.col("warnings"),
            GTDBtk_version = pl.lit(path + "/config.yaml").map_elements(get_gtdbtk_version),
            )
    )

    return bin_info

def get_genome_info(df):
    bin_info = (
        pl.concat([
            read_bin_info(p) for p in df.unique("recover_dir").get_column("recover_dir").to_list()
            ])
    )

    original_name = pl.col("genome_path").str.replace(r".*/", "").str.replace(r".[^.]+$", "")

    output = (
        df
        .with_columns(
            original_name = original_name,
            binner = original_name.str.extract(r"(^[^_.]+_)").str.replace(r"_", ""),
            )
        .join(bin_info, on=["recover_dir", "original_name"], how="left")
        .select(
            "name", "original_name", "coassembly", "recover_dir", "genome_path", "binner",
            "genome_size", "contigs_total", "contig_N50", "contig_max", "gc",
            "coding_density", "translation_table", "completeness", "contamination", "CheckM_model", "CheckM_version",
            "taxonomy", "red_value", "GTDBtk_method", "GTDBtk_version", "GTDBtk_note", "GTDBtk_warnings",
        )
    )

    return output

def get_genome_novelty(genome_info, min_quality=50):
    novelties = ["unclassified", "phylum", "class", "order", "family", "genus", "species", "strain"]
    novelty = {
        n: i for i, n in enumerate(novelties)
    }

    all_novelties = (
        pl.DataFrame({"domain": ["d__Archaea", "d__Bacteria"]})
        .join(
            pl.DataFrame({"novelty": novelties[1:]}),
            how="cross",
            )
    )

    output = (
        genome_info
        .filter(pl.col("completeness") - 5 * pl.col("contamination") >= min_quality)
        .filter(pl.col("name") == pl.col("95ANI_representative"))
        .select("taxonomy")
        .with_columns(pl.col("taxonomy").str.split(";"))
        .with_columns(
            pl.col("taxonomy").list.get(0).str.replace("Unclassified", "").alias("domain"),
            pl.col("taxonomy").list.contains("p__").alias("phylum"),
            pl.col("taxonomy").list.contains("c__").alias("class"),
            pl.col("taxonomy").list.contains("o__").alias("order"),
            pl.col("taxonomy").list.contains("f__").alias("family"),
            pl.col("taxonomy").list.contains("g__").alias("genus"),
            pl.col("taxonomy").list.contains("s__").alias("species"),
            pl.col("taxonomy").list.contains("Unclassified").alias("unclassified"),
            )
        .with_columns(
            novelty = 
                pl.when(pl.col("unclassified")).then(pl.lit("unclassified"))
                .when(pl.col("phylum")).then(pl.lit("phylum"))
                .when(pl.col("class")).then(pl.lit("class"))
                .when(pl.col("order")).then(pl.lit("order"))
                .when(pl.col("family")).then(pl.lit("family"))
                .when(pl.col("genus")).then(pl.lit("genus"))
                .when(pl.col("species")).then(pl.lit("species"))
                .otherwise(pl.lit("strain")),
            umbrella = 
                pl.when(pl.col("unclassified")).then(pl.lit(""))
                .when(pl.col("phylum")).then(pl.col("taxonomy").list.get(0))
                .when(pl.col("class")).then(pl.col("taxonomy").list.get(1))
                .when(pl.col("order")).then(pl.col("taxonomy").list.get(2))
                .when(pl.col("family")).then(pl.col("taxonomy").list.get(3))
                .when(pl.col("genus")).then(pl.col("taxonomy").list.get(4))
                .when(pl.col("species")).then(pl.col("taxonomy").list.get(5))
                .otherwise(pl.lit("")),
        )
        .group_by("domain", "novelty", "umbrella", maintain_order=True)
        .count()
        .group_by("domain", "novelty", maintain_order=True)
        .agg(
            distinct = pl.count(),
            total = pl.sum("count"),
            )
        .join(all_novelties, on=["domain", "novelty"], how="outer")
        .fill_null(0)
        .sort("domain", pl.col("novelty").map_dict(novelty))
    )

    return output

def get_genome_taxa_recovery(genome_info, coassemble_info, min_quality=50):
    targets = (
        coassemble_info
        .filter(pl.col("target").is_not_null())
        .get_column("target")
        .unique()
        .to_list()
    )

    all_results = (
        pl.DataFrame({"domain": ["d__Archaea", "d__Bacteria"]})
        .join(
            pl.DataFrame({"target": targets}).with_columns(pl.col("target").cast(str)),
            how="cross",
            )
    )

    output = (
        genome_info
        .filter(pl.col("completeness") - 5 * pl.col("contamination") >= min_quality)
        .filter(pl.col("name") == pl.col("95ANI_representative"))
        .with_columns(target = pl.col("taxonomy").str.extract(r"p__([^;]+)"))
        .filter(pl.col("target").is_in(targets))
        .select("target", "taxonomy")
        .with_columns(pl.col("taxonomy").str.split(";"))
        .with_columns(
            pl.col("taxonomy").list.get(0).str.replace("Unclassified", "").alias("domain"),
            pl.col("taxonomy").list.contains("p__").alias("phylum"),
            pl.col("taxonomy").list.contains("c__").alias("class"),
            pl.col("taxonomy").list.contains("o__").alias("order"),
            pl.col("taxonomy").list.contains("f__").alias("family"),
            pl.col("taxonomy").list.contains("g__").alias("genus"),
            pl.col("taxonomy").list.contains("s__").alias("species"),
            pl.col("taxonomy").list.contains("Unclassified").alias("unclassified"),
            )
        .with_columns(
            novelty = 
                pl.when(pl.col("unclassified")).then(pl.lit("unclassified"))
                .when(pl.col("phylum")).then(pl.lit("phylum"))
                .when(pl.col("class")).then(pl.lit("class"))
                .when(pl.col("order")).then(pl.lit("order"))
                .when(pl.col("family")).then(pl.lit("family"))
                .when(pl.col("genus")).then(pl.lit("genus"))
                .when(pl.col("species")).then(pl.lit("species"))
                .otherwise(pl.lit("strain")),
            )
        .group_by("target", "domain")
        .agg(
            (pl.col("novelty") == "phylum").sum().alias("phyla"),
            (pl.col("novelty") == "class").sum().alias("class"),
            (pl.col("novelty") == "order").sum().alias("order"),
            (pl.col("novelty") == "family").sum().alias("family"),
            (pl.col("novelty") == "genus").sum().alias("genus"),
            (pl.col("novelty") == "species").sum().alias("species"),
            (pl.col("novelty") == "strain").sum().alias("strain"),
            )
        .join(all_results, on=["target", "domain"], how="outer")
        .fill_null(0)
        .sort("target", "domain")
    )

    target_domains = (
        output
        .filter(
            pl.sum_horizontal("phyla", "class", "order", "family", "genus", "species", "strain") > 0
            )
        .group_by("target")
        .agg(selected_domain = pl.col("domain"))
    )

    return (
        output
        .join(target_domains, on="target", how="outer")
        .filter(
            pl.col("selected_domain").is_null() | pl.col("domain").is_in(pl.col("selected_domain"))
            )
        .drop("selected_domain")
    )

def update_metadata(genome, cluster, gtdbtk, gene_counts, gunc, rrna, trna):
    cluster = (
        cluster
        .with_columns(
            pl.col("name").str.replace(r".*/", "").str.replace(r"\.fna$", ""),
            pl.col("95ANI_representative").str.replace(r".*/", "").str.replace(r"\.fna$", ""),
            )
    )

    taxonomy = (
        gtdbtk
        .select(
            name = pl.col("user_genome"),
            taxonomy = pl.col("classification"),
            red_value = pl.col("red_value").cast(float),
            GTDBtk_method = pl.col("classification_method"),
            GTDBtk_version = pl.lit("v2.3.0"),
            GTDBtk_note = pl.col("note"),
            GTDBtk_warnings = pl.col("warnings"),
            )
    )

    gunc = (
        gunc
        .select(
            name = pl.col("genome").str.replace(r"\.fna$", ""),
            GUNC_tax_level = pl.col("taxonomic_level"),
            GUNC_contamination = pl.col("contamination_portion"),
            GUNC_representation = pl.col("reference_representation_score"),
            GUNC_pass = pl.col("pass.GUNC"),
        )
    )

    rrna = (
        rrna
        .with_columns(rrna = pl.col("attribute").str.extract(r"Name=(\d+S_rRNA)"))
        .filter(pl.col("rrna").is_in(["5S_rRNA", "16S_rRNA", "23S_rRNA"]))
        .group_by("name", "rrna")
        .count()
        .vstack(
            pl.DataFrame({
                "name": ["dummy", "dummy", "dummy"],
                "rrna": ["5S_rRNA", "16S_rRNA", "23S_rRNA"],
                "count": [0, 0, 0],
                })
            .with_columns(pl.col("count").cast(pl.UInt32))
            )
        .pivot(index="name", columns="rrna", values="count", aggregate_function=None)
        .filter(pl.col("name") != "dummy")
        .fill_null(0)
        .select(
            "name",
            rRNA_5S = "5S_rRNA",
            rRNA_16S = "16S_rRNA",
            rRNA_23S = "23S_rRNA",
            )
    )

    common_trnas = [
        # Standard 20 amino acid tRNAs
        "Ala", "Arg", "Asn", "Asp", "Cys",
        "Gln", "Glu", "Gly", "His", "Ile",
        "Leu", "Lys", "Met", "Phe", "Pro",
        "Ser", "Thr", "Trp", "Tyr", "Val",
        # Suppressor tRNA (stop-codon recoding)
        # "Sup",
        # Selenocysteine tRNA
        # "SeC",
        # Initiator methionine tRNA
        # "iMet",
        # N-formylmethionine tRNA
        # "fMet",
        # Isoleucine AUA-codon tRNA
        # "Ile2",
        # Unknown tRNA
        # "Undet",
    ]

    trna = (
        trna
        .filter(pl.col("type").is_in(common_trnas))
        .group_by("name", "type")
        .count()
        .group_by("name")
        .count()
        .select("name", tRNAs = "count")
    )

    return (
        genome
        .drop(
            # External paths
            "recover_dir", "genome_path",
            # Old GTDBtk columns
            "taxonomy", "red_value", "GTDBtk_method", "GTDBtk_version", "GTDBtk_note", "GTDBtk_warnings"
            )
        .with_columns(genome_path = pl.concat_str("coassembly", pl.lit("/genomes/"), "name", pl.lit(".fna")))
        .select("name", "original_name", "coassembly", "genome_path", pl.exclude("name", "original_name", "coassembly", "genome_path"))
        .join(cluster, on="name", how="left")
        .join(taxonomy, on="name", how="left")
        .join(gene_counts, on="name", how="left")
        .join(gunc, on="name", how="left")
        .join(rrna, on="name", how="left")
        .join(trna, on="name", how="left")
    )
