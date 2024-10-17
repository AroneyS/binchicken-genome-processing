"""
Post-processing of Bin chicken runs

snakemake \
  --snakefile scripts/post_processing.smk \
  --config \
  input_dir=[input directory config] \
  stage_output=[repeated output folder] \
  coassembly_input=[previous coassembly metadata] \
  genome_input=[previous genome metadata] \
  --directory [continuous output folder] \
  --use-conda --cores 64

Input dir:
    binchicken_dir: Bin chicken output directory
    style: project-specific or taxa-targeting
    scope: [project] or global
    target: null or [target]
    appraise_level: Appraise level (global % cutoff or taxa-level)
    aviary_version: vN.N.N

Coassemble metadata
    name: new_name in form of either [project]_co1 or [target]_co1  - generated
    style: project-specific or taxa-targeting                       - input
    scope: [project] or global                                      - input
    target: null or [target]                                        - input
    binchicken_version: binchicken version                          - binchicken config
    appraise_level: Appraisal level (global % cutoff or taxa-level) - input
    samples: coassembly samples                                     - elusive_clusters.tsv
    recover_samples: samples for diff. abund. binning               - elusive_clusters.tsv
    unmapping_identity: % identity cutoff for unmapping             - binchicken config
    unmapping_alignment: % alignment cutoff for unmapping           - binchicken config
    aviary_version: Aviary version                                  - input
    assemble_dir: Aviary assemble dir                               - search for assemble within coassemble
    assembler: metaspades/megahit                                   - "use_megahit" in assemble config
    assembly_size                                                   - "Total Contig Length":"All" in assembly_stats.txt
    contigs_total                                                   - "Number of Contigs":"All" in assembly_stats.txt
    contig_N50                                                      - "Main genome contig N/L50" in assembly_stats.txt
    gc                                                              - "GC" in assembly_stats.txt
    recover_dir: Aviary recover dir                                 - search for recover within coassemble and match with elusive_clusters.tsv

Genome metadata
    name: bin new_name in form of [co name]_1 - generated
    original_name: Bin Id                     - "Bin Id" in bin_info.tsv
    coassembly: coassembly new_name           - key
    binner: source binner                     - within original bin name
    genome_size                               - "Genome_Size" in bin_info.tsv
    contigs_total                             - "Total_Contigs" in bin_info.tsv
    contig_N50                                - "Contig_N50" in bin_info.tsv
    contig_max                                - "Max_Contig_Length" in bin_info.tsv
    gc                                        - "GC_Content" in bin_info.tsv
    coding_density                            - "Coding_Density" in bin_info.tsv
    translation_table                         - "Translation_Table_Used" in bin_info.tsv
    completeness: CheckM completeness         - "Completeness (CheckM2)" in bin_info.tsv
    contamination: CheckM contamination       - "Contamination (CheckM2)" in bin_info.tsv
    CheckM_model: completeness model          - "Completeness_Model_Used" in bin_info.tsv
    CheckM_version: CheckM2vN.N.N             - set to CheckM2v1.0.2
    taxonomy: GTDBtk taxonomy                 - "classification" in bin_info.tsv
    red_value                                 - "red_value" in bin_info.tsv
    GTDBtk_method:                            - "classification_method" in bin_info.tsv
    GTDBtk_version: GTDBtk version            - "gtdbtk_folder" in recover config
    GTDBtk_note                               - "note" in bin_info.tsv
    GTDBtk_warnings                           - "warnings" in bin_info.tsv
"""

import os
import polars as pl
import sys
# To enable snakemake to load local modules
sys.path.insert(0, Path(workflow.basedir).parent.as_posix())
from post_processing import TAXA_TARGETING, PROJECT_SPECIFIC, get_coassemblies, name_coassemblies, get_coassembly_info, get_genomes, name_genomes, get_genome_info, get_genome_novelty, get_genome_taxa_recovery, update_metadata
METAPACKAGE = "/home/aroneys/m/db/singlem/S3.2.0.GTDB_r214.metapackage_20230428.smpkg"
GTDBTK_MASH = "/work/microbiome/db/gtdb/gtdb_release220/mashdb"
GTDB_ARCHAEAL_TAXONOMY = "/work/microbiome/db/gtdb/gtdb_release220/ar53_taxonomy_r220.tsv"
GTDB_BACTERIAL_TAXONOMY = "/work/microbiome/db/gtdb/gtdb_release220/bac120_taxonomy_r220.tsv"
NON_GTDB_GENOMES = "/work/microbiome/ibis/SRA/inputs/spire_genomes2.txt"
NON_GTDB_CHECKM2 = "/work/microbiome/ibis/SRA/data/spire_checkm2_qualities.tsv"
GUNC_DB = "/home/aroneys/s/aroneys/db/GUNC/gunc_db_progenomes2.1.dmnd"
EGGNOG_DATA = "/work/microbiome/db/eggnog-mapper/2.1.3"

wildcard_constraints:
    genome="[a-zA-Z0-9_]+"

######################
### Setup metadata ###
######################
try:
    cumulative_names = pl.read_csv("compiled/coassembly_input_metadata.tsv", separator="\t")
except KeyError:
    cumulative_names = pl.DataFrame(schema = {"name": pl.Utf8})

INPUT_DIR_COLUMNS={
    "binchicken_dir": str,
    "style": str,
    "scope": str,
    "target": str,
    "appraise_level": str,
    "aviary_version": str,
    }

os.makedirs(config["stage_output"], exist_ok=True)

try:
    input_dir = pl.read_csv(config["input_dir"], separator="\t", dtypes=INPUT_DIR_COLUMNS)

    new_coassembly_input_metadata = (
        input_dir
        .pipe(get_coassemblies)
        .pipe(name_coassemblies, old=cumulative_names)
        .pipe(get_coassembly_info)
    )

    new_coassemblies = (
        new_coassembly_input_metadata
        .select("name")
    )

    new_genome_input_metadata = (
        new_coassembly_input_metadata
        .select("name", "recover_dir")
        .pipe(get_genomes)
        .pipe(name_genomes)
        .pipe(get_genome_info)
    )

    (
        new_coassemblies
        .write_csv(os.path.join(config["stage_output"], "new_coassembly_list.tsv"), separator="\t")
    )

    (
        new_genome_input_metadata
        .with_columns(pl.col("name").alias("95ANI_representative"))
        .pipe(get_genome_novelty)
        .write_csv(os.path.join(config["stage_output"], "new_genome_novelty.tsv"), separator="\t")
    )

    (
        new_genome_input_metadata
        .with_columns(pl.col("name").alias("95ANI_representative"))
        .pipe(get_genome_taxa_recovery, coassemble_info=new_coassembly_input_metadata)
        .write_csv(os.path.join(config["stage_output"], "new_genome_recovery.tsv"), separator="\t")
    )

    coassembly_input_metadata = pl.concat([
            pl.read_csv("compiled/coassembly_input_metadata.tsv", separator="\t").with_columns(pl.col("appraise_level").cast(str)),
            new_coassembly_input_metadata,
        ])

    genome_input_metadata = pl.concat([
            pl.read_csv("compiled/genome_input_metadata.tsv", separator="\t"),
            new_genome_input_metadata,
        ])

    coassembly_input_metadata.sort("name").write_csv("compiled/coassembly_input_metadata.tsv", separator="\t")
    genome_input_metadata.sort("name").write_csv("compiled/genome_input_metadata.tsv", separator="\t")
except KeyError:
    coassembly_input_metadata = pl.read_csv("compiled/coassembly_input_metadata.tsv", separator="\t").with_columns(pl.col("appraise_level").cast(str))
    new_coassembly_input_metadata = None
    try:
        new_coassemblies = pl.read_csv(config["stage_output"] + "/new_coassembly_list.tsv", separator="\t")
    except FileNotFoundError:
        new_coassemblies = pl.DataFrame(schema = {"name": pl.Utf8})

    genome_input_metadata = pl.read_csv("compiled/genome_input_metadata.tsv", separator="\t")
    new_genome_input_metadata = genome_input_metadata.filter(pl.col("coassembly").is_in(new_coassemblies["name"]))

manually_disqualified = pl.read_csv("compiled/manually_disqualified.tsv", separator="\t")

#################
### Functions ###
#################
def get_genomes(coassembly, min_quality=50, genome_metadata=genome_input_metadata):
    return (
        genome_metadata
        .filter(pl.col("coassembly") == coassembly)
        .filter(pl.col("completeness") - 5 * pl.col("contamination") >= min_quality)
        .join(manually_disqualified, on="name", how="anti")
        .with_columns(name = pl.concat_str(pl.lit(coassembly + "/genomes/"), "name", pl.lit(".fna")))
        .get_column("name")
        .to_list()
    )

def get_genome_names(coassembly, min_quality=50, genome_metadata=genome_input_metadata):
    return (
        genome_metadata
        .filter(pl.col("coassembly") == coassembly)
        .filter(pl.col("completeness") - 5 * pl.col("contamination") >= min_quality)
        .join(manually_disqualified, on="name", how="anti")
        .get_column("name")
        .to_list()
    )

def get_all_genomes(min_quality=50, genome_metadata=genome_input_metadata):
    return (
        genome_metadata
        .filter(pl.col("completeness") - 5 * pl.col("contamination") >= min_quality)
        .join(manually_disqualified, on="name", how="anti")
        .with_columns(name = pl.concat_str("coassembly", pl.lit("/genomes/"), "name", pl.lit(".fna")))
        .get_column("name")
        .to_list()
    )

def get_domain_genomes(domain, min_quality=50, genome_metadata=genome_input_metadata):
    return (
        genome_metadata
        .filter(pl.col("taxonomy").str.contains("d__" + domain))
        .filter(pl.col("completeness") - 5 * pl.col("contamination") >= min_quality)
        .join(manually_disqualified, on="name", how="anti")
        .select("genome_path", "name", "translation_table")
    )

def get_genome_domain(wildcards, genome_metadata=genome_input_metadata):
    archaeal_genomes = (
        genome_metadata
        .filter(pl.col("name") == wildcards.genome)
        .filter(pl.col("taxonomy").str.contains("d__Archaea"))
    )

    return "arc" if archaeal_genomes.height > 0 else "bac"

####################
### Global rules ###
####################
rule all:
    input:
        expand("{coassembly}/done/all_done", coassembly = coassembly_input_metadata["name"].to_list()),
        "compiled/coassembly_metadata.tsv",
        "compiled/genome_metadata.tsv",
        "compiled/genome_list.tsv",
        "compiled/exclude_coassemblies.tsv",
        "compiled/genome_clusters.tsv",
        "compiled/genome_novelty.tsv",
        "compiled/genome_recovery.tsv",
        "compiled/archaeal_gtdbtk/done",
        "compiled/bacterial_gtdbtk/done",
        "compiled/global_clusters.tsv",
        "compiled/evaluate.tsv",
        "compiled/eggnog.tsv",
        "compiled/eggnog_filtered.tsv",
        "compiled/new_genome.done" if new_genome_input_metadata.height > 0 else [],
        "compiled/new_genome_singlem.done"  if new_coassemblies.height > 0 else [],

rule compile_done:
    input:
        "{coassembly}/done/cleanup.done",
        "{coassembly}/done/bins.done",
        "{coassembly}/done/assembly.done",
        "{coassembly}/done/evaluate.done",
        "{coassembly}/done/transcripts.done",
        "{coassembly}/done/proteins.done",
        "{coassembly}/done/singlem.done",
        "{coassembly}/done/gene_counts.csv",
        "{coassembly}/done/gunc.done",
        "{coassembly}/done/trna.tsv",
        "{coassembly}/done/rrna.tsv",
        "{coassembly}/done/eggnog.tsv",
    output:
        done = "{coassembly}/done/all_done",
    localrule: True
    shell:
        "touch {output.done}"

rule compile_evaluate:
    input:
        expand("{coassembly}/done/evaluate.done", coassembly = coassembly_input_metadata["name"].to_list()),
    output:
        "compiled/evaluate.tsv",
    params:
        matched_hits = expand("{coassembly}/evaluate/evaluate/evaluate/matched_hits.tsv", coassembly = coassembly_input_metadata["name"].to_list()),
    localrule: True
    run:
        queries = []
        for f in params.matched_hits:
            if os.path.exists(f):
                q = pl.scan_csv(f, separator="\t").with_columns(coassembly = pl.lit(f).str.extract(r"([^/]+)/evaluate/evaluate/evaluate/matched_hits.tsv"))
                queries.append(q)

        (
            pl.concat(queries)
            .sink_csv(output[0], separator="\t")
        )

rule compile_genomes:
    input:
        expand("{coassembly}/done/bins.done", coassembly = coassembly_input_metadata["name"].to_list()),
    output:
        "compiled/genome_list.tsv",
    params:
        genomes = ["$PWD/" + g for g in get_all_genomes(50)],
    localrule: True
    shell:
        "echo {params.genomes} | tr ' ' '\n' > {output}"

rule compile_coassemblies:
    input:
        "compiled/genome_list.tsv",
    output:
        "compiled/exclude_coassemblies.tsv",
    params:
        samples = coassembly_input_metadata.get_column("samples").unique().sort().to_list(),
    localrule: True
    shell:
        "echo {params.samples} | tr ' ' '\n' > {output}"

rule coassembly_metadata:
    input:
        "compiled/genome_list.tsv",
    output:
        "compiled/coassembly_metadata.tsv",
    threads: 32
    localrule: True
    run:
        (
            coassembly_input_metadata
            .drop("binchicken_dir", "assemble_dir", "recover_dir")
            .write_csv(output[0], separator="\t")
        )

rule genome_checkm:
    input:
        "compiled/genome_list.tsv",
    output:
        "compiled/genome_checkm.tsv",
    threads: 32
    localrule: True
    run:
        (
            genome_input_metadata
            .select(
                "name",
                pl.lit(None).alias("col1"),
                pl.lit(None).alias("col2"),
                pl.lit(None).alias("col3"),
                pl.lit(None).alias("col4"),
                pl.lit(None).alias("col5"),
                pl.lit(None).alias("col6"),
                pl.lit(None).alias("col7"),
                pl.lit(None).alias("col8"),
                pl.lit(None).alias("col9"),
                pl.lit(None).alias("col10"),
                "completeness",
                "contamination",
                pl.lit(0.0).alias("strain_heterogeneity"),
                )
            .write_csv(output[0], separator="\t")
        )

rule genome_clusters:
    input:
        genomes = "compiled/genome_list.tsv",
        checkm = "compiled/genome_checkm.tsv",
    output:
        "compiled/genome_clusters.tsv",
    threads: 64
    params:
        precluster_ani = 90,
        ani = 95,
    resources:
        mem_mb = 500*1000,
        runtime = "24h",
    log:
        "compiled/genome_clusters.log"
    conda:
        "coverm.yml"
    shell:
      "coverm cluster "
      "--genome-fasta-list {input.genomes} "
      "--checkm-tab-table {input.checkm} "
      "--precluster-ani {params.precluster_ani} "
      "--ani {params.ani} "
      "--output-cluster-definition {output} "
      "--threads {threads} "
      "&> {log} "

rule compile_gene_counts:
    input:
        expand("{coassembly}/done/gene_counts.csv", coassembly = coassembly_input_metadata["name"].to_list()),
    output:
        "compiled/gene_counts.csv",
    localrule: True
    shell:
        "cat {input} > {output} "

rule compile_gunc:
    input:
        expand("{coassembly}/done/gunc.done", coassembly = coassembly_input_metadata["name"].to_list()),
    output:
        "compiled/GUNC.tsv",
    localrule: True
    shell:
        "awk 'NR == 1 || FNR > 1' "
        "$(find . -name 'GUNC.progenomes_2.1.maxCSS_level.tsv') "
        "> {output} "

rule compile_rrna:
    input:
        expand("{coassembly}/done/rrna.tsv", coassembly = coassembly_input_metadata["name"].to_list()),
    output:
        "compiled/rRNA.tsv",
    localrule: True
    shell:
        "awk 'NR == 1 || FNR > 1' "
        "{input} "
        "> {output} "

rule compile_trna:
    input:
        expand("{coassembly}/done/trna.tsv", coassembly = coassembly_input_metadata["name"].to_list()),
    output:
        "compiled/tRNA.tsv",
    localrule: True
    shell:
        "awk 'NR == 1 || FNR > 1' "
        "{input} "
        "> {output} "

rule compile_eggnog:
    input:
        expand("{coassembly}/done/eggnog.tsv", coassembly = coassembly_input_metadata["name"].to_list()),
    output:
        "compiled/eggnog.tsv",
    localrule: True
    shell:
        "awk 'NR == 1 || FNR > 1' "
        "{input} "
        "> {output} "

rule filter_eggnog:
    input:
        "compiled/eggnog.tsv",
    output:
        "compiled/eggnog_filtered.tsv",
    threads: 32
    localrule: True
    run:
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
            "genome": str,
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

        (
            pl.read_csv(input[0], separator="\t", dtypes=EGGNOG_COLUMNS)
            .filter(
                pl.col("KEGG_ko")
                .str.contains("|".join(kos_of_interest))
                )
            .write_csv(output[0], separator="\t")
        )

rule genome_metadata:
    input:
        archaea = "compiled/archaeal_gtdbtk/classify",
        bacteria = "compiled/bacterial_gtdbtk/classify",
        clusters = "compiled/genome_clusters.tsv",
        gene_counts = "compiled/gene_counts.csv",
        gunc = "compiled/GUNC.tsv",
        trna = "compiled/tRNA.tsv",
        rrna = "compiled/rRNA.tsv",
    output:
        "compiled/genome_metadata.tsv",
    threads: 32
    params:
        coassemblies = coassembly_input_metadata["name"].to_list(),
    localrule: True
    run:
        cluster = pl.read_csv(input.clusters, separator="\t", has_header=False, new_columns=["95ANI_representative", "name"])

        gtdbtk_columns = {
            "user_genome": str,
            "classification": str,
            "fastani_reference": str,
            "fastani_reference_radius": str,
            "fastani_taxonomy": str,
            "fastani_ani": str,
            "fastani_af": str,
            "closest_placement_reference": str,
            "closest_placement_radius": str,
            "closest_placement_taxonomy": str,
            "closest_placement_ani": str,
            "closest_placement_af": str,
            "pplacer_taxonomy": str,
            "classification_method": str,
            "note": str,
            "other_related_references(genome_id,species_name,radius,ANI,AF)": str,
            "msa_percent": str,
            "translation_table": str,
            "red_value": str,
            "warnings": str,
        }

        trna_columns = {
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

        gtdbtk = (
            pl.concat([
                pl.read_csv(input.archaea + "/gtdbtk.ar53.summary.tsv", separator="\t", null_values="N/A", dtypes=gtdbtk_columns),
                pl.read_csv(input.bacteria + "/gtdbtk.bac120.summary.tsv", separator="\t", null_values="N/A", dtypes=gtdbtk_columns),
                ])
        )

        gene_counts = pl.read_csv(input.gene_counts, has_header=False, new_columns=["name", "gene_count"])
        gunc = pl.read_csv(input.gunc, separator="\t", null_values="nan")
        rrna = pl.read_csv(input.rrna, separator="\t")
        trna = pl.read_csv(input.trna, separator="\t", dtypes=trna_columns)

        (
            genome_input_metadata
            .pipe(update_metadata, cluster=cluster, gtdbtk=gtdbtk, gene_counts=gene_counts, gunc=gunc, rrna=rrna, trna=trna)
            .write_csv(output[0], separator="\t")
        )

rule genome_novelty:
    input:
        "compiled/genome_metadata.tsv",
    output:
        "compiled/genome_novelty.tsv",
    threads: 32
    localrule: True
    run:
        (
            pl.read_csv(input[0], separator="\t")
            .pipe(get_genome_novelty)
            .write_csv(output[0], separator="\t")
        )

rule genome_recovery:
    input:
        "compiled/genome_metadata.tsv",
    output:
        "compiled/genome_recovery.tsv",
    threads: 32
    localrule: True
    run:
        (
            pl.read_csv(input[0], separator="\t")
            .pipe(get_genome_taxa_recovery, coassemble_info=coassembly_input_metadata)
            .write_csv(output[0], separator="\t")
        )

##############
### GTDBtk ###
##############
rule domain_genomes:
    input:
        "compiled/genome_list.tsv",
    output:
        "compiled/{domain}_genomes.tsv",
    threads: 32
    params:
        domain = lambda wildcards: "Archaea" if wildcards.domain == "archaeal" else "Bacteria",
    localrule: True
    run:
        (
            get_domain_genomes(params.domain)
            .write_csv(output[0], separator="\t", has_header=False)
        )

rule gtdbtk_classify:
    input:
        "compiled/{domain}_genomes.tsv",
    output:
        directory("compiled/{domain}_gtdbtk/classify"),
    threads: 64
    params:
        pplacer_threads = 32,
    resources:
        mem_mb = 500*1000,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    log:
        "compiled/{domain}_gtdbtk/classify.log"
    conda:
        "gtdbtk.yaml"
    shell:
        "gtdbtk classify_wf "
        "--mash_db {GTDBTK_MASH} "
        "--batchfile {input} "
        "--pplacer_cpus {params.pplacer_threads} "
        "--out_dir {output} "
        "--cpus {threads} "
        "&> {log} "

rule domain_representatives:
    input:
        "compiled/genome_metadata.tsv",
    output:
        "compiled/{domain}_representatives.tsv",
    threads: 32
    params:
        domain = lambda wildcards: "Archaea" if wildcards.domain == "archaeal" else "Bacteria",
    localrule: True
    run:
        genome_metadata = (
            pl.read_csv(input[0], separator="\t")
            .filter(pl.col("name") == pl.col("95ANI_representative"))
            )

        (
            get_domain_genomes(params.domain, genome_metadata=genome_metadata)
            .select("genome_path", "name")
            .write_csv(output[0], separator="\t", has_header=False)
        )

rule global_clusters:
    input:
        "compiled/archaeal_spire_novel.tsv",
        "compiled/bacterial_spire_novel.tsv",
        "compiled/gtdbr220_genome_reps.tsv",
        "compiled/archaeal_representatives.tsv",
        "compiled/bacterial_representatives.tsv",
    output:
        "compiled/global_clusters.tsv",
    threads: 64
    params:
        precluster_ani = 90,
        ani = 95,
    resources:
        mem_mb = 500*1000,
        runtime = lambda wildcards, attempt: 24*60*attempt,
    log:
        "compiled/global_clusters.log"
    conda:
        "coverm.yml"
    shell:
      "coverm cluster "
      "--genome-fasta-list <(cut -f1 {input}) "
      "--precluster-ani {params.precluster_ani} "
      "--ani {params.ani} "
      "--output-cluster-definition {output} "
      "--threads {threads} "
      "&> {log} "

rule domain_comb_reps:
    input:
        "compiled/global_clusters.tsv",
        "compiled/{domain}_representatives.tsv",
        "compiled/{domain}_spire_novel.tsv",
    output:
        "compiled/{domain}_comb_reps.tsv",
    threads: 1
    localrule: True
    run:
        global_clusters = (
            pl.read_csv(input[0], separator="\t", has_header=False, new_columns=["95ANI_representative", "name"])
            .with_columns(nonbc = pl.col("name").str.contains("_co").not_())
            .with_columns(nonbc_cluster = pl.col("nonbc").any().over("95ANI_representative"))
            .filter(pl.col("nonbc_cluster").not_())
            .filter(pl.col("95ANI_representative") == pl.col("name"))
            .select(path = pl.col("name"))
        )
        binchicken_novel = (
            pl.read_csv(input[1], separator="\t", has_header=False, new_columns=["path", "name"])
            .join(global_clusters, on="path", how="inner")
        )
        spire_novel = pl.read_csv(input[2], separator="\t", has_header=False, new_columns=["path", "name"])

        (
            pl.concat([binchicken_novel, spire_novel])
            .write_csv(output[0], separator="\t", has_header=False)
        )

rule gtdbtk_denovo:
    input:
        genomes = "compiled/{domain}_comb_reps.tsv",
        classifications = "compiled/{domain}_gtdbtk/classify",
    output:
        dir = directory("compiled/{domain}_gtdbtk/de_novo"),
    threads: 64
    params:
        summary = lambda wildcards: "gtdbtk.ar53.summary.tsv" if wildcards.domain == "archaeal" else "gtdbtk.bac120.summary.tsv",
        outgroup = lambda wildcards: "p__Altiarchaeota" if wildcards.domain == "archaeal" else "p__Fusobacteriota",
        domain = lambda wildcards: "--archaea" if wildcards.domain == "archaeal" else "--bacteria",
    resources:
        mem_mb = 500*1000,
        runtime = lambda wildcards, attempt: 48*60*attempt,
    log:
        "compiled/{domain}_gtdbtk/denovo.log"
    conda:
        "gtdbtk.yaml"
    shell:
        "gtdbtk de_novo_wf "
        "--batchfile {input.genomes} "
        "--gtdbtk_classification_file {input.classifications}/{params.summary} "
        "{params.domain} "
        "--outgroup_taxon {params.outgroup} "
        "--out_dir {output.dir} "
        "--cpus {threads} "
        "&> {log} "

rule phylogenetic_diversity:
    input:
        tree = "compiled/{domain}_gtdbtk/de_novo",
        taxa = "compiled/{domain}_representatives.tsv",
    output:
        pd = "compiled/{domain}_gtdbtk/phylogenetic_diversity.tsv",
        taxa = "compiled/{domain}_gtdbtk/pd_taxa.tsv",
        per_taxa = "compiled/{domain}_gtdbtk/pd_per_taxa.tsv",
    threads: 64
    params:
        tree = lambda wildcards: "gtdbtk.ar53.decorated.tree" if wildcards.domain == "archaeal" else "gtdbtk.bac120.decorated.tree",
    localrule: True
    log:
        "compiled/{domain}_gtdbtk/phylogenetic_diversity.log"
    conda:
        "genometreetk.yaml"
    shell:
        "cut -f2 {input.taxa} > {output.taxa} && "
        "genometreetk pd "
        "{input.tree}/{params.tree} "
        "{output.taxa} "
        "--per_taxa_pg_file {output.per_taxa} "
        "> {output.pd} "
        "2> {log} "

rule phylogenetic_diversity_clades:
    input:
        tree = "compiled/{domain}_gtdbtk/de_novo",
        taxa = "compiled/{domain}_gtdbtk/pd_taxa.tsv",
    output:
        pd = "compiled/{domain}_gtdbtk/phylogenetic_diversity_clades.tsv",
        done = "compiled/{domain}_gtdbtk/clades.done",
    threads: 64
    params:
        tree = lambda wildcards: "gtdbtk.ar53.decorated.tree" if wildcards.domain == "archaeal" else "gtdbtk.bac120.decorated.tree",
    localrule: True
    log:
        "compiled/{domain}_gtdbtk/phylogenetic_diversity_clades.log"
    conda:
        "genometreetk.yaml"
    shell:
        "genometreetk pd_clade "
        "{input.tree}/{params.tree} "
        "{input.taxa} "
        "{output.pd} "
        "&> {log} "
        "&& touch {output.done} "

rule phylogenetic_growth_summary:
    input:
        pd = "compiled/{domain}_gtdbtk/phylogenetic_diversity_clades.tsv",
        done = "compiled/{domain}_gtdbtk/clades.done",
    output:
        "compiled/{domain}_gtdbtk/phylogenetic_growth_summary.tsv",
    threads: 32
    localrule: True
    log:
        "compiled/{domain}_gtdbtk/phylogenetic_growth_summary.log"
    run:
        (
            pl.read_csv(input.pd, separator="\t")
            .filter(pl.col("Clade").str.contains("d__|p__"))
            .select(
                clade = 
                    pl.when(pl.col("Clade").str.contains("d__"))
                    .then(pl.col("Clade").str.extract(r"(d__[^;]+)"))
                    .otherwise(pl.col("Clade").str.extract(r"(p__[^;]+)")),
                prior = pl.col("Out Taxa"),
                binchicken = pl.col("In Taxa"),
                pg = pl.col("In Percent PG").round(1),
                sg = (pl.col("In Taxa") / pl.col("Out Taxa") * 100).round(1),
                )
            .sort(pl.col("clade").str.extract(r"^(.)"), "pg", descending=[False, True])
            .write_csv(output[0], separator="\t")
        )

rule phylorank_red:
    input:
        tree = "compiled/{domain}_gtdbtk/de_novo",
    output:
        directory("compiled/{domain}_gtdbtk/phylorank"),
    threads: 64
    params:
        tree = lambda wildcards: "gtdbtk.ar53.decorated.tree" if wildcards.domain == "archaeal" else "gtdbtk.bac120.decorated.tree",
        taxonomy = lambda wildcards: GTDB_ARCHAEAL_TAXONOMY if wildcards.domain == "archaeal" else GTDB_BACTERIAL_TAXONOMY,
    localrule: True
    log:
        "compiled/{domain}_gtdbtk/phylorank.log"
    conda:
        "phylorank.yml"
    shell:
        "phylorank outliers "
        "{input.tree}/{params.tree} "
        "{params.taxonomy} "
        "{output} "
        "&> {log} "

rule tree2tbl:
    input:
        dir = "compiled/{domain}_gtdbtk/phylorank",
    output:
        tree = "compiled/{domain}_gtdbtk/phylorank_tree.tsv",
    threads: 1
    params:
        domain = lambda wildcards: "archaea" if wildcards.domain == "archaeal" else "bacteria",
        extra_genomes = NON_GTDB_GENOMES,
        input_tree = lambda wildcards: "gtdbtk.ar53.decorated.red_decorated.tree" if wildcards.domain == "archaeal" else "gtdbtk.bac120.decorated.red_decorated.tree",
    localrule: True
    log:
        "compiled/{domain}_gtdbtk/tree2tbl.log"
    conda:
        "r-treedataverse"
    script:
        "tree2tbl.R"

rule genome_checkm2:
    input:
        genomes = "compiled/genome_list.tsv",
        checkm = "compiled/genome_checkm.tsv",
    output:
        "compiled/genome_checkm2.tsv",
    threads: 1
    params:
        spire = NON_GTDB_CHECKM2
    localrule: True
    shell:
        "cut -f 1,2,3 {params.spire} > {output} && "
        "cut -f1,12,13 {input.checkm} | tail -n+2 >> {output}"

rule name_clades:
    input:
        tree = "compiled/{domain}_gtdbtk/phylorank_tree.tsv",
        checkm = "compiled/genome_checkm2.tsv",
        # clusters = "compiled/global_clusters.tsv",
    output:
        directory("compiled/{domain}_gtdbtk/name_clades"),
    threads: 64
    params:
        script = "name_clades.py",
        domain = lambda wildcards: "d__Archaea" if wildcards.domain == "archaeal" else "d__Bacteria",
        taxonomy = lambda wildcards: GTDB_ARCHAEAL_TAXONOMY if wildcards.domain == "archaeal" else GTDB_BACTERIAL_TAXONOMY,
    localrule: True
    log:
        "compiled/{domain}_gtdbtk/name_clades.log"
    shell:
        "python {params.script} "
        "--tree-df {input.tree} "
        "--metadata {input.checkm} "
        "--gtdb {params.taxonomy} "
        "--domain {params.domain} "
        # "--clusters {input.clusters} "
        "--output {output} "
        "&> {log} "

rule clade_novelty:
    input:
        "compiled/{domain}_comb_reps.tsv",
        "compiled/global_clusters.tsv",
        "compiled/{domain}_gtdbtk/name_clades",
    output:
        "compiled/{domain}_gtdbtk/clade_novelty.tsv",
        "compiled/{domain}_gtdbtk/distinct_novelty.tsv",
    threads: 16
    localrule: True
    log:
        "compiled/{domain}_gtdbtk/clade_novelty.log"
    run:
        levels = ["p", "c", "o", "f", "g", "s"]
        magsets = ["GTDB", "UHGG", "GEM", "OceanDNA", "SPIRE", "SMAG", "Tengchong", "binchicken"]

        representatives = (
            pl.read_csv(input[0], separator="\t", has_header=False, new_columns=["path", "name"])
            .select(pl.col("path").alias("95ANI_representative"))
        )

        global_clusters = (
            pl.read_csv(input[1], separator="\t", has_header=False, new_columns=["95ANI_representative", "name"])
            .with_columns(magset = pl.col("name").str.extract(r"/work/microbiome/db/([^/]+)"))
            .with_columns(
                magset = 
                    pl.when(pl.col("magset") == "uhgg_v2")
                        .then(pl.lit("UHGG"))
                    .when(pl.col("name").str.contains(r"ocean_mags"))
                        .then(pl.lit("OceanDNA"))
                    .when(pl.col("name").str.contains(r"Tengchong_genomes"))
                        .then(pl.lit("Tengchong"))
                    .when(pl.col("name").str.contains(r"_co"))
                        .then(pl.lit("binchicken"))
                    .when(pl.col("magset").is_not_null())
                        .then(pl.col("magset"))
                    .otherwise(pl.lit("GTDB")),
                genome = pl.col("name").str.extract(r".*/([^/]+).fn?a$")
                )
            .with_columns(order = pl.col("magset").map_elements(lambda x: magsets.index(x)))
        )

        species_counts = (
            global_clusters
            .sort("order")
            .group_by("95ANI_representative", maintain_order=True)
            .first()
            .join(representatives, on="95ANI_representative", how="inner")
            .group_by("magset")
            .count()
            .with_columns(level = pl.lit("s"))
            .select("level", "magset", "count")
        )

        node_names = pl.read_csv(input[2] + "/node_names.tsv", separator="\t")

        available_columns = ["level", "binchicken", "SPIRE", "GEM"]
        if species_counts.filter(pl.col("magset") == "Tengchong").height > 0:
            available_columns.append("Tengchong")
        available_columns.extend(["OceanDNA", "SMAG"])
        if species_counts.filter(pl.col("magset") == "UHGG").height > 0:
            available_columns.append("UHGG")

        (
            pl.concat([
                node_names
                    .select("clade")
                    .unique()
                    .filter(pl.col("clade").str.contains(r"s__").not_())
                    .with_columns(
                        level = pl.col("clade").str.extract(r"(^.)__[^_]+"),
                        magset = pl.col("clade").str.extract(r"^.__([^_]+)")
                        )
                    .group_by("level", "magset")
                    .count(),
                species_counts
                ])
            .pivot(index="level", columns="magset", values="count")
            .fill_null(0)
            .with_columns(order = pl.col("level").map_elements(lambda x: levels.index(x)))
            .sort("order")
            .select(available_columns)
            .write_csv(output[0], separator="\t")
        )

        (
            pl.read_csv(input[2] + "/genome_taxonomy.tsv", separator="\t")
            .join(global_clusters.select("genome", "magset"), how="inner", on="genome")
            .filter(pl.col("magset") == "binchicken")
            .with_columns(pl.col("taxonomy").str.split(";"))
            .explode("taxonomy")
            .join(node_names, how="inner", left_on="taxonomy", right_on="clade")
            .with_columns(level = pl.col("taxonomy").str.extract(r"(^.)__"))
            .unique(["taxonomy", "level"])
            .group_by("level")
            .count()
            .with_columns(order = pl.col("level").map_elements(lambda x: levels.index(x)))
            .sort("order")
            .select("level", "count")
            .write_csv(output[1], separator="\t")
        )

rule comb_taxonomy:
    input:
        names = "compiled/{domain}_gtdbtk/name_clades",
    output:
        "compiled/{domain}_gtdbtk/comb_named_taxonomy.tsv",
    threads: 1
    params:
        taxonomy = lambda wildcards: GTDB_ARCHAEAL_TAXONOMY if wildcards.domain == "archaeal" else GTDB_BACTERIAL_TAXONOMY,
    localrule: True
    shell:
        "cat {input.names}/genome_taxonomy.tsv {params.taxonomy} > {output} "

rule decorate_named_tree:
    input:
        tree = "compiled/{domain}_gtdbtk/de_novo",
        comb = "compiled/{domain}_gtdbtk/comb_named_taxonomy.tsv",
    output:
        "compiled/{domain}_gtdbtk/named_decorated.tree",
    threads: 64
    params:
        tree = lambda wildcards: "gtdbtk.ar53.decorated.tree" if wildcards.domain == "archaeal" else "gtdbtk.bac120.decorated.tree",
    localrule: True
    log:
        "compiled/{domain}_gtdbtk/named_decorated.log"
    conda:
        "gtdbtk.yaml"
    shell:
        "gtdbtk decorate "
        "--input_tree {input.tree}/{params.tree} "
        "--custom_taxonomy_file {input.comb} "
        "--output_tree {output} "
        "&> {log} "

rule compile_de_novo:
    input:
        "compiled/{domain}_gtdbtk/phylogenetic_diversity.tsv",
        "compiled/{domain}_gtdbtk/phylogenetic_growth_summary.tsv",
        "compiled/{domain}_gtdbtk/named_decorated.tree",
        "compiled/{domain}_gtdbtk/clade_novelty.tsv",
    output:
        done = "compiled/{domain}_gtdbtk/done",
    localrule: True
    shell:
        "touch {output.done} "

#####################
### Stage outputs ###
#####################
rule compile_new_genomes:
    input:
        expand("{coassembly}/done/bins.done", coassembly = new_coassemblies["name"].to_list()),
    output:
        "compiled/new_genome.done",
    params:
        genomes = ["$PWD/" + g for g in get_all_genomes(50, genome_metadata=new_genome_input_metadata)] if new_genome_input_metadata.height > 0 else [],
        stage_output = config["stage_output"],
        output = "new_genome_list.tsv",
    localrule: True
    shell:
        "echo {params.genomes} | tr ' ' '\n' "
        "> {params.stage_output}/{params.output} "
        "&& touch {output}"

rule compile_new_singlem:
    input:
        expand("{coassembly}/done/evaluate.done", coassembly = new_coassemblies["name"].to_list()),
    output:
        "compiled/new_genome_singlem.done",
    params:
        singlem = expand("{coassembly}/evaluate/evaluate/summarise", coassembly = new_coassemblies["name"].to_list()) if new_coassemblies.height > 0 else [],
        stage_output = config["stage_output"],
        output = "new_genome.otu_table.tsv",
    localrule: True
    shell:
        "awk 'NR == 1 || FNR > 1' "
        "$(find {params.singlem} -name '*.otu_table.tsv') "
        r"| sed 's/coassembly_[[:digit:]]\+-//;s/_transcripts//' "
        "> {params.stage_output}/{params.output} "
        "&& touch {output}"

########################
### Coassembly rules ###
########################
rule remove_reads:
    output:
        done = "{coassembly}/done/cleanup.done",
    params:
        binchicken_dir = lambda wildcards: coassembly_input_metadata.filter(pl.col("name") == wildcards.coassembly).get_column("binchicken_dir").to_list()[0],
    localrule: True
    shell:
        "rm -f {params.binchicken_dir}/coassemble/mapping/* {params.binchicken_dir}/coassemble/sra/* {params.binchicken_dir}/coassemble/sra_qc/*fastq.gz "
        "&& touch {output.done}"

rule copy_assembly:
    output:
        fasta = "{coassembly}/{coassembly}.fasta",
        done = "{coassembly}/done/assembly.done",
    params:
        assemble_dir = lambda wildcards: coassembly_input_metadata.filter(pl.col("name") == wildcards.coassembly).get_column("assemble_dir").to_list()[0],
    localrule: True
    shell:
        "cp {params.assemble_dir}/assembly/final_contigs.fasta {output.fasta} "
        "&& touch {output.done}"

rule copy_genomes:
    output:
        folder = directory("{coassembly}/genomes"),
        done = "{coassembly}/done/bins.done",
    threads: 32
    params:
        old_paths = lambda wildcards: genome_input_metadata.filter(pl.col("coassembly") == wildcards.coassembly).get_column("genome_path").to_list(),
        new_names = lambda wildcards: genome_input_metadata.filter(pl.col("coassembly") == wildcards.coassembly).get_column("name").to_list(),
    localrule: True
    shell:
        "mkdir -p {output.folder} "
        "&& parallel -j {threads} "
        "cp {{1}} {output.folder}/{{2}}.fna "
        "::: {params.old_paths} :::+ {params.new_names} "
        "&& touch {output.done}"

rule run_evaluate:
    input:
        "{coassembly}/done/bins.done",
    output:
        done = "{coassembly}/done/evaluate.done",
        folder = directory("{coassembly}/evaluate"),
    threads: 16
    params:
        binchicken_dir = lambda wildcards: coassembly_input_metadata.filter(pl.col("name") == wildcards.coassembly).get_column("binchicken_dir").to_list()[0],
        coassembly = lambda wildcards: coassembly_input_metadata.filter(pl.col("name") == wildcards.coassembly).get_column("original_name").to_list()[0],
        new_genomes = lambda wildcards: get_genomes(wildcards.coassembly, 50),
        metapackage = METAPACKAGE,
    resources:
        mem_mb = 125*1000,
        runtime = "24h",
    log:
        "{coassembly}/done/evaluate.log",
    conda:
        "binchicken.yml"
    shell:
        "test -z '{params.new_genomes}' && mkdir -p {output.folder} && touch {output.done} && exit 0 || "
        "binchicken evaluate "
        "--output {output.folder} "
        "--coassemble-output {params.binchicken_dir} "
        "--coassembly-run {params.coassembly} "
        "--new-genomes {params.new_genomes} "
        "--singlem-metapackage {params.metapackage} "
        "--cores {threads} "
        "&> {log} "
        "&& touch {output.done} "

rule collate_transcripts:
    input:
        "{coassembly}/evaluate",
    output:
        done = "{coassembly}/done/transcripts.done",
        folder = directory("{coassembly}/transcripts"),
    threads: 32
    localrule: True
    shell:
        "mkdir -p {output.folder} && "
        "find {input}/evaluate/transcripts -name '*_transcripts.fna' | "
        r"parallel -j{threads} ln -sr {{}} {output.folder}/'{{= s=.*/==;s/coassembly_\d+-// =}}' "
        "&& touch {output.done}"
        "|| touch {output.done}"

rule collate_proteins:
    input:
        "{coassembly}/evaluate",
    output:
        done = "{coassembly}/done/proteins.done",
        folder = directory("{coassembly}/proteins"),
    threads: 32
    localrule: True
    shell:
        "mkdir -p {output.folder} && "
        "find {input}/evaluate/transcripts -name '*_transcripts.faa' | "
        r"parallel -j{threads} ln -sr {{}} {output.folder}/'{{= s=.*/==;s/coassembly_\d+-// =}}' "
        "&& touch {output.done}"
        "|| touch {output.done}"

rule collate_singlem:
    input:
        "{coassembly}/evaluate",
    output:
        done = "{coassembly}/done/singlem.done",
        folder = directory("{coassembly}/singlem"),
    threads: 32
    localrule: True
    shell:
        "mkdir -p {output.folder} && "
        "find {input}/evaluate/pipe -name '*.otu_table.tsv' | "
        r"parallel -j{threads} ln -sr {{}} {output.folder}/'{{= s=.*/==;s/coassembly_\d+-// =}}' "
        "&& touch {output.done}"
        "|| touch {output.done}"

rule run_gene_counting:
    input:
        done = "{coassembly}/done/proteins.done",
    output:
        "{coassembly}/done/gene_counts.csv",
    params:
        coassembly = "{coassembly}",
    threads: 1
    localrule: True
    shell:
        "find {params.coassembly}/proteins -name '*.faa' | "
        "parallel -j {threads} "
        "'echo {{= s:.*/::;s:_transcripts.*:: =}},$(grep \"^>\" {{}} | wc -l)' "
        "> {output} "

rule run_gunc:
    input:
        "{coassembly}/done/bins.done",
    output:
        done = "{coassembly}/done/gunc.done",
        folder = directory("{coassembly}/gunc"),
    threads: 16
    params:
        database = GUNC_DB,
        genomes = lambda wildcards: get_genomes(wildcards.coassembly, 50),
    resources:
        mem_mb = 125*1000,
        runtime = "24h",
    log:
        "{coassembly}/done/gunc.log",
    conda:
        "gunc.yml"
    shell:
        "mkdir -p {output.folder} && "
        "echo '{params.genomes}' | tr ' ' '\n' > {output.folder}/genome_list.txt && "
        "gunc run "
        "--input_file {output.folder}/genome_list.txt "
        "--db_file {params.database} "
        "--out_dir {output.folder} "
        "--temp_dir {output.folder} "
        "--threads {threads} "
        "&> {log} "
        "&& touch {output.done} "
        "|| touch {output.done} "

rule run_barrnap:
    input:
        "{coassembly}/done/bins.done",
    output:
        "{coassembly}/rna/{genome}.gff",
    threads: 4
    params:
        genome = "{coassembly}/genomes/{genome}.fna",
        kingdom = get_genome_domain,
    resources:
        mem_mb = 32*1000,
        runtime = "24h",
    log:
        "{coassembly}/rna/{genome}.log",
    conda:
        "barrnap.yml"
    shell:
        "barrnap "
        "--kingdom {params.kingdom} "
        "--threads {threads} "
        "{params.genome} "
        "> {output} "
        "2> {log} "

rule collect_barrnap:
    input:
        gff = lambda wildcards: expand("{coassembly}/rna/{genome}.gff", genome=get_genome_names(wildcards.coassembly, 50), allow_missing=True)
    output:
        "{coassembly}/done/rrna.tsv",
    threads: 32
    localrule: True
    run:
        gff_columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        gff_schema = {
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

        queries = []
        if input.gff:
            for file in input.gff:
                with open(file) as f:
                    if f.read() == "##gff-version 3\n":
                        continue
                queries.append(
                    pl.scan_csv(file, has_header=False, comment_char="#", separator="\t", new_columns=gff_columns, dtypes=gff_schema)
                    .with_columns(path = pl.lit(file))
                )

        if queries:
            (
                pl.concat(pl.collect_all(queries))
                .with_columns(
                    name = pl.col("path").str.extract(r"([^/]+)\.gff$"),
                    )
                .select(["name"] + gff_columns)
                .write_csv(output[0], separator="\t")
            )
        else:
            with open(output[0], "w") as f:
                f.write("\t".join(gff_columns) + "\n")

rule run_trnascan:
    input:
        "{coassembly}/done/bins.done",
    output:
        "{coassembly}/trna/{genome}.out",
    threads: 1
    params:
        genome = "{coassembly}/genomes/{genome}.fna",
        domain = lambda wildcards: "-B" if get_genome_domain(wildcards) == "bac" else "-A",
    resources:
        mem_mb = 8*1000,
        # Sometimes stalls and fill log file with garbage
        # A retry should fix this
        runtime = "24h",
    log:
        "{coassembly}/trna/{genome}.log",
    conda:
        "trnascan.yml"
    shell:
        "tRNAscan-SE "
        "{params.domain} "
        "-o {output} "
        "{params.genome} "
        "&> {log} "

rule collect_trnascan:
    input:
        out = lambda wildcards: expand("{coassembly}/trna/{genome}.out", genome=get_genome_names(wildcards.coassembly, 50), allow_missing=True)
    output:
        "{coassembly}/done/trna.tsv",
    threads: 32
    localrule: True
    run:
        out_columns = ["seqname", "tRNA_no", "begin", "end", "type", "codon", "intron_begin", "intron_end", "score", "note"]
        out_schema = {
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

        queries = []
        if input.out:
            for file in input.out:
                with open(file) as f:
                    if f.read() != "":
                        queries.append(
                            pl.scan_csv(file, has_header=False, skip_rows=3, separator="\t", new_columns=out_columns, dtypes=out_schema)
                            .with_columns(path = pl.lit(file))
                        )

        if queries:
            (
                pl.concat(pl.collect_all(queries))
                .with_columns(
                    name = pl.col("path").str.extract(r"([^/]+)\.out$"),
                    )
                .select(["name"] + out_columns)
                .write_csv(output[0], separator="\t")
            )
        else:
            with open(output[0], "w") as f:
                f.write("\t".join(out_columns) + "\n")

rule run_eggnog:
    input:
        "{coassembly}/done/proteins.done",
    output:
        "{coassembly}/eggnog/{genome}.emapper.annotations",
        "{coassembly}/eggnog/{genome}.emapper.hits",
        "{coassembly}/eggnog/{genome}.emapper.seed_orthologs",
    threads: 5
    params:
        output = "{coassembly}/eggnog/{genome}",
        protein = "{coassembly}/proteins/{genome}_transcripts.faa",
        eggnog_data = os.path.basename(EGGNOG_DATA)
    resources:
        mem_mb = 5*8*1000,
        runtime = "4h",
        extra_mqsub_args = "--scratch-data " + EGGNOG_DATA,
    log:
        "{coassembly}/eggnog/{genome}.log",
    conda:
        "eggnog.yml"
    shell:
        "EGGNOG_DATA_DIR=$TMPDIR/{params.eggnog_data} "
        "emapper.py "
        "-i {params.protein} "
        "-m diamond "
        "--target_orthologs one2one "
        "--query_cover 50.0 "
        "--evalue 0.0000001 "
        "--cpu {threads} "
        "-o {params.output} "
        "&> {log} "

rule collect_eggnog:
    input:
        annotations = lambda wildcards: expand("{coassembly}/eggnog/{genome}.emapper.annotations", genome=get_genome_names(wildcards.coassembly, 50), allow_missing=True)
    output:
        "{coassembly}/done/eggnog.tsv",
    threads: 16
    localrule: True
    run:
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

        queries = []
        if input.annotations:
            for file in input.annotations:
                with open(file) as f:
                    genome = os.path.basename(file).split(".")[0]
                    queries.append(
                        pl.scan_csv(
                            file,
                            separator="\t",
                            has_header=False,
                            new_columns=list(EGGNOG_COLUMNS.keys()),
                            dtypes=EGGNOG_COLUMNS,
                            comment_char="#",
                            )
                        .with_columns(genome = pl.lit(genome))
                    )

        if queries:
            (
                pl.concat(pl.collect_all(queries))
                .write_csv(output[0], separator="\t")
            )
        else:
            with open(output[0], "w") as f:
                f.write("\t".join(list(EGGNOG_COLUMNS) + ["genome"]) + "\n")
