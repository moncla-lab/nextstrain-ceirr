import csv
import json

from ceirr import SEGMENTS
from ceirr import phenotypic_characterization_annotation
from ceirr import split_phenotypes_excel
from ceirr import genoflu_dataflow
from ceirr import genoflu_postprocess

wildcard_constraints:
  segment="[^/]+"

NUMBER_OF_GENOTYPES = 15
GENOTYPES_TO_INCLUDE = ['D1.1', 'D1.2']

rule all:
    input:
        auspice_json = expand("data/ml/h5nx-{segment}-ceirr.json", segment=SEGMENTS)

rule files:
    params:
        input_metadata = "data/h5nx/metadata-with-clade.tsv",
        reference = "config/reference_sequence_{segment}_A_goose_CR_2021.gb", #H3N8 from 1997
        vaccine_strains = "config/vaccine_strains.json"

files = rules.files.params

rule unzip_h5nx:
    input:
        "h5-data-updates/h5nx.zip"
    output:
        expand("data/h5nx/{segment}/sequences.fasta", segment=SEGMENTS),
        files.input_metadata
    shell:
        """
            unzip -o h5-data-updates/h5nx.zip -d data/
        """

rule extract_excel:
    input:
        "data/Phenotypic characterizations.xlsx"
    output:
        phenotypes_tsv = 'data/phenotypes.tsv',
        sources_tsv = 'data/sources.tsv'
    run:
        split_phenotypes_excel(input[0], output.phenotypes_tsv, output.sources_tsv)

rule genoflu_dataflow:
    input:
        expand("data/h5nx/{segment}/sequences.fasta", segment=SEGMENTS)
    output:
        expand("data/genoflu/{segment}.fasta", segment=SEGMENTS)
    run:
        genoflu_dataflow()

rule genoflu_run:
    input:
        rules.genoflu_dataflow.output
    output:
        'data/genoflu/results/results.tsv'
    shell:
        '''
            # this avoids a quirk of the GenoFlu package... avoids UnboundLocalError related to excel_stats
            rm -rf data/genoflu/temp/
            python GenoFLU-multi/bin/genoflu-multi.py -n 12 -f data/genoflu
        '''

rule genoflu_postprocess:
    input:
        metadata=files.input_metadata,
        genoflu=rules.genoflu_run.output[0]
    output:
        metadata='data/metadata.tsv',
        counts='data/genoflu/results/counts.tsv'
    run:
        genoflu_postprocess(
            input.metadata, input.genoflu, output.metadata, output.counts,
            NUMBER_OF_GENOTYPES, GENOTYPES_TO_INCLUDE
        )

def min_length(w):
    len_dict = {"pb2": 2100, "pb1": 2100, "pa": 2000, "ha":1600, "np":1400, "na":1270, "mp":900, "ns":800}
    length = len_dict[w.segment]
    return(length)

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - excluding strains in {input.exclude}
          - samples with missing region and country metadata
          - excluding strains prior to {params.min_date}
        """
    input:
        sequences = "data/h5nx/{segment}/sequences.fasta",
        metadata = rules.genoflu_postprocess.output.metadata,
        include = "config/include_strains.txt",
        exclude = "config/exclude_strains.txt"
    output:
        sequences = "data/results/filtered_{segment}.fasta"
    params:
        group_by = "month host region genoflu_bin", #month host location
        sequences_per_group = 15,
        min_date = 2021,
        min_length = min_length,  # instead of specifying one parameter value, we can use a function to specify minimum lengths that are unique to each segment
        exclude_where = "host=laboratoryderived host=ferret host=unknown host=other host=host country=? region=?"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --include {input.include} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --exclude-where {params.exclude_where} \
            --min-length {params.min_length} \
            --non-nucleotide
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "data/results/aligned_{segment}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference \
            --nthreads 1
        """


rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "data/results/tree-raw_{segment}.nwk"
    params:
        method = "iqtree"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method} \
            --nthreads 1
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.genoflu_postprocess.output.metadata
    output:
        tree = "data/results/tree_{segment}.nwk",
        node_data = "data/results/branch-lengths_{segment}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "data/results/nt-muts_{segment}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}\
            --keep-ambiguous
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "data/results/aa-muts_{segment}.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.genoflu_postprocess.output.metadata
    output:
        node_data = "data/results/traits_{segment}.json",
    params:
        columns = lambda wc: f"host region country division flyway Domestic_Status genoflu_bin Genotype genoflu_{wc.segment}_lineage",
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule configs:
    input:
        "config/auspice_config.json"
    output:
        "config/auspice_{segment}_config.json"
    params:
        lambda wildcards: f'GenoFlu {wildcards.segment.upper()} lineage'
    shell:
        '''
            jq '.colorings += [
                {{
                    "key": "genoflu_{wildcards.segment}_lineage",
                    "title": "{params}",
                    "type": "categorical"
                }}
            ]' {input} > {output}
        '''

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.genoflu_postprocess.output.metadata,
        node_data = [rules.refine.output.node_data,rules.traits.output.node_data,rules.ancestral.output.node_data,rules.translate.output.node_data,files.vaccine_strains],
        auspice_config = rules.configs.output[0],
        colors = "config/colors.tsv",
        description = "config/description.md",
    output:
        auspice_json = "data/auspice/h5nx_{segment}.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data}\
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --include-root-sequence \
            --colors {input.colors} \
            --output {output.auspice_json}
        """

rule annotate:
    input:
        auspice = rules.export.output.auspice_json,
        phenotypes = rules.extract_excel.output.phenotypes_tsv,
        sources = rules.extract_excel.output.sources_tsv
    output:
        "data/ml/h5nx-{segment}-ceirr.json"
    run:
        phenotypic_characterization_annotation(
            input.auspice, input.phenotypes, input.sources, output[0]
        )

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
