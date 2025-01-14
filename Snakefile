import csv
import json


SEGMENTS = ["pb1","pb2","na","pa","ha","np","mp","ns"]

rule all:
    input:
        auspice_json = expand("data/ml/h5nx_{segment}.json", segment=SEGMENTS)

rule unzip_h5nx:
    input:
        "h5-data-updates/h5nx.zip"
    output:
        expand("data/h5nx/{segment}/sequences.fasta", segment=SEGMENTS)
    shell:
        """
            unzip -o h5-data-updates/h5nx.zip -d data/
        """

rule files:
    params:
        input_metadata = "data/h5nx/metadata.tsv",
        reference = "config/reference_sequence_{segment}_A_goose_CR_2021.gb", #H3N8 from 1997
        vaccine_strains = "config/vaccine_strains.json"

files = rules.files.params

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
        metadata = files.input_metadata,
        include = "config/include_strains.txt",
        exclude = "config/exclude_strains.txt"
    output:
        sequences = "data/results/filtered_{segment}.fasta"
    params:
        group_by = "month host region", #month host location
        sequences_per_group = 25, #test changing from 25 to 2
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
        metadata = files.input_metadata
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
        metadata = files.input_metadata
    output:
        node_data = "data/results/traits_{segment}.json",
    params:
        columns = "host region country division flyway Domestic_Status",
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = files.input_metadata,
        node_data = [rules.refine.output.node_data,rules.traits.output.node_data,rules.ancestral.output.node_data,rules.translate.output.node_data,files.vaccine_strains],
        auspice_config = "config/auspice_config.json",
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
        phenotypes = 'data/phenotypes.tsv',
        sources = 'data/sources.tsv'
    output:
        "data/ml/h5nx_{segment}.json"
    run:
        with open(input.auspice) as auspice_file:
            auspice = json.load(auspice_file)
        with open(input.phenotypes) as phenotypes_file:
            phenotypes = list(csv.DictReader(phenotypes_file, delimiter='\t'))
        with open(input.sources) as sources_file:
            sources = list(csv.DictReader(sources_file, delimiter='\t'))
        def clean_key(key):
            return key.replace("\u200b", "").replace("\xa0", "").strip()
        annotation_dict = {
            clean_key(row['Nextstrain Strain']): row
            for row in phenotypes 
            if row != ''
        }
        has_annotation = set(annotation_dict.keys())
        newcol_to_pcclass = {
            'Receptor binding': 'invitro',
            'Human airway replication': 'invitro',
            'pH of Inactivation/ pH of fusion': 'invitro',
            'Antiviral sensitvity': 'antiviral',
            'Mouse pathogenesis': 'invivo',
            'Ferret pathogenesis': 'invivo',
            'Ferret transmission':  'invivo',
            'Swine Pathogenesis': 'invivo',
            'Swine Transmission':  'invivo'
        }
        valid_keys = set(newcol_to_pcclass.keys())
        ordinals = [
            'first', 'second', 'third', 'fourth', 'fifth',
            'sixth', 'seventh', 'eighth', 'ninth', 'tenth'
        ]
        def traverse(node):
            pc_attrs = {
                'phenotypic_characterization': {'value': 'No'},
                'invivo_characterization': {'value': 'No'},
                'invitro_characterization': {'value': 'No'},
                'antiviral_characterization': {'value': 'No'},
            }
            node_attrs = {}
            strain = node['name']
            if strain in has_annotation:
                strain_annotations = annotation_dict[strain]
                for key, value in strain_annotations.items():
                    if not key in valid_keys:
                        continue
                    if value == '' or value is None or value == '\xa0':
                        continue
                        
                    pc = newcol_to_pcclass[key]
                    pc_key = f'{pc}_characterization'
                    if value != None and str(value) != "nan":
                        node_attrs[key] = {'value': value}
                        pc_attrs['phenotypic_characterization']['value'] = 'Yes'
                        pc_attrs[pc_key]['value'] = 'Yes'
                        
                for i, note_str in enumerate(strain_annotations['SourcesByIndex'].split(',')):
                    note_index = int(note_str)
                    key = ordinals[i].capitalize() + ' source'
                    node_attrs[key] = {'value': sources[note_index]['Title']}
                    if sources[note_index]['URL'] != '':
                        node_attrs[key]['url'] = sources[note_index]['URL']

                    
            node['node_attrs'].update(node_attrs)
            node['node_attrs'].update(pc_attrs)
            # Recursively traverse each child if 'children' is present and is a list
            for child in node.get('children', []):
                traverse(child)
        traverse(auspice['tree'])

        auspice['meta']['colorings'].append({
            'key': 'phenotypic_characterization',
            'title': 'Phenotypic characterization',
            'type': 'categorical'
        })
        auspice['meta']['colorings'].append({
            'key': 'invivo_characterization',
            'title': 'In vivo phenotypic characterization',
            'type': 'categorical'
        })
        auspice['meta']['colorings'].append({
            'key': 'invitro_characterization',
            'title': 'In vitro phenotypic characterization',
            'type': 'categorical'
        })
        auspice['meta']['colorings'].append({
            'key': 'antiviral_characterization',
            'title': 'Antiviral phenotypic characterization',
            'type': 'categorical'
        })

        with open(output[0], 'w') as json_file:
            json.dump(auspice, json_file, indent=2)

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
