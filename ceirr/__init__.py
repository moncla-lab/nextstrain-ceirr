import json
import csv

from Bio import SeqIO
import pandas as pd


SEGMENTS = ["pb1","pb2","na","pa","ha","np","mp","ns"]


def split_phenotypes_excel(
    input_excel, phenotypes_output, sources_output
):
    pd.read_excel(input_excel).to_csv(
        phenotypes_output, sep='\t', index=False
    )
    pd.read_excel(
        input_excel, sheet_name='Sources'
    ).to_csv(
        sources_output, sep='\t', index=False
    )



def phenotypic_characterization_annotation(
    auspice_path, phenotypes_path, sources_path, output_path
):
        with open(auspice_path) as auspice_file:
            auspice = json.load(auspice_file)
        with open(phenotypes_path) as phenotypes_file:
            phenotypes = list(csv.DictReader(phenotypes_file, delimiter='\t'))
        with open(sources_path) as sources_file:
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

        with open(output_path, 'w') as json_file:
            json.dump(auspice, json_file, indent=2)


def genoflu_dataflow():
    seq_dicts = {}
    all_headers = set()

    for segment in SEGMENTS:
        fasta_file = f'data/h5nx/{segment}/sequences.fasta'
        seq_dict = {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}
        seq_dicts[segment] = seq_dict
        if not all_headers:
            all_headers = set(seq_dict.keys())
        else:
            all_headers.intersection_update(seq_dict.keys())

    sorted_headers = sorted(all_headers)

    for segment in SEGMENTS:
        output_file = f"data/genoflu/{segment}.fasta"
        with open(output_file, "w") as out_f:
            for header in sorted_headers:
                SeqIO.write(seq_dicts[segment][header], out_f, "fasta")


def parse_genoflu_genotypes_list(annotation):
    result = {segment: None for segment in SEGMENTS}
    if pd.isna(annotation):
        return result
    try:
        for entry in annotation.split(', '):
            genoflu_gene_key, value = entry.split(':')
            ml_gene_key = genoflu_gene_key.strip().lower()
            if ml_gene_key in result:
                result[ml_gene_key] = value.strip()
    except:
        import pdb; pdb.set_trace()

    return result


def genoflu_postprocess(
        input_tsv, genoflu_tsv, counts_tsv, number_of_genotypes=9, included_genotypes=[]
    ):
    df = pd.read_csv(input_tsv, sep='\t')
    df['was_assigned'] = df.Genotype.apply(
        lambda col: col[0:min(len('Not assigned:'), len(col))] != 'Not assigned:'
    )
    df.loc[~df.was_assigned, 'Genotype'] = 'Unassigned'
    counts = df['Genotype'].value_counts()
    print('Top Genoflu genotypes:', counts.to_string())
    number_to_bin = number_of_genotypes - len(included_genotypes)
    parsed_genotype_list = [
        parse_genoflu_genotypes_list(row)
        for row in df['Genotype List Used >=98%']
    ]
    for segment in SEGMENTS:
        df[f'Genoflu {segment} variant'] = [
            row[f'{segment}'] 
            for row in parsed_genotype_list 
        ]
    top = counts.head(number_to_bin).index
    desired_genotypes = list(top) + included_genotypes
    df["genoflu_bin"] = df["Genotype"].where(
        df["Genotype"].isin(desired_genotypes), "Not dominant genotype"
    )
    df.rename(columns={"Genotype List Used, >=98%": "Genotype List Used >=98%"}, inplace=True)
    counts.to_csv(counts_tsv, sep="\t", header=False)
    df.to_csv(genoflu_tsv, index=False, sep='\t')


def merge_metadata(input_metadata, input_genoflu, output_tsv):
    metadata_df = pd.read_csv(input_metadata, sep='\t')
    genoflu_df = pd.read_csv(input_genoflu, sep='\t')
    merged_df = metadata_df.merge(
        genoflu_df, left_on="strain", right_on="Strain", how="left"
    )
    merged_df.to_csv(output_tsv, index=False, sep='\t')