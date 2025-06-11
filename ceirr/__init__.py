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

        auspice['meta']['colorings'] = [{
            'key': 'phenotypic_characterization',
            'title': 'Phenotypic characterization',
            'type': 'categorical'
        },
        {
            'key': 'invivo_characterization',
            'title': 'In vivo phenotypic characterization',
            'type': 'categorical'
        },
        {
            'key': 'invitro_characterization',
            'title': 'In vitro phenotypic characterization',
            'type': 'categorical'
        },
        {
            'key': 'antiviral_characterization',
            'title': 'Antiviral phenotypic characterization',
            'type': 'categorical'
        }] + auspice['meta']['colorings']

        with open(output_path, 'w') as json_file:
            json.dump(auspice, json_file, indent=2)


