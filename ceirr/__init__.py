import json
import csv

from Bio import SeqIO
import pandas as pd


SEGMENTS = ["pb1", "pb2", "na", "pa", "ha", "np", "mp", "ns"]

# Column name for strain in new spreadsheet format
STRAIN_COLUMN = "Virus strain"

# Column name for source string in new spreadsheet format
SOURCE_COLUMN = "Source Preference/Publication (PMID or PMC #) or CEIRR Center"

# Mapping of phenotype columns to characterization categories
NEWCOL_TO_PCCLASS = {
    "Receptor Binding Preference": "invitro",
    "Receptor Binding Preference - Method ": "invitro",
    "HA hemagglutination activity": "invitro",
    "pH of Inactivation/ pH of fusion (Please define the type of assay used and any comparative viruses used in the assay)": "invitro",
    "NA activity (please include assay type/substrate used and the comparative virus used).": "invitro",
    "Replication capacity in relevant three-dimensional airway cultures. ": "invitro",
    "Antiviral Sensitivity": "antiviral",
    "Pathogenesis - Animal Species Tested": "invivo",
    "Pathological Findings - Experimental design/details": "invivo",
    "Transmission - Animal Species Tested": "invivo",
    "Transmission Findings - Experimental design/details": "invivo",
}


def split_phenotypes_excel(input_excel, phenotypes_output):
    """Extract the CEIRR RAP H5 sheet from the Excel file to TSV."""
    pd.read_excel(input_excel, sheet_name="CEIRR RAP H5").to_csv(
        phenotypes_output, sep="\t", index=False
    )


def clean_key(key):
    """Clean strain names by removing invisible characters."""
    if key is None or (isinstance(key, float) and pd.isna(key)):
        return ""
    return str(key).replace("\u200b", "").replace("\xa0", "").strip()


def load_strain_crossref(crossref_path):
    """
    Load strain cross-reference table.
    Returns dict mapping google_sheet_metadata_strain -> ml_database_strain
    """
    crossref = {}
    with open(crossref_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sheet_strain = clean_key(row["google_sheet_metadata_strain"])
            db_strain = clean_key(row["ml_database_strain"])
            if sheet_strain and db_strain:
                crossref[sheet_strain] = db_strain
    return crossref


def load_source_strings(source_strings_path):
    """
    Load source string mapping table.
    Returns dict mapping source_string -> list of source_ids
    """
    source_strings = {}
    with open(source_strings_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            source_string = row["source_string"]
            ids_str = row.get("source_ids", "")
            if source_string and ids_str:
                source_ids = [int(x.strip()) for x in ids_str.split(",") if x.strip()]
                source_strings[source_string] = source_ids
    return source_strings


def load_sources(sources_path):
    """
    Load source details table.
    Returns dict mapping source_id -> {name, url}
    """
    sources = {}
    with open(sources_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            source_id = int(row["source_id"])
            sources[source_id] = {
                "name": row["name"],
                "url": row.get("url", "")
            }
    return sources


def validate_strain_crossref(spreadsheet_strains, crossref):
    """
    Validate strain cross-reference and print warnings.
    Returns set of matched spreadsheet strains.
    """
    spreadsheet_set = set(spreadsheet_strains)
    crossref_set = set(crossref.keys())

    matched = spreadsheet_set & crossref_set
    unmatched = spreadsheet_set - crossref_set
    stale = crossref_set - spreadsheet_set

    print("\n=== Strain Cross-Reference Validation ===")
    print(f"MATCHED ({len(matched)}/{len(spreadsheet_set)}):")
    for strain in sorted(matched)[:5]:
        print(f"  {strain} → {crossref[strain]}")
    if len(matched) > 5:
        print(f"  ... and {len(matched) - 5} more")

    if unmatched:
        print(f"\nUNMATCHED ({len(unmatched)}/{len(spreadsheet_set)}) - add to maintenance_data/strain_crossref.tsv:")
        for strain in sorted(unmatched):
            print(f"  {strain}")

    if stale:
        print(f"\nSTALE ({len(stale)} entries) - in crossref but not in spreadsheet:")
        for strain in sorted(stale):
            print(f"  {strain}")

    print()
    return matched


def validate_source_strings(spreadsheet_sources, source_strings_map):
    """
    Validate source string mappings and print warnings.
    Returns set of matched source strings.
    """
    spreadsheet_set = set(s for s in spreadsheet_sources if s and str(s) != "nan")
    mapping_set = set(source_strings_map.keys())

    matched = spreadsheet_set & mapping_set
    unmatched = spreadsheet_set - mapping_set
    stale = mapping_set - spreadsheet_set

    print("=== Source String Validation ===")
    print(f"MATCHED ({len(matched)}/{len(spreadsheet_set)}):")
    for source in sorted(matched)[:3]:
        print(f"  {source[:60]}{'...' if len(source) > 60 else ''}")
    if len(matched) > 3:
        print(f"  ... and {len(matched) - 3} more")

    if unmatched:
        print(f"\nUNMATCHED ({len(unmatched)}/{len(spreadsheet_set)}) - add to maintenance_data/source_strings.tsv:")
        for source in sorted(unmatched):
            print(f"  {source}")

    if stale:
        print(f"\nSTALE ({len(stale)} entries) - in source_strings.tsv but not in spreadsheet:")
        for source in sorted(stale):
            print(f"  {source[:60]}{'...' if len(source) > 60 else ''}")

    print()
    return matched


def phenotypic_characterization_annotation(
    auspice_path,
    phenotypes_path,
    strain_crossref_path,
    source_strings_path,
    sources_path,
    output_path
):
    """
    Annotate Auspice JSON with phenotypic characterization data.

    Uses cross-reference tables for strain matching and source citation.
    """
    # Load Auspice JSON
    with open(auspice_path) as auspice_file:
        auspice = json.load(auspice_file)

    # Load phenotypes data
    with open(phenotypes_path) as phenotypes_file:
        phenotypes = list(csv.DictReader(phenotypes_file, delimiter="\t"))

    # Load cross-reference tables
    crossref = load_strain_crossref(strain_crossref_path)
    source_strings_map = load_source_strings(source_strings_path)
    sources_map = load_sources(sources_path)

    # Get all strains from spreadsheet
    spreadsheet_strains = [clean_key(row.get(STRAIN_COLUMN, "")) for row in phenotypes]
    spreadsheet_strains = [s for s in spreadsheet_strains if s]

    # Get all source strings from spreadsheet
    spreadsheet_sources = [row.get(SOURCE_COLUMN, "") for row in phenotypes]

    # Validate and print warnings
    validate_strain_crossref(spreadsheet_strains, crossref)
    validate_source_strings(spreadsheet_sources, source_strings_map)

    # Build annotation dict keyed by ml_database_strain
    annotation_dict = {}
    for row in phenotypes:
        sheet_strain = clean_key(row.get(STRAIN_COLUMN, ""))
        if sheet_strain and sheet_strain in crossref:
            db_strain = crossref[sheet_strain]
            annotation_dict[db_strain] = row

    has_annotation = set(annotation_dict.keys())
    valid_keys = set(NEWCOL_TO_PCCLASS.keys())

    ordinals = [
        "first", "second", "third", "fourth", "fifth",
        "sixth", "seventh", "eighth", "ninth", "tenth",
    ]

    def traverse(node):
        pc_attrs = {
            "phenotypic_characterization": {"value": "No"},
            "invivo_characterization": {"value": "No"},
            "invitro_characterization": {"value": "No"},
            "antiviral_characterization": {"value": "No"},
        }
        node_attrs = {}
        strain = node["name"]

        if strain in has_annotation:
            strain_annotations = annotation_dict[strain]

            # Process phenotype columns
            for key, value in strain_annotations.items():
                if key not in valid_keys:
                    continue
                if value == "" or value is None or value == "\xa0":
                    continue
                if str(value) == "nan":
                    continue

                pc = NEWCOL_TO_PCCLASS[key]
                pc_key = f"{pc}_characterization"
                node_attrs[key] = {"value": value}
                pc_attrs["phenotypic_characterization"]["value"] = "Yes"
                pc_attrs[pc_key]["value"] = "Yes"

            # Process source citations
            source_string = strain_annotations.get(SOURCE_COLUMN, "")
            if source_string and str(source_string) != "nan":
                source_ids = source_strings_map.get(source_string, [])
                for i, sid in enumerate(source_ids):
                    if sid in sources_map:
                        source = sources_map[sid]
                        key = ordinals[i].capitalize() + " source"
                        node_attrs[key] = {"value": source["name"]}
                        if source["url"]:
                            node_attrs[key]["url"] = source["url"]

        node["node_attrs"].update(node_attrs)
        node["node_attrs"].update(pc_attrs)

        # Recursively traverse children
        for child in node.get("children", []):
            traverse(child)

    traverse(auspice["tree"])

    # Add colorings for phenotypic characterization categories
    auspice["meta"]["colorings"] = [
        {
            "key": "phenotypic_characterization",
            "title": "Phenotypic characterization",
            "type": "categorical",
        },
        {
            "key": "invivo_characterization",
            "title": "In vivo phenotypic characterization",
            "type": "categorical",
        },
        {
            "key": "invitro_characterization",
            "title": "In vitro phenotypic characterization",
            "type": "categorical",
        },
        {
            "key": "antiviral_characterization",
            "title": "Antiviral phenotypic characterization",
            "type": "categorical",
        },
    ] + auspice["meta"]["colorings"]

    with open(output_path, "w") as json_file:
        json.dump(auspice, json_file, indent=2)


def create_segment_config(
    input_config_path, url_file_path, output_config_path, segment, ceirr_url_order
):
    """
    Create segment-specific auspice config with CEIRR URL and GenoFlu coloring.

    Takes as input config path, output config path, and information on the ceirr url.
    This includes the url location of the fillable ceirr spreadsheet (url_file_path),
    and the ordering of the url in the auspice config. The ceirr_url_order should be
    set to a 0-indexed integer that relates which maintainer link ceirr_url_order
    should reference.
    """
    # Read the base config
    with open(input_config_path) as f:
        config = json.load(f)

    # Read the URL from file and add to third maintainer
    try:
        with open(url_file_path) as f:
            ceirr_url = f.read().strip()
        if len(config.get("maintainers", [])) > 2:
            config["maintainers"][ceirr_url_order]["url"] = ceirr_url
    except FileNotFoundError:
        # If URL file doesn't exist, continue without adding URL
        pass

    # Add the segment-specific genoflu coloring
    segment_key = segment.upper()
    genoflu_coloring = {
        "key": f"genoflu_{segment_key}",
        "title": f"GenoFlu {segment_key} lineage",
        "type": "categorical",
    }

    config["colorings"] = config.get("colorings", []) + [genoflu_coloring]

    # Write the modified config
    with open(output_config_path, "w") as f:
        json.dump(config, f, indent=2)
