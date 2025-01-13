import pandas as pd


phenotypes_excel = 'data/Phenotypic characterizations.xlsx'
phenotypes_tsv = 'data/phenotypes.tsv'
sources_tsv = 'data/sources.tsv'

pd.read_excel(phenotypes_excel).to_csv(
    phenotypes_tsv, sep='\t', index=False
)
pd.read_excel(
    phenotypes_excel, sheet_name='Sources'
).to_csv(
    sources_tsv, sep='\t', index=False
)
