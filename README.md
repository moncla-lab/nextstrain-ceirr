# nextstrain-ceirr

Nextstrain pipelines for [CEIRR](https://www.ceirr-network.org/).

## Usage

Usage instructions assume you've successfully followed the installation instructions and have a basic understanding of [Nextstrain](https://nextstrain.org/).

Run:

```
snakemake -j $NUMBER_OF_JOBS all
```

View:

```
auspice view --datasetDir data/ml
```


## Installation

[Install and configure Bioconda](https://bioconda.github.io/).

Install dependencies:

```
conda create -n nextstrain-ceirr pandas biopython blast openpyxl snakemake nextstrain jq
```

Obtain link phenotypes spreadsheet through the link in the CEIRR hosted build. Place at `data/Phenotypic characterizations.xlsx`.

### Downloading the metadata

When viewing the URL, from the main menu do **File > Create a Copy > Download a Copy**.

Pull down this repository and dependent repositories:

```
git clone https://github.com/moncla-lab/nextstrain-ceirr
cd nextstrain-ceirr
git submodule update --init
```

## Submodule Management

### Updating submodules when data changes:
```bash
git submodule update --remote    # then add, commit, and push like usual
```

### Submodule overview for maintainers:
- **h5-data-updates**: Contains shared Genoflu analysis functions. Changes here affect both CEIRR and North America pipelines.
- **GenoFLU-multi**: External tool for influenza genotyping. Update only when new versions are released.
- **nextstrain_hpai_north_america**: Contains logic for preprocessing, running, and processing GenoFlu.