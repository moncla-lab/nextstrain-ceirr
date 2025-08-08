# nextstrain-ceirr

Nextstrain pipelines for [CEIRR](https://www.ceirr-network.org/).

## Usage

Usage instructions assume you've successfully followed the installation instructions and have a basic understanding of [Nextstrain](https://nextstrain.org/).

### Downloading the metadata

Obtain a link to the phenotypes spreadsheet through the URL in the CEIRR hosted build. Place at `data/Phenotypic characterizations.xlsx`.

To download a copy of the spreadsheet, from the main menu in Excel do: **File > Create a Copy > Download a Copy**.

Run:

```
snakemake -j $NUMBER_OF_JOBS all
```

View:

```
nextstrain view data/ml
```


## Installation

[Install and configure Bioconda](https://bioconda.github.io/).

Install dependencies:

```
conda create -n nextstrain-ceirr pandas biopython blast openpyxl snakemake nextstrain
```

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

## Uploading data

Visit the [CEIRR web application](https://app.ceirr-network.org/) and navigate to the the Nextstrain tab. Privileged users will see options to upload new datasets and new narratives at the bottom of the Auspice page.