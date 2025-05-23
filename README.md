# nextstrain-ceirr

Nextstrain pipelines for [CEIRR](https://www.ceirr-network.org/).

## Usage

Usage instructions assume you've successfully followed the installation instructions and have a basic understanding of [Nextstrain](https://nextstrain.org/).

Run:

```
nextstrain build . -j $NUMBER_OF_JOBS all
```

View:

```
nextstrain view data/ml
```


## Installation

[Install and configure Bioconda](https://bioconda.github.io/).

Install dependencies:

```
conda create -n nextstrain-ceirr pandas biopython blast openpyxl snakemake nextstrain jq
```

Contact Stephen for access to private data. Place phenotypic spreadsheet at `data/Phenotypic characterizations.xlsx`.

Pull down this repository and dependent repositories:

```
git clone --recurse-submodules https://github.com/moncla-lab/nextstrain-ceirr
```

### Moncla lab instructions 
If you are a member of the Moncla lab, you should have access to the `h5-data-updates` repo and the phenotypes spreadsheet. If so, do the following: 