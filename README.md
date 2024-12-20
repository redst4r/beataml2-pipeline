# beataml2-pipeline

[![Powered by Kedro](https://img.shields.io/badge/powered_by-kedro-ffc900?logo=kedro)](https://kedro.org)

## Overview

A simple Kedro pipeline to turn the [BeatAML multiomics datasets](https://biodev.github.io/BeatAML2/) into a single [muon](https://muon.scverse.org/) object.
Beware that processing might be done a little bit different than in the original publication and [code](https://github.com/biodev/beataml2_manuscript):
In particular, we use the `labId` as the unique identfier for a sample. 
- A patient might have multiple `labIds`: Collection of PBMCs and/or Bone marrow, or multiple timepoints (initial diagnosis and relapse)
- A `labId` has at maximum one DNAseq, one RNAseq and one drug-screen done.
In the final `muon` object, each row corresponds to a `labId` (see `.obs.index`).

## Usage
Kedro is a pipeline-framework that simply turns input files (the excel/txt files on https://biodev.github.io/BeatAML2/) into some ouput (`.h5mu`) via a series of transformations.
This pipeline either 
- pulls the raw files from the website each time you run the pipeline,
- or you can download the files manually and advise kedro to use those local files.

### Kedro setup:
All required dependencies are listed in the `pyproject.toml` (including kedro), just install those into a virtual environment with your favorite package manager. Make sure to activate the environment.

Using `uv` it's a simple as
```bash
uv run  # just a dummy run to create a venv in `.venv` and install all packages
source .venv/bin/activate 
```
### Remote files
Probably the easiest way to get going: This simply pulls the required files of the web (see the specifications in `conf/base/catalog.yml`)
```bash
kedro run -e base
```
This will pull the raw data from the web, and store the final `.h5mu` in `data/03_primary/beataml2.h5mu`.

*Note*: This can get time-consuming if the pipeline is run multiple times (downloading ~500Mb each time).
Use the `manual download` option in those cases (which persists those files to disk)

### Manually downloaded
A bit more manual work: Download all the files from [this website](https://biodev.github.io/BeatAML2/) and [this spreadsheet from PMC](https://pmc.ncbi.nlm.nih.gov/articles/instance/6280667/bin/NIHMS1504008-supplement-Supplementary_Tables_S1-S22.xlsx) into `data/01_raw`, making sure that the names are correct (compare with `conf/local/catalog.yaml`).

*Note*: this can be semi-automated, running the script `download_data.sh`.

```bash
kedro run -e local
```
This will use manually downloaded data on disk (in `data/01_raw/`), and store the final `.h5mu` in `data/03_primary/beataml2.h5mu`

### Analysing the h5mu
```python
import muon as mu
mu.read('data/03_primary/beataml2.h5mu')
```

## TODO
- [ ] Nov2024: In the original code, they overwrite some of the WES-data with labs done for specific mutations. E.g. even though the sequencing data didnt show any NPM1 mutation, the clinical data might have NPM1 annotated as mutated (via some other assay). They'd correct the `sample x mutation` matrix to include that mutation (changing the entry from 0 to 1)

