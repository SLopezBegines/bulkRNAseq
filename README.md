# bulkRNAseq

> Short description of the project.

## Structure

```
bulkRNAseq/
├── rawdata/        # Immutable raw input data
├── data/           # Processed/intermediate data
├── code/           # Scripts and modules
│   └── utils/      # Shared helper functions
├── notebooks/      # Exploratory notebooks
│   ├── output/     # Rendered reports
│   └── old/        # Archived notebooks
├── output/         # Final results
│   ├── figures/
│   ├── tables/
│   └── reports/
├── docs/           # Documentation and notes
├── logs/           # Pipeline logs (git-ignored)
├── config/         # Parameters and config files
└── envs/           # Conda/pip environment files
```

## Setup

```bash
conda env create -f envs/environment.yml
conda activate bulkRNAseq
```

## Usage

_Describe how to run the analysis here._
