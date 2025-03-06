# MRI acquisition and reconstruction cookbook: recipes for reproducibility, served with real-world flavour

This repository contains the code to reproduce the results of the article **MRI acquisition and reconstruction cookbook: recipes for reproducibility, served with real-world flavour*.*

## Requirements

This repository has various requirements to run the experiments. Key requirements are:
+ Fig 3/4:
  + pypulseq (version 1.4.2)
  + MRzeroCore (version 0.3.12)
  + A conda environment file is provided in `environment.yml` to install all required packages.
+ Fig 6/7:
  + bart (version 0.9.00) c.f. https://github.com/mrirecon/bart
  + view (version 0.3.00) c.f. https://github.com/mrirecon/view
  + For interaction of BART and python (Fig 7), the BART_TOOLBOX_PATH environment variable must be set to the path of the BART toolbox.

For composing of figures to final PDFs, inkscape is required.

The code was tested on Debian 12.

## Data

k-Space data to reproduce the reconstructions (Figures 6 and 7) is hosted on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14497769.svg)](https://doi.org/10.5281/zenodo.14497769).
The data can be downloaded with the script in `Data/download_data.sh`.

## Reproducing
All experiments can be reproduced by running
```bash
export BART_TOOLBOX_PATH=/path/to/bart/
./run_all.sh
```

## How to cite

Tamir, J.I., Blumenthal, M., Wang, J. et al. MRI acquisition and reconstruction cookbook: recipes for reproducibility, served with real-world flavour. Magn Reson Mater Phy (2025). https://doi.org/10.1007/s10334-025-01236-4

