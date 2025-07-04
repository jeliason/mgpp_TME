# MGPP modeling for Tissue Data

This repository contains code and data for the paper:

*Investigating Ecological Interactions in the Tumor Microenvironment using Joint Species Distribution Models for Point Patterns* (Eliason & Rao 2025)

We present a statistically rigorous framework for analyzing spatial cell–cell interactions in the tumor microenvironment (TME) using multitype Gibbs point process (MGPP) models. These models are implemented as joint species distribution models (JSDMs) to quantify spatial attraction and repulsion among annotated cell types in multiplexed colorectal cancer images.

---

## Purpose

This repository enables full reproducibility of the analyses and figures in the manuscript. It includes:

- Preprocessed CODEX multiplexed imaging data for 35 colorectal cancer patients.
- MGPP model fitting functions using a fork of the `ppjsdm` R package.
- An annotated Quarto notebook reproducing all figures and tables.
- An HPC-compatible batch script for parallel model fitting across patients.

---

## Getting Started

### Prerequisites

- R version ≥ 4.2  
- Packages: `tidyverse`, `spatstat.core`, `future`, `future.apply`, `remotes`, `quarto`

### Installation

Clone the repository and install required packages:

```bash
git clone https://github.com/jeliason/mgpp_TME.git
cd mgpp_TME
# Install required R packages
Rscript -e "install.packages(c('remotes', 'quarto', 'tidyverse', 'spatstat.core', 'future', 'future.apply'))"
Rscript -e "remotes::install_github('jeliason/ppjsdm')"
```

---

## Example Analysis

To reproduce the full example workflow and regenerate all manuscript figures and tables:

```r
quarto::quarto_render("mgpp_example_workflow.qmd")
```

The HTML output will be saved in `_output/`.

---

## Full Model Fitting on HPC

To run the full MGPP model across all patients in parallel:

1. Modify `hpc.R` for your computing environment (paths).
2. Submit the job on Slurm using `rslurm`:

```bash
Rscript hpc.R
```

---

## Citation

Please cite the manuscript if you use this code or data:

```bibtex
@article{eliason2024,
  title = {Investigating {{Ecological Interactions}} in the {{Tumor Microenvironment Using Joint Species Distribution Models}} for {{Point Patterns}}},
  author = {Eliason, Joel and Rao, Arvind},
  year = {2024},
  month = jun,
  journal = {The New England Journal of Statistics in Data Science},
  volume = {2},
  number = {3},
  pages = {296--310},
  publisher = {New England Statistical Society},
  issn = {2693-7166},
  doi = {10.51387/24-NEJSDS66},
}

```

---

## License

This repository is released under the MIT License. See LICENSE for details.

---

## Contact

For questions, feedback, or bug reports, please open an issue or contact:

**Joel Eliason**  
joelne@umich.edu
