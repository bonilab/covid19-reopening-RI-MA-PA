# Public repository for the Boni modeling group's paper: "SARS-CoV-2 epidemic after social and economic reopening in three US states reveals shifts in age structure and clinical characteristics."
*Authors: Nathan B Wikle, Thu Nguyen-Anh Tran, Ephraim M Hanks, Maciej F Boni, and 12 others.*

This repository contains source code and additional supplementary materials from our manuscript, "SARS-CoV-2 epidemic after social and economic reopening in three US states reveals shifts in age structure and clinical characteristics." This manuscript has been accepted for publication in *Science Advances*; a preprint is available on [medRxiv](https://www.medrxiv.org/content/10.1101/2020.11.17.20232918v3).

The following instructions provide details on how to run the source code underlying the analysis, including replication of the main figures and results.

## Software Requirements

### R 

The code has been tested with R version 3.6.0 (2019-04-26) -- "Planting of a Tree."  The following R packages (version number in parentheses) must be installed before the code will run successfully:

- `mc2d` (0.1.18)
- `mvtnorm` (1.1.1)
- `readxl` (1.3.1)
- `fda` (5.1.9)
- `snow` (0.4.3)
- `doParallel` (1.0.16)
- `foreach` (1.5.1)
- `msm` (1.6.8)
- `splines2` (0.4.1)
- `reticulate` (1.18)
- `doMC` (1.3.7)
- `extrafont` (0.17)
- `scales` (1.1.1)
- `dplyr` (1.0.7)
- `forcats` (0.5.1)
- `viridis` (0.6.1)
- `hrbrthemes` (0.8.0)
- `ggplot2` (3.3.5)
- `ggridges` (0.5.2)

In addition, the Ubunto Mono font family is used in the production of several figures. It can be downloaded [here](https://fonts.google.com/specimen/Ubuntu+Mono).

### C++


## Instructions

Before running any code, make sure the required R packages have been installed. [...]

