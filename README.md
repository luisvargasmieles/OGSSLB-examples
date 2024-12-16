# OGSSLB-examples

This repository contains R code to reproduce the results presented in the paper *"Outcome-Guided Spike-and-Slab Lasso Biclustering: A Novel Approach for Enhancing Biclustering Techniques for Gene Expression Analysis"* (Vargas-Mieles, Kirk, Wallace, 2024).

The `OGSSLB` R package is available on GitHub:  
[OGSSLB GitHub Repository](https://github.com/luisvargasmieles/OGSSLB).

## Instructions

Follow these steps to run the code:

1. Install the `SSLB` and `OGSSLB` R packages using the `devtools` package. Run the following commands in your R environment:

    ```R
    install.packages("devtools")
    library(devtools)
    install_github("gemoran/SSLB")
    install_github("luisvargasmieles/OGSSLB")
    ```

2. This repository includes examples comparing OGSSLB (Vargas-Mieles et al., 2024) to the SSLB method (Moran et al., 2021).

## Repository Structure

- **`OGSSLB-functions.R`**: Contains helper functions required for all the R scripts.
- **`sim1`**: This folder includes code for reproducing the simulation study:
    - **`OGSSLB-sim1.R`**: Runs the SSLB and OGSSLB biclustering methods and plots Consensus Scores.
- **`ImmuNexUT`**: This folder includes:
    1. Code for processing the dataset from Ota et al. (2021).
    2. Scripts for reproducing the results from the Immune Cell Gene Expression Atlas (University of Tokyo) numerical experiment section, as shown in Vargas-Mieles et al. (2024).

## References

1. Luis A. Vargas-Mieles, Paul D. W. Kirk, Chris Wallace, *"Outcome-Guided Spike-and-Slab Lasso Biclustering: A Novel Approach for Enhancing Biclustering Techniques for Gene Expression Analysis"*, arXiv preprint arXiv:2412.08416, 2024.
2. Gemma E. Moran, Veronika Ročková, Edward I. George, *"Spike-and-Slab Lasso Biclustering"*, The Annals of Applied Statistics, 15(1), 148-173, 2021.
3. Mineto Ota et al., *"Dynamic landscape of immune cell-specific gene regulation in immune-mediated diseases"*, Cell, 184(11), 3006-3021, 2021.
