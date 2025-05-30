# README: Code for _Measuring the growth of infectious disease modelling publications and their impact on policymaking: a machine learning review_ manuscript

This repository contains the code and data used to generate the figures and analyses for the peer-reviewed manuscript titled **Measuring the growth of infectious disease modelling publications and their impact on policymaking: a machine learning review**. The code is organized into several R Markdown files, each corresponding to a specific figure or analysis in the manuscript. Below is a guide to help you navigate the repository and understand the structure of the code.

## Repository Structure

### Files and Folders

- **`data/`**: This folder contains the raw data files used in the analysis. The data files are not included in this repository due to size or privacy concerns, but the code assumes the data is stored in this directory.
  
- **`plots/`**: This folder contains the output plots generated by the R Markdown files. The plots are saved as image files (e.g., `.png`, `.jpg`) and are used in the manuscript.

- **`results/`**: This folder contains any intermediate or final results generated by the analysis, such as summary tables or processed data.

- **`Figure 1.Rmd`**: This file contains the code for generating **Figure 1** in the manuscript, which visualizes the timeline of publications, policy citations, and disease outbreaks.

- **`Figure 2B 3.Rmd`**: This file contains the code for generating **Figure 2B and 2C** in the manuscript, which analyze the cross-referencing of policy documents and abstracts by country and income group.

- **`regression.Rmd`**: This file contains the code for performing piecewise Poisson regression analysis on publication and citation data, identifying changepoints in the trends over time.

## Data Sources

The data used in this analysis comes from the following sources:

1. **Publication Data**: The file `490k_141224_1234_dedup.csv` contains metadata for publications, including titles, years, and inclusion/exclusion decisions.
   
2. **Policy Citation Data**: The file `overton_results_expanded.csv` contains data on policy documents that cite infectious disease modeling (IDM) literature, including publication dates and country information.

3. **Outbreak Timeline Data**: The file `Epidemics Timeline.xlsx` contains data on the occurrence of disease outbreaks over time.

4. **World Map Data**: The file `clean_idm_lit_with_locations.csv` contains geographic information for IDM literature, including country codes and locations.

## Code Overview

### Figure 1: Timeline of Publications, Policy Citations, and Outbreaks

- **Objective**: This figure visualizes the trends in publication counts, policy citations, and disease outbreaks over time.
- **Key Steps**:
  - Load and clean publication and policy citation data.
  - Merge the datasets to align the time ranges.
  - Plot the timeline using `ggplot2`, with separate axes for publication counts and policy citations.
  - Overlay disease outbreak data on the timeline.

### Figure 2B and 2C: Cross-Referencing Policy Documents and Abstracts

- **Objective**: These figures analyze how policy documents cite IDM literature, both from the same country and from different countries.
- **Key Steps**:
  - Standardize country codes across datasets using ISO3 codes.
  - Cross-reference policy documents and abstracts by country.
  - Generate world maps showing the distribution of citations, both within and across countries.
  - Calculate and visualize the deviance between self-citing and foreign citations.

### Regression Analysis: Identifying Changepoints in Publication and Citation Trends

- **Objective**: This analysis identifies changepoints in the trends of publication and citation counts over time using piecewise Poisson regression.
- **Key Steps**:
  - Fit piecewise Poisson regression models to publication and citation data.
  - Identify changepoints using the `segmented` package.
  - Compare models using AIC (Akaike Information Criterion) to determine the best fit.
  - Visualize the fitted models alongside the raw data.

## Dependencies

The code relies on the following R packages:

- `ggplot2` for data visualization.
- `dplyr` and `tidyr` for data manipulation.
- `janitor` for cleaning column names.
- `segmented` for piecewise regression analysis.
- `sf` and `rnaturalearth` for geographic mapping.
- `ggrepel` for better label placement in plots.

To install these packages, run the following code in R:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "janitor", "segmented", "sf", "rnaturalearth", "ggrepel"))
```

## How to Run the Code

1. **Clone the Repository**: Clone this repository to your local machine.
2. **Download the Data**: Ensure that the data files are placed in the `data/` folder.
3. **Open R Markdown Files**: Open the `.Rmd` files in RStudio.
4. **Run the Code**: Knit the R Markdown files to generate the plots and analyses. The output will be saved in the `plots/` and `results/` folders.

## Contact Information

For any questions or issues related to the code, please contact:

- **Name**: Paula Christen
- **Email**: paula.christen16@imperial.ac.uk 
- **Institution**: Imperial College London / University of Nairobi


---

Thank you for your interest in our research! We hope this code is useful for your own work or for reproducing the results of our manuscript.