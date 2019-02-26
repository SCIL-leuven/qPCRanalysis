# qPCRanalysis
Functions for easy analysis of qPCR data

## WHY USE THIS PACKAGE
* Simplified workflow
* Calculate delta CT and delta delta CT

## Requirements

Install these packages in Rstudio:

  `install.packages(c("devtools", "tidyverse", "lazyeval"))`

## Install

Installation of this package

  `devtools::install_github("SCIL-leuven/qPCRanalysis")`

## Data file

Import all the raw data as an excel file without any selections. We are only interested in the results, that is the third sheet. In the basic template on the Viia7, genes are **Targets** and samples **Samples**. You can find them under the Sample Name and Target Name column.


Sample | Target | CT  
-------|--------|----
Sample 1    | Gene 1    | 21.23  
Sample 2    | Gene 2    | 24.48  
Sample 3    | Gene 3    | 18.98   


## Workflow

### Calculate Delta CT

To calculate delta CT we use the `calculate_DCT()` function. This function requires four arguments:
- df : dataframe structured like the proposed data file
- hkg : name of housekeeping gene or genes that you want to use to normalize against
- sample_col : name of the sample column
- gene_col : name of the gene column

It will pass a dataframe with two added columns: the Delta CT (DCT column) and the relative expression to hkg (RE column)

2. Delta Delta CT
