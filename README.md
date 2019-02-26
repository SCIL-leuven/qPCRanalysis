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


Well | Well Position | Sample Name | Target Name | CT  
-----|---------------|-------------|-------------|----
1    | A1            | Sample 1    | Gene 1    | 21.23  
2    | A2            | Sample 2    | Gene 1    | 24.48  
3    | A3            | Sample 1    | Gene 2    | 18.98   


## Annotation file

To plot groups, you need to create an annotation file (this is the easiest to do in excel but it can also be a .csv or .txt file). Here is an example of how it should look:

Sample Name | Condition | Mouse | ...
----------|-----------|-----------|------
Sample 1  | healthy   | a  | ...
Sample 2  | healthy | b  | ...
Sample 3  | diseased | a | ...

## Workflow

1. Delta CT
2. Delta Delta CT
