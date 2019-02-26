# qPCRanalysis
Functions for easy analysis of qPCR data

## WHY USE THIS PACKAGE
* Easy to use
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

### Load packages

```
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(qPCRanalysis)
library(ggpubr)
library(lazyeval)
```

### Load data

```
qpcr <- read_excel("data/fig1_skeletal_muscle/2016-02-09 174725_Chemokines_WT_Sgca_skMuscle_Jordi-ViiA7-export.xlsx", sheet = 3, skip = 35, col_names = TRUE)
qpcr$CT <- as.numeric(qpcr$CT)
head(qpcr)
```
```
## # A tibble: 6 x 35
##    Well `Well Position` Omit  `Sample Name` `Target Name` Task  Reporter
##   <dbl> <chr>           <lgl> <chr>         <chr>         <chr> <chr>   
## 1     1 A1              FALSE WT1           Tnf           UNKN~ SYBR    
## 2     2 A2              FALSE WT1           Il1b          UNKN~ SYBR    
## 3     3 A3              FALSE WT1           Ccl3          UNKN~ SYBR    
## 4     4 A4              FALSE WT1           Ccl4          UNKN~ SYBR    
## 5     5 A5              FALSE WT1           Ccl5          UNKN~ SYBR    
## 6     6 A6              FALSE WT1           Ccl6          UNKN~ SYBR    
## # ... with 28 more variables: Quencher <chr>, CT <dbl>, `Ct Mean` <dbl>,
## #   `Ct SD` <lgl>, Quantity <lgl>, `Quantity Mean` <lgl>, `Quantity
## #   SD` <lgl>, `Automatic Ct Threshold` <lgl>, `Ct Threshold` <dbl>,
## #   `Automatic Baseline` <lgl>, `Baseline Start` <dbl>, `Baseline
## #   End` <dbl>, Comments <lgl>, `Y-Intercept` <lgl>, `R(superscript
## #   2)` <lgl>, Slope <lgl>, Tm1 <dbl>, Tm2 <dbl>, Tm3 <lgl>,
## #   Custom1 <lgl>, Custom2 <lgl>, Custom3 <lgl>, Custom4 <lgl>,
## #   Custom5 <lgl>, Custom6 <lgl>, EXPFAIL <chr>, THOLDFAIL <chr>,
## #   MTP <chr>
```



### Calculate Delta CT

To calculate delta CT we use the `calculate_DCT()` function. This function requires four arguments:
- **df** : dataframe structured like the proposed data file
- **hkg** : name of housekeeping gene or genes that you want to use to normalize against
- **sample_col** : name of the sample column
- **gene_col** : name of the gene column

It will pass a dataframe with two added columns: 
- **DCT** : Delta CT values
- **RE** : relative expression to hkg

```
qpcr <- calculate_DCT(df = qpcr, hkg = c("Rab35", "Rpl13a", "Psma3"), sample_col = "Sample", gene_col = "Target")
```
```
## # A tibble: 230 x 7
## # Groups:   Sample Name [10]
##    `Sample Name` `Target Name`    CT genotype CT_hkg Delta_ct rel_expr
##    <chr>         <chr>         <dbl> <chr>     <dbl>    <dbl>    <dbl>
##  1 WT1           Tnf            33.3 Healthy    23.0   -10.3     1257.
##  2 WT1           Il1b           34.0 Healthy    23.0   -11.0     2065.
##  3 WT1           Ccl3           40   Healthy    23.0   -17.0   132718.
##  4 WT1           Ccl4           40   Healthy    23.0   -17.0   132718.
##  5 WT1           Ccl5           32.4 Healthy    23.0    -9.45     702.
##  6 WT1           Ccl6           29.7 Healthy    23.0    -6.72     106.
##  7 WT1           Ccl7           34.6 Healthy    23.0   -11.6     3169.
##  8 WT1           Ccl8           30.5 Healthy    23.0    -7.48     179.
##  9 WT1           Ccl9           29.7 Healthy    23.0    -6.70     104.
## 10 WT1           Ccl12          32.3 Healthy    23.0    -9.36     659.
## # ... with 220 more rows
```

### Calculate Delta Delta CT

To calculate Delta Delta CT use the `calculate_DDCT()` function. This function  can only be run after the `calculate_DCT()` function is used and requires five argeuments:
- **df**: dataframe structured like the proposed data file, has to contain DCT and RE column
- **gene_col** : name of the gene column
- **sample_col** : name of the sample column
- **var_col** : column name of variables to normalize your control against
- **control**: name of variable to use as control

It will pass a dataframe with four added columns
- **DDCTavg** : average Delta Delta CT values
- **DDCTsem** : standard error to the mean of Delta Delta CT
- **DDCTmin** : minimum sem value
- **DDCTmax** : maximum sem value

```
qpcr <- calculate_DDCT(df = qpcr, gene_col = "Target", sample_col = "Sample", var = "Celltype", control = "Fibroblast")
```

## Examples

Example vignettes will be available soon

