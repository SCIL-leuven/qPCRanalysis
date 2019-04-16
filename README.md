# qPCRanalysis
Functions for easy analysis of qPCR data

## WHY USE THIS PACKAGE
* Easy to use
* Calculate delta CT and delta delta CT
* Integrated statistics

## Requirements

Install these packages in Rstudio:

  `install.packages(c("devtools", "tidyverse", "lazyeval"))`

## Install

Installation of this package

  `devtools::install_github("SCIL-leuven/qPCRanalysis")`

## Load packages

```
library(readxl)
library(tidyr)
library(dplyr)
<<<<<<< HEAD
library(ggplot2)
library(qPCRanalysis)
library(ggpubr)
library(qPCRanalysis)
library(lazyeval)
```

## Load data
We load the data with the `read_excel()` function of the readxl package. The data should consist of at least 4 columns:

* __Sample__ : invidividual names per sample
* __Variable__ : the biological variability in your data set like a disease, different timepoints, perturbation ...
* __Target__ : the gene names of your primers
* __CT__ : the raw CT values

```
qpcr <- read_excel("qPCR_data.xlsx", col_names = TRUE)
qpcr$CT <- as.numeric(qpcr$CT)
head(qpcr)
```
```
## # A tibble: 6 x 4
##   Sample Genotype Target    CT
##   <chr>  <chr>    <chr>  <dbl>
## 1 WT1    Healthy  Tnf     33.3
## 2 WT1    Healthy  Il1b    34.0
## 3 WT1    Healthy  Ccl3    NA  
## 4 WT1    Healthy  Ccl4    NA  
## 5 WT1    Healthy  Ccl5    32.4
## 6 WT1    Healthy  Ccl6    29.7
```

## Calculate Delta CT
To calculate delta CT we use the `calculate_DCT()` function. This function requires four arguments:

* __df__ : dataframe structured like the proposed data file
* __hkg__ : name of housekeeping gene or genes that you want to use to normalize against
* __sample_col__ : name of the sample column
* __gene_col__ : name of the gene column

It will pass a dataframe with three added columns: 

* __CT_hkg__ : average CT value of housekeeping genes
* __DCT__ : Delta CT values
* __RE__ : relative expression to hkg

```
qpcr <- calculate_DCT(df = qpcr, hkg = c("Rab35", "Rpl13a", "PSma3"), sample_col = "Sample", gene_col = "Target")
```
```
## # A tibble: 264 x 7
## # Groups:   Sample [12]
##    Sample Genotype Target    CT CT_hkg    DCT       RE
##    <chr>  <chr>    <chr>  <dbl>  <dbl>  <dbl>    <dbl>
##  1 WT1    Healthy  Tnf     33.3   24.9  -8.42  0.00291
##  2 WT1    Healthy  Il1b    34.0   24.9  -9.14  0.00177
##  3 WT1    Healthy  Ccl3    NA     24.9  NA    NA      
##  4 WT1    Healthy  Ccl4    NA     24.9  NA    NA      
##  5 WT1    Healthy  Ccl5    32.4   24.9  -7.58  0.00522
##  6 WT1    Healthy  Ccl6    29.7   24.9  -4.85  0.0346 
##  7 WT1    Healthy  Ccl7    34.6   24.9  -9.76  0.00116
##  8 WT1    Healthy  Ccl8    30.5   24.9  -5.61  0.0205 
##  9 WT1    Healthy  Ccl9    29.7   24.9  -4.83  0.0351 
## 10 WT1    Healthy  Ccl12   32.3   24.9  -7.49  0.00555
## # ... with 254 more rows
```

## Calculate Delta Delta CT
To calculate Delta Delta CT use the `calculate_DDCT()` function. This function  can only be run after the `calculate_DCT()` function is used and requires five argeuments:

* __df__: dataframe structured like the proposed data file, has to contain DCT and RE column
* __gene_col__ : name of the gene column
* __sample_col__ : name of the sample column
* __var_col__ : column name of variables to normalize your control against
* __control__: name of variable to use as control

It will pass a dataframe with seven added columns

* __DCTavg__ : average Delta CT value
* __DCTsem__ : s.e.m. of Delta CT
* __DCTsemperc__ : percentage of s.e.m.
* __DDCTavg__ : average Delta Delta CT value
* __DDCTsem__ : standard error to the mean of Delta Delta CT
* __DDCTmin__ : minimum sem value
* __DDCTmax__ : maximum sem value

```
qpcr <- calculate_DDCT(df = qpcr, 
                        gene_col = "Target", 
                         sample_col = "Sample", 
                         var_col = "Genotype", 
                         control = "Healthy")
head(qpcr)
```
```
## # A tibble: 6 x 9
##   Target Genotype  DCTavg DCTsem DCTsemperc DDCTavg DDCTsem DDCTmin DDCTmax
##   <chr>  <chr>      <dbl>  <dbl>      <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
## 1 Ccl12  Dystrop~   -6.03  0.319      -5.30    15.6   2.04   13.6     17.7 
## 2 Ccl12  Healthy   -10.00  1.19      -11.9      1     0.169   0.831    1.17
## 3 Ccl12  NA        NaN    NA          NA      NaN    NA     NaN       NA   
## 4 Ccl17  Dystrop~   -8.36  0.416      -4.97    14.7   1.37   13.3     16.0 
## 5 Ccl17  Healthy   -12.2   0.963      -7.87     1     0.111   0.889    1.11
## 6 Ccl17  NA        NaN    NA          NA      NaN    NA     NaN       NA
```

## Examples

Example analysis available in vignette folder

