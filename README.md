qPCRanalysis
================
Jordi Camps
May 10 2019

## Installation

Install following packages

``` r
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(ggpubr)) install.packages("ggpubr") 
if(!require(lazyeval)) install.packages("lazyeval")
if(!require(WriteXLS)) install.packages("WriteXLS")
if(!require(devtools)) install.packages("devtools")
if(!require(qPCRanalysis)) devtools::install_github("SCIL-leuven/qPCRanalysis")
```

## Load packages

``` r
library(tidyverse)
library(readxl)
library(qPCRanalysis)
library(ggpubr)
library(lazyeval)
```

## Load data

We load the data with the **read\_excel()** function of the readxl
package. The data should consist of at least 4 columns:

  - **Sample** : invidividual names per sample
  - **Gene** : the gene names of your primers
  - **CT** : the raw CT values
  - **A grouping variable** : the biological variability in your data
    set like a disease, different timepoints, perturbation
…

<!-- end list -->

``` r
qpcr <- read_excel("vignette/qPCR_data.xlsx", col_names = TRUE, sheet = 3, skip = 35)
head(qpcr)
```

    ## # A tibble: 6 x 35
    ##    Well `Well Position` Omit  `Sample Name` `Target Name` Task  Reporter
    ##   <dbl> <chr>           <lgl> <chr>         <chr>         <chr> <chr>   
    ## 1     1 A1              FALSE WT_1          Tnf           UNKN~ SYBR    
    ## 2     2 A2              FALSE WT_1          Il1b          UNKN~ SYBR    
    ## 3     3 A3              FALSE WT_1          Ccl3          UNKN~ SYBR    
    ## 4     4 A4              FALSE WT_1          Ccl4          UNKN~ SYBR    
    ## 5     5 A5              FALSE WT_1          Ccl5          UNKN~ SYBR    
    ## 6     6 A6              FALSE WT_1          Ccl6          UNKN~ SYBR    
    ## # ... with 28 more variables: Quencher <chr>, CT <chr>, `Ct Mean` <dbl>,
    ## #   `Ct SD` <lgl>, Quantity <lgl>, `Quantity Mean` <lgl>, `Quantity
    ## #   SD` <lgl>, `Automatic Ct Threshold` <lgl>, `Ct Threshold` <dbl>,
    ## #   `Automatic Baseline` <lgl>, `Baseline Start` <dbl>, `Baseline
    ## #   End` <dbl>, Comments <lgl>, `Y-Intercept` <lgl>, `R(superscript
    ## #   2)` <lgl>, Slope <lgl>, Tm1 <dbl>, Tm2 <dbl>, Tm3 <lgl>,
    ## #   Custom1 <lgl>, Custom2 <lgl>, Custom3 <lgl>, Custom4 <lgl>,
    ## #   Custom5 <lgl>, Custom6 <lgl>, EXPFAIL <chr>, THOLDFAIL <chr>,
    ## #   MTP <chr>

## Prepare data

### Select columns

Now we need to clean the data frame. We will select only the columns
that we need and change their names because R does not like spaces
inside column names.

``` r
qpcr <- select(qpcr, c("Sample Name", "Target Name", "CT"))
colnames(qpcr) <- c("Sample", "Gene", "CT")
head(qpcr)
```

    ## # A tibble: 6 x 3
    ##   Sample Gene  CT                
    ##   <chr>  <chr> <chr>             
    ## 1 WT_1   Tnf   33.277999877929688
    ## 2 WT_1   Il1b  33.993999481201172
    ## 3 WT_1   Ccl3  Undetermined      
    ## 4 WT_1   Ccl4  Undetermined      
    ## 5 WT_1   Ccl5  32.437000274658203
    ## 6 WT_1   Ccl6  29.705999374389648

### Switch CT column to numeric object

When we check the structure of the data frame we see that the CT column
is not numeric, we will change this. In addition, all undetermined
values are switched to NA. If desired, these values can be changed to
40.

``` r
#Change CT column to numeric values
qpcr$CT <- as.numeric(qpcr$CT)
```

    ## Warning: NAs introduced by coercion

``` r
#Change NA values in CT column to 40
qpcr[is.na(qpcr$CT), "CT"] <- 40
head(qpcr)
```

    ## # A tibble: 6 x 3
    ##   Sample Gene     CT
    ##   <chr>  <chr> <dbl>
    ## 1 WT_1   Tnf    33.3
    ## 2 WT_1   Il1b   34.0
    ## 3 WT_1   Ccl3   40  
    ## 4 WT_1   Ccl4   40  
    ## 5 WT_1   Ccl5   32.4
    ## 6 WT_1   Ccl6   29.7

### Create a grouping variable

In this example we want to see the difference between different
genotypes. We will split the Sample Name column in two to create a
**sample** column and a **replicate** column.

``` r
qpcr <- qpcr %>%
  #Copy Sample column
  mutate(temp = Sample) %>%
  #Split Sample column to create genotype column
  separate(col = temp, into = c("Genotype", "temp"), sep = "_") %>%
  #remove temp column
  select(-temp)
```

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 24 rows
    ## [265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279,
    ## 280, 281, 282, 283, 284, ...].

``` r
head(qpcr)
```

    ## # A tibble: 6 x 4
    ##   Sample Gene     CT Genotype
    ##   <chr>  <chr> <dbl> <chr>   
    ## 1 WT_1   Tnf    33.3 WT      
    ## 2 WT_1   Il1b   34.0 WT      
    ## 3 WT_1   Ccl3   40   WT      
    ## 4 WT_1   Ccl4   40   WT      
    ## 5 WT_1   Ccl5   32.4 WT      
    ## 6 WT_1   Ccl6   29.7 WT

### Remove blanc rows

As a last step in cleaning the data we will still remove the blanc rows.

``` r
qpcr <- filter(qpcr, Sample != "Blanc")
```

## Calculate Delta CT

To calculate delta CT we use the `calculate_DCT()` function. This
function requires four arguments:

  - **df** : dataframe structured like the proposed data file
  - **hkg** : name of housekeeping gene or genes that you want to use to
    normalize against
  - **sample\_col** : name of the sample column
  - **gene\_col** : name of the gene column

It will pass a dataframe with three added columns:

  - **hkg** : names of housekeeping genes used for normalizing
  - **CT\_hkg** : average CT value of housekeeping genes
  - **DCT** : Delta CT values
  - **RE** : relative expression to hkg

<!-- end list -->

``` r
qpcr <- calculate_DCT(df = qpcr, 
                      hkg = c("Rab35", "Rpl13a", "Psma3"), 
                      sample_col = "Sample", 
                      gene_col = "Gene")
```

    ## Joining, by = "Sample"

    ## # A tibble: 231 x 8
    ## # Groups:   Sample [11]
    ##    Sample Gene     CT Genotype hkg                CT_hkg    DCT        RE
    ##    <chr>  <chr> <dbl> <chr>    <chr>               <dbl>  <dbl>     <dbl>
    ##  1 WT_1   Tnf    33.3 WT       Rab35_Rpl13a_Psma3   24.2  -9.06 0.00187  
    ##  2 WT_1   Il1b   34.0 WT       Rab35_Rpl13a_Psma3   24.2  -9.78 0.00114  
    ##  3 WT_1   Ccl3   40   WT       Rab35_Rpl13a_Psma3   24.2 -15.8  0.0000177
    ##  4 WT_1   Ccl4   40   WT       Rab35_Rpl13a_Psma3   24.2 -15.8  0.0000177
    ##  5 WT_1   Ccl5   32.4 WT       Rab35_Rpl13a_Psma3   24.2  -8.22 0.00335  
    ##  6 WT_1   Ccl6   29.7 WT       Rab35_Rpl13a_Psma3   24.2  -5.49 0.0222   
    ##  7 WT_1   Ccl7   34.6 WT       Rab35_Rpl13a_Psma3   24.2 -10.4  0.000741 
    ##  8 WT_1   Ccl8   30.5 WT       Rab35_Rpl13a_Psma3   24.2  -6.25 0.0131   
    ##  9 WT_1   Ccl9   29.7 WT       Rab35_Rpl13a_Psma3   24.2  -5.47 0.0225   
    ## 10 WT_1   Ccl12  32.3 WT       Rab35_Rpl13a_Psma3   24.2  -8.13 0.00356  
    ## # ... with 221 more rows

### Statistics

Perform statistical test between two groups and add information to qpcr
data frame. We will do this with the ggpubr package. We want to perform
a t-test between healthy and dystrophic samples. Therefore, we put *DCT*
as our numeric variable and *Genotype* as our grouping variable. We
Specify the method as *t.test*, specify that the data in unpaired and
group by
*Gene*.

``` r
stat <- compare_means(DCT ~ Genotype, qpcr, method = "t.test", paired = FALSE, group.by = "Gene")
```

    ## Adding missing grouping variables: `Sample`

``` r
head(stat)
```

    ## # A tibble: 6 x 9
    ##   Gene  .y.   group1 group2     p p.adj p.format p.signif method
    ##   <chr> <chr> <chr>  <chr>  <dbl> <dbl> <chr>    <chr>    <chr> 
    ## 1 Tnf   DCT   WT     Sgca   0.566     1 0.57     ns       T-test
    ## 2 Il1b  DCT   WT     Sgca   0.630     1 0.63     ns       T-test
    ## 3 Ccl3  DCT   WT     Sgca   0.214     1 0.21     ns       T-test
    ## 4 Ccl4  DCT   WT     Sgca   0.268     1 0.27     ns       T-test
    ## 5 Ccl5  DCT   WT     Sgca   0.772     1 0.77     ns       T-test
    ## 6 Ccl6  DCT   WT     Sgca   0.541     1 0.54     ns       T-test

### Plot

``` r
ggplot(qpcr, aes(x = Genotype, y = DCT, col = Genotype)) +
  geom_boxplot() +
  facet_wrap(~Gene, scales = "free_y") +
  stat_compare_means(method = "t.test", label = "p.signif", label.x = 1.5)
```

![](README_figs/README-unnamed-chunk-11-1.png)<!-- -->

## Calculate Delta Delta CT

If you have unpaired data, the only way to calculate a delta delta CT or
fold change is to take an average of your sample and compate the fold
change to another group that you are referring to, plus also correctly
propagate your error.

To calculate Delta Delta CT use the `calculate_DDCT()` function. This
function can only be run after the `calculate_DCT()` function is used
and requires five argeuments:

  - **df**: dataframe structured like the proposed data file, has to
    contain DCT and RE column
  - **gene\_col** : name of the gene column
  - **sample\_col** : name of the sample column
  - **var\_col** : column name of variables to normalize your control
    against
  - **control**: name of variable to use as control

It will pass a dataframe with seven added columns

  - **DDCTavg** : average Delta Delta CT value
  - **DDCTsem** : standard error to the mean of Delta Delta CT
  - **DDCTmin** : minimum sem value
  - **DDCTmax** : maximum sem value

<!-- end list -->

``` r
ddct <- calculate_DDCT(df = qpcr, 
                       gene_col = "Gene", 
                       sample_col = "Sample", 
                       var_col = "Genotype", 
                       control = "WT")
head(ddct)
```

    ## # A tibble: 6 x 6
    ##   Gene  Genotype DDCTavg DDCTsem DDCTmin DDCTmax
    ##   <chr> <chr>      <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 Ccl12 Sgca       12.9    3.48    9.47    16.4 
    ## 2 Ccl12 WT          1      0.371   0.629    1.37
    ## 3 Ccl17 Sgca        3.36   0.729   2.63     4.09
    ## 4 Ccl17 WT          1      0.298   0.702    1.30
    ## 5 Ccl22 Sgca        2.84   0.652   2.19     3.49
    ## 6 Ccl22 WT          1      0.302   0.698    1.30

### Plot

``` r
ggplot(ddct, aes(x = Genotype, y = DDCTavg, fill = Genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = DDCTmin, ymax = DDCTmax), width = 0.1, data = ddct) +
  facet_wrap(~Gene, scales = "free_y")
```

![](README_figs/README-unnamed-chunk-13-1.png)<!-- -->

## Export

At any given point you can export the data frame to an excel file with
this function. Here we will export the data in an excel with three
sheets. The first sheet will contain your raw values and normalized
delta CT values. The second sheet will contain your statistics. The
third sheet will contain all the values normalized to your control
group.

``` r
library(WriteXLS)
WriteXLS(c("qpcr", "stat", "ddct"), "vignette/qpcr.xlsx")
```
