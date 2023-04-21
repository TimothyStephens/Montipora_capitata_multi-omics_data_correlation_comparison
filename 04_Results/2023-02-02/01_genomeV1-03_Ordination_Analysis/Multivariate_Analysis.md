---
title: "M. capitata genome v1 Ordination Analysis"
author: "Timothy Stephens"
date: "02/03/2023"
output: 
  html_document:
    keep_md: yes
---



# Setup

Setup R env. Load packages and set default image export formats, size and resolution.


```r
knitr::opts_chunk$set(echo = TRUE,
                      fig.height = 8, 
                      fig.width = 12, 
                      dev = c("png", "pdf"),
                      dpi = 1000)
library(lubridate)
```

```
## 
## Attaching package: 'lubridate'
```

```
## The following objects are masked from 'package:base':
## 
##     date, intersect, setdiff, union
```

```r
#library(tidyverse)
library(seacarb)
```

```
## Loading required package: oce
```

```
## Loading required package: gsw
```

```
## Loading required package: SolveSAPHE
```

```r
library(matrixStats)
library(vegan)
```

```
## Loading required package: permute
```

```
## Loading required package: lattice
```

```
## This is vegan 2.6-4
```

```r
library(lme4)
```

```
## Loading required package: Matrix
```

```r
library(ape)
library(emmeans)
library(gridExtra)
library(multcompView)
library(plotrix)
```

```
## 
## Attaching package: 'plotrix'
```

```
## The following object is masked from 'package:oce':
## 
##     rescale
```

```r
library(reshape2)
library(ggpubr)
```

```
## Loading required package: ggplot2
```

```
## 
## Attaching package: 'ggpubr'
```

```
## The following object is masked from 'package:ape':
## 
##     rotate
```

```r
library(sva, warn.conflicts = FALSE)
```

```
## Loading required package: mgcv
```

```
## Loading required package: nlme
```

```
## 
## Attaching package: 'nlme'
```

```
## The following object is masked from 'package:lme4':
## 
##     lmList
```

```
## This is mgcv 1.8-41. For overview type 'help("mgcv-package")'.
```

```
## Loading required package: genefilter
```

```
## 
## Attaching package: 'genefilter'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     rowSds, rowVars
```

```
## Loading required package: BiocParallel
```

```r
library(mixOmics)
```

```
## Loading required package: MASS
```

```
## 
## Attaching package: 'MASS'
```

```
## The following object is masked from 'package:genefilter':
## 
##     area
```

```
## 
## Loaded mixOmics 6.18.1
## Thank you for using mixOmics!
## Tutorials: http://mixomics.org
## Bookdown vignette: https://mixomicsteam.github.io/Bookdown
## Questions, issues: Follow the prompts at http://mixomics.org/contact-us
## Cite us:  citation('mixOmics')
```

```r
library(ggforce)
library(plyr)
```

```
## 
## Attaching package: 'plyr'
```

```
## The following object is masked from 'package:ggpubr':
## 
##     mutate
```

```
## The following object is masked from 'package:matrixStats':
## 
##     count
```

```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:plyr':
## 
##     arrange, count, desc, failwith, id, mutate, rename, summarise,
##     summarize
```

```
## The following object is masked from 'package:MASS':
## 
##     select
```

```
## The following object is masked from 'package:nlme':
## 
##     collapse
```

```
## The following object is masked from 'package:gridExtra':
## 
##     combine
```

```
## The following object is masked from 'package:ape':
## 
##     where
```

```
## The following object is masked from 'package:matrixStats':
## 
##     count
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(plotly)
```

```
## 
## Attaching package: 'plotly'
```

```
## The following objects are masked from 'package:plyr':
## 
##     arrange, mutate, rename, summarise
```

```
## The following object is masked from 'package:MASS':
## 
##     select
```

```
## The following object is masked from 'package:ggplot2':
## 
##     last_plot
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

```
## The following object is masked from 'package:graphics':
## 
##     layout
```

```r
library(gapminder)
library(htmlwidgets)
options(scipen = 999) #Prevent scientific notation
set.seed(54321)

#source("ggbiplot/R/ggbiplot.r")
```



Plot settings.

```r
options(repr.plot.width = 2, repr.plot.height = 3)
label.decimal.rounding <- 1
plot.text.size  <- 12
plot.title.size <- 10

Colors <- c("T1-Amb" = "#00ebf7", "T3-Amb" = "#01ff00", "T5-Amb" = "#006600", "Field" = "#660099", "T1-HiT" = "#fad700", "T3-HiT" = "#ff7c00", "T5-HiT" = "#ff0000")
```



# Functions

Functions for plotting and analyzing data.


```r
#' Run PERMANOVA analysis
#' 
#' @param data expression/accumulation dataframe (rownames=gene_names; colnames=sample_names)
#' @param samplesInfo sample_names metadata (required: "TimePoint", "Treatment", and "Tank")
#' @param out.prefix prefix of output csv file with PERMANOVA results
#' @param permutations number of permutations to perform (default: 999)
#' @return PERMANOVA results from vegan::adonis2
run_PERMANOVA <- function(data, samplesInfo, out.prefix, permutations=999){
  # Transpose df (samples as rownames)
  t.data <- t(data)
  t.data <- as.data.frame(t.data)

  # Identify factors
  #  - Will only return overlap between data and samplesInfo
  #  - Make the "by" column (which was used for merging) rownames.
  t.data.sampleInfo <- merge(samplesInfo, t.data, by = 0) %>%
    tibble::column_to_rownames(var = "Row.names")
  t.data.sampleInfo <- t.data.sampleInfo[row.names(t.data), ]
  
  # Use square root or proportions to minimize influence of most abundant groups
  t.data.sqrt <- sqrt(t.data)
  
  # Create a dissimilarity matrix (vegan)
  t.data.dist <- vegdist(t.data.sqrt, method='bray')
  
  # Run perMANOVA (vegan)
  mod.data <- adonis2(t.data.dist ~  TimePoint * Treatment, 
                      data = t.data.sampleInfo, 
                      permutations = permutations, 
                      strata = t.data.sampleInfo$Tank, 
                      method = 'bray')
  
  # Write PERMANOVA results to file
  write.csv(mod.data, paste(out.prefix,".PERMANOVA_results.csv", sep=''))
  
  # Return PERMANOVA results from function
  return(mod.data)
}
```




```r
#' Plot PCA
#' 
#' @param data expression/accumulation dataframe (rownames=gene_names; colnames=sample_names)
#' @param samplesInfo sample_names metadata (required: "Condition", "Treatment", and "Plug.ID")
#' @param out.prefix prefix of output html file with PCA plot
#' @param plot.title plot title
#' @return PCA plot as a ggplot2 object
plot_PCA <- function(data, samplesInfo, out.prefix, plot.title){
  # Transpose df (samples as rownames)
  t.data <- t(data)
  t.data <- as.data.frame(t.data)
  
  # Run PCA
  pca.out <- prcomp(t.data,
                    center = FALSE,
                    scale. = FALSE)
  pca.out.summary <- summary(pca.out)
  #print(pca.out.summary)
  
  # Create plotting dataframe by pulling out axes and joining to metadata
  PC1.df <- data.frame(pca.out$x[,1])
  names(PC1.df)[1] <- 'PC1'
  
  PC2.df <- data.frame(pca.out$x[,2])
  names(PC2.df)[1] <- 'PC2'
  
  # Merge PCA results with metadata sampleInfo file (used for easy access to variables for plotting)
  MyMerge <- function(x, y){
    return(merge(x, y, by=0, all.x=F, all.y= F) %>% 
             tibble::column_to_rownames(var = "Row.names"))
  }
  pca.data <- Reduce(MyMerge, list(samplesInfo, PC1.df, PC2.df))
  #print(pca.data)
  
  # Calculate the percent of variance explained by first two axes and format for plot labels
  PCA.axis1.var <- pca.out.summary$importance[2,1]
  PCA.axis2.var <- pca.out.summary$importance[2,2]
  PCA.axis1.var.label <- paste("PC1 (",round(PCA.axis1.var*100,label.decimal.rounding),"%)", sep='')
  PCA.axis2.var.label <- paste("PC2 (",round(PCA.axis2.var*100,label.decimal.rounding),"%)", sep='')
  
  # Plot PCA
  pca.plot <- ggplot(pca.data, aes(x = PC1, y = PC2)) +
              geom_point(aes(color = Condition, shape = Treatment, size = 2)) +
              scale_color_manual(values = Colors) +
              geom_mark_ellipse(aes(color = Condition),
                                expand = unit(1,"mm")) +
              labs(title=plot.title,
                   x=PCA.axis1.var.label,
                   y=PCA.axis2.var.label) +
              geom_text(label=pca.data$Plug.ID,
                        check_overlap=F) +
              theme(text       = element_text(family="Helvetica", size = plot.text.size),
                    plot.title = element_text(family="Helvetica", size = plot.title.size, hjust = 0.5),
                    aspect.ratio=1)
  
  # Write PCA plot as interactive HTML file
  saveWidget(suppressWarnings(ggplotly(pca.plot)), file=paste(out.prefix,".PCA.html", sep=''))
  
  # Return PCA plot from function
  return(pca.plot)
}
```




```r
#' Plot PCoA
#' 
#' @param data expression/accumulation dataframe (rownames=gene_names; colnames=sample_names)
#' @param samplesInfo sample_names metadata (required: "Condition", "Treatment", and "Plug.ID")
#' @param out.prefix prefix of output html file with PCoA plot
#' @param plot.title plot title
#' @return PCoA plot as a ggplot2 object
plot_PCoA <- function(data, samplesInfo, out.prefix, plot.title){
  # Transpose df (samples as rownames)
  t.data <- t(data)
  t.data <- as.data.frame(t.data)
  
  # Identify factors
  #  - Will only return overlap between data and samplesInfo
  #  - Make the "by" column (which was used for merging) rownames.
  t.data.sampleInfo <- merge(samplesInfo, t.data, by = 0) %>%
    tibble::column_to_rownames(var = "Row.names")
  
  # Use square root or proportions to minimize influence of most abundant groups
  t.data.sqrt <- sqrt(t.data)
  
  # Create a dissimilarity matrix (vegan)
  t.data.dist <- vegdist(t.data.sqrt, method='bray')
  
  # Run PCoA (vegan)
  pcoa <- wcmdscale(t.data.dist, k = 2, eig = TRUE, add = "cailliez")
  
  # Create plotting dataframe by pulling out axes and joining to metadata
  axis1 <- pcoa$points[,1]
  axis1.df <- data.frame(axis1)
  axis2 <- pcoa$points[,2]
  axis2.df <- data.frame(axis2)
  
  # Merge PCA results with metadata sampleInfo file (used for easy access to variables for plotting)
  MyMerge <- function(x, y){
    return(merge(x, y, by=0, all.x=F, all.y= F) %>% 
             tibble::column_to_rownames(var = "Row.names"))
  }
  pcoa.data <- Reduce(MyMerge, list(samplesInfo, axis1.df, axis2.df))
  #print(pcoa.data)
  
  # Calculate the percent of variance explained by first two axes.
  pcoa.axis1.var <- sum((as.vector(pcoa$eig)/sum(pcoa$eig))[1])
  pcoa.axis2.var <- sum((as.vector(pcoa$eig)/sum(pcoa$eig))[2])
  pcoa.axis1.var.label <- paste("PCoA 1 (",round(pcoa.axis1.var*100,label.decimal.rounding),"%)", sep='')
  pcoa.axis2.var.label <- paste("PCoA 2 (",round(pcoa.axis2.var*100,label.decimal.rounding),"%)", sep='')
  
  pcoa.plot <- ggplot(pcoa.data, aes(x = axis1, y = axis2)) +
    geom_point(aes(color = Condition, shape = Treatment, size = 2)) +
    scale_color_manual(values = Colors) +
    geom_mark_ellipse(aes(color = Condition),
                      expand = unit(1,"mm")) +
    labs(title=plot.title,
         x=pcoa.axis1.var.label,
         y=pcoa.axis2.var.label) +
    geom_text(label=pcoa.data$Plug.ID,
              check_overlap=F) +
    theme(text       = element_text(family="Helvetica", size = plot.text.size),
          plot.title = element_text(family="Helvetica", size = plot.title.size, hjust = 0.5),
          aspect.ratio=1)

  # Write PCoA plot as interactive HTML file
  saveWidget(suppressWarnings(ggplotly(pcoa.plot)), file=paste(out.prefix,".PCoA.html", sep=''))
  
  # Return PCoA plot from function
  return(pcoa.plot)
}
```




```r
#' Plot PLS-DA
#' 
#' @param data expression/accumulation dataframe (rownames=gene_names; colnames=sample_names)
#' @param samplesInfo sample_names metadata (required: "Group", "Condition", "Treatment", and "Plug.ID"). Where "Group" is Genotype+TimePoint+Treatment
#' @param out.prefix prefix of output html file with PLS-DA plot
#' @param plot.title plot title
#' @return PLS-DA plot as a ggplot2 object
plot_PLSDA <- function(data, samplesInfo, out.prefix, plot.title){
  # Transpose df (samples as rownames)
  t.data <- t(data)
  t.data <- as.data.frame(t.data)
  
  # Set "group" as vector
  samplesInfo <- samplesInfo %>% filter(row.names(.) %in% colnames(data))
  sample.info <- samplesInfo$Group
  sample.info <- as.factor(sample.info)

  # Determine the number of components to retain (ncomp -> K−1, where K is the number of classes)
  # In this case the components are the treatemnt "Groups" i.e., 
  #      - MC-289_Field MC-289_T1-Amb MC-289_T1-HiT MC-289_T3-Amb MC-289_T3-HiT MC-289_T5-Amb MC-289_T5-HiT
  # Total 7 levels (6 ncomp's)
  sample.class.number <- length(levels(sample.info))
  k <- sample.class.number-1
  #print(paste("Running with k=",k, sep=''))
  
  # Check data for correct dimensions
  #print(summary(sample.info))
  #print(length(sample.info))
  #print(dim(t.data))
  
  # Create PLSDA object
  plsda <- splsda(t.data, sample.info, ncomp = k, scale = FALSE, near.zero.var = TRUE) 
  
  # Pull PLS-DA 1 & 2
  plsda.df <- data.frame(plsda$variates$X[,1:2])
  names(plsda.df)[1] <- 'PLSDA1'
  names(plsda.df)[2] <- 'PLSDA2'
  
  # Merge PCA results with metadata sampleInfo file (used for easy access to variables for plotting)
  MyMerge <- function(x, y){
    return(merge(x, y, by=0, all.x=F, all.y= F) %>% 
             tibble::column_to_rownames(var = "Row.names"))
  }
  plsda.data <- Reduce(MyMerge, list(samplesInfo, plsda.df))
  #print(plsda.data)
  
  # Calculate the percent of variance explained by first two axes.
  plsda$prop_expl_var$X
  plsda.axis1.var <- plsda$prop_expl_var$X[1]
  plsda.axis2.var <- plsda$prop_expl_var$X[2]
  plsda.axis1.var.label <- paste("PLS-DA 1 (",round(plsda.axis1.var*100,label.decimal.rounding),"%)", sep='')
  plsda.axis2.var.label <- paste("PLS-DA 2 (",round(plsda.axis2.var*100,label.decimal.rounding),"%)", sep='')
  
  # Plot initial PLSDA object
  plsda.plot <- ggplot(plsda.data, aes(x = PLSDA1, y = PLSDA2)) +
    geom_point(aes(color = Condition, shape = Treatment, size = 2)) +
    scale_color_manual(values = Colors) +
    geom_mark_ellipse(aes(color = Condition),
                      expand = unit(1,"mm")) +
    labs(title=plot.title,
         x=plsda.axis1.var.label,
         y=plsda.axis2.var.label) +
    geom_text(label=plsda.data$Plug.ID,
              check_overlap=F) +
    theme(text       = element_text(family="Helvetica", size = plot.text.size),
          plot.title = element_text(family="Helvetica", size = plot.title.size, hjust = 0.5),
          aspect.ratio=1)
  
  # Write PLS-DA plot as interactive HTML file
  saveWidget(suppressWarnings(ggplotly(plsda.plot)), file=paste(out.prefix,".PLS-DA.html", sep=''))
  
  # Return PLS-DA plot from function
  return(plsda.plot)
}
```





# Sample metadata

Load sample metadata to use for all xpression/accumulation datasets.

```r
# Make "Sample.ID" row names
samplesInfo <- read.table("samples_info.txt", h = T, row.names = NULL, sep = "\t", check.names = F) %>% 
  tibble::column_to_rownames(var = "Sample.ID")

# Convert to character variables (so they are treated as descrete and not continuious) and factors.
samplesInfo$Treatment <- as.factor(as.character(samplesInfo$Treatment))
samplesInfo$TimePoint <- as.factor(as.character(samplesInfo$TimePoint))
samplesInfo$Tank      <- as.factor(as.character(samplesInfo$Tank))

samplesInfo
```

```
##                            Group Treatment Genotype Condition TimePoint Plug.ID
## MC-289_Field_289-1  MC-289_Field     Field   MC-289     Field         0   289-1
## MC-289_Field_289-2  MC-289_Field     Field   MC-289     Field         0   289-2
## MC-289_Field_289-3  MC-289_Field     Field   MC-289     Field         0   289-3
## MC-289_T1-Amb_1535 MC-289_T1-Amb   Ambient   MC-289    T1-Amb         1    1535
## MC-289_T1-Amb_1603 MC-289_T1-Amb   Ambient   MC-289    T1-Amb         1    1603
## MC-289_T1-Amb_1609 MC-289_T1-Amb   Ambient   MC-289    T1-Amb         1    1609
## MC-289_T1-HiT_2023 MC-289_T1-HiT      High   MC-289    T1-HiT         1    2023
## MC-289_T1-HiT_2750 MC-289_T1-HiT      High   MC-289    T1-HiT         1    2750
## MC-289_T1-HiT_2878 MC-289_T1-HiT      High   MC-289    T1-HiT         1    2878
## MC-289_T3-Amb_1595 MC-289_T3-Amb   Ambient   MC-289    T3-Amb         3    1595
## MC-289_T3-Amb_2530 MC-289_T3-Amb   Ambient   MC-289    T3-Amb         3    2530
## MC-289_T3-Amb_2741 MC-289_T3-Amb   Ambient   MC-289    T3-Amb         3    2741
## MC-289_T3-HiT_2058 MC-289_T3-HiT      High   MC-289    T3-HiT         3    2058
## MC-289_T3-HiT_2183 MC-289_T3-HiT      High   MC-289    T3-HiT         3    2183
## MC-289_T3-HiT_2998 MC-289_T3-HiT      High   MC-289    T3-HiT         3    2998
## MC-289_T5-Amb_1426 MC-289_T5-Amb   Ambient   MC-289    T5-Amb         5    1426
## MC-289_T5-Amb_1721 MC-289_T5-Amb   Ambient   MC-289    T5-Amb         5    1721
## MC-289_T5-Amb_2874 MC-289_T5-Amb   Ambient   MC-289    T5-Amb         5    2874
## MC-289_T5-HiT_1262 MC-289_T5-HiT      High   MC-289    T5-HiT         5    1262
## MC-289_T5-HiT_1341 MC-289_T5-HiT      High   MC-289    T5-HiT         5    1341
## MC-289_T5-HiT_2998 MC-289_T5-HiT      High   MC-289    T5-HiT         5    2998
##                    Tank
## MC-289_Field_289-1    7
## MC-289_Field_289-2    7
## MC-289_Field_289-3    7
## MC-289_T1-Amb_1535    1
## MC-289_T1-Amb_1603    6
## MC-289_T1-Amb_1609    4
## MC-289_T1-HiT_2023    2
## MC-289_T1-HiT_2750    3
## MC-289_T1-HiT_2878    5
## MC-289_T3-Amb_1595    4
## MC-289_T3-Amb_2530    1
## MC-289_T3-Amb_2741    6
## MC-289_T3-HiT_2058    5
## MC-289_T3-HiT_2183    2
## MC-289_T3-HiT_2998    3
## MC-289_T5-Amb_1426    4
## MC-289_T5-Amb_1721    1
## MC-289_T5-Amb_2874    6
## MC-289_T5-HiT_1262    5
## MC-289_T5-HiT_1341    3
## MC-289_T5-HiT_2998    3
```





# Analyze proteomic data


```r
# Load expression data:
#   - Make "Names" column rownames
#   - Make "NA" values "0"
proteomic.data <- read.table("Mcapitata_V1_proteomic_data.tsv", h = T, row.names = NULL, sep = "\t", check.names = F) %>% 
  tibble::column_to_rownames(var = "Name") %>% 
  replace(is.na(.), 0)
dim(proteomic.data)
```

```
## [1] 4036   14
```

```r
# Set file names
out.prefix <- "Mcapitata_V1_proteomic_data"

run_PERMANOVA(proteomic.data, samplesInfo, out.prefix, permutations=how(nperm=190, minperm=190))
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Blocks:  strata 
## Permutation: free
## Number of permutations: 190
## 
## adonis2(formula = t.data.dist ~ TimePoint * Treatment, data = t.data.sampleInfo, permutations = permutations, method = "bray", strata = t.data.sampleInfo$Tank)
##                     Df  SumOfSqs      R2      F Pr(>F)
## TimePoint            3 0.0089362 0.32257 1.9920 0.1047
## Treatment            1 0.0046062 0.16627 3.0803 0.1466
## TimePoint:Treatment  2 0.0036936 0.13333 1.2350 0.1518
## Residual             7 0.0104674 0.37784              
## Total               13 0.0277034 1.00000
```

```r
proteomic.PCA       <- plot_PCA(     proteomic.data, samplesInfo, out.prefix, "Proteomic data (PCA)")
proteomic.PCoA      <- plot_PCoA(    proteomic.data, samplesInfo, out.prefix, "Proteomic data (PCoA)")
proteomic.PLSDA     <- plot_PLSDA(   proteomic.data, samplesInfo, out.prefix, "Proteomic data (sPLS-DA)")
```





# Analyze transcriptomic data


```r
# Load expression data:
#   - Make "Names" column rownames
#   - Make "NA" values "0"
#   - Ignore rows with a combined (summed) TPM across samples <=100
transcriptomic.data <- read.table("Mcapitata_V1_transcriptomic_data.tsv.TPM.txt", h = T, row.names = NULL, sep = "\t", check.names = F) %>% 
  tibble::column_to_rownames(var = "Name") %>% 
  replace(is.na(.), 0) %>% 
  filter(rowSums(.) > 100)
dim(transcriptomic.data)
```

```
## [1] 16046    21
```

```r
# Set file names
out.prefix <- "Mcapitata_V1_transcriptomic_data"

run_PERMANOVA(transcriptomic.data, samplesInfo, out.prefix)
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Blocks:  strata 
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = t.data.dist ~ TimePoint * Treatment, data = t.data.sampleInfo, permutations = permutations, method = "bray", strata = t.data.sampleInfo$Tank)
##                     Df SumOfSqs      R2      F Pr(>F)  
## TimePoint            3 0.040207 0.28861 2.8343  0.057 .
## Treatment            1 0.020644 0.14819 4.3658  0.104  
## TimePoint:Treatment  2 0.012260 0.08800 1.2963  0.256  
## Residual            14 0.066200 0.47520                
## Total               20 0.139311 1.00000                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
transcriptomic.PCA       <- plot_PCA(     transcriptomic.data, samplesInfo, out.prefix, "Transcriptomic data (PCA)")
transcriptomic.PCoA      <- plot_PCoA(    transcriptomic.data, samplesInfo, out.prefix, "Transcriptomic data (PCoA)")
transcriptomic.PLSDA     <- plot_PLSDA(   transcriptomic.data, samplesInfo, out.prefix, "Transcriptomic data (sPLS-DA)")
```

```
## Warning in Check.entry.pls(X, Y, ncomp, keepX, keepY, mode = mode, scale = scale, : Zero- or near-zero variance predictors.
##  Reset
##             predictors matrix to not near-zero variance predictors.
## 
##             See $nzv for problematic predictors.
```





# Analyze transcriptomic data that overlaps with proteomic data


```r
# Load expression data:
#   - Make "Names" column rownames
#   - Make "NA" values "0"
#   - Ignore rows with a combined (summed) TPM across samples <=100
#   - **Keep only rows in proteomic data**
trans_prot_overlap.data <- transcriptomic.data %>%
  filter(row.names(.) %in% rownames(proteomic.data))
dim(trans_prot_overlap.data)
```

```
## [1] 3638   21
```

```r
# Set file names
out.prefix <- "Mcapitata_V1_transcriptomic-proteomic_overlap_data"

run_PERMANOVA(trans_prot_overlap.data, samplesInfo, out.prefix)
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Blocks:  strata 
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = t.data.dist ~ TimePoint * Treatment, data = t.data.sampleInfo, permutations = permutations, method = "bray", strata = t.data.sampleInfo$Tank)
##                     Df SumOfSqs      R2      F Pr(>F)  
## TimePoint            3 0.031890 0.27902 2.7330  0.060 .
## Treatment            1 0.017760 0.15539 4.5661  0.112  
## TimePoint:Treatment  2 0.010190 0.08916 1.3099  0.259  
## Residual            14 0.054452 0.47643                
## Total               20 0.114292 1.00000                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
trans_prot_overlap.PCA       <- plot_PCA(     trans_prot_overlap.data, samplesInfo, out.prefix, "Transcripts with proteomic evidence (PCA)")
trans_prot_overlap.PCoA      <- plot_PCoA(    trans_prot_overlap.data, samplesInfo, out.prefix, "Transcripts with proteomic evidence (PCoA)")
trans_prot_overlap.PLSDA     <- plot_PLSDA(   trans_prot_overlap.data, samplesInfo, out.prefix, "Transcripts with proteomic evidence (sPLS-DA)")
```





# Combined plots


```r
ggarrange(proteomic.PCoA, transcriptomic.PCoA, trans_prot_overlap.PCoA,
          ncol = 2, nrow = 2,
          labels = "AUTO",
          font.label = list(size = 24, color = "black", face = "bold", family = "Helvetica"),
          hjust = -3,
          common.legend = TRUE,
          legend = "right")
```

```
## Warning: Using the `size` aesthetic in this geom was deprecated in ggplot2 3.4.0.
## ℹ Please use `linewidth` in the `default_aes` field and elsewhere instead.
```

![](Multivariate_Analysis_files/figure-html/combined_plot_main-1.png)<!-- -->



```r
ggarrange(proteomic.PLSDA, transcriptomic.PLSDA, trans_prot_overlap.PLSDA,
          proteomic.PCA,   transcriptomic.PCA,   trans_prot_overlap.PCA,
          ncol = 3, nrow = 2,
          labels = "AUTO",
          font.label = list(size = 24, color = "black", face = "bold", family = "Helvetica"),
          hjust = -0.2,
          common.legend = TRUE,
          legend = "right")
```

![](Multivariate_Analysis_files/figure-html/combined_plot_supp-1.png)<!-- -->





# Session Info


```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-apple-darwin18.7.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /usr/local/Cellar/openblas/0.3.18/lib/libopenblasp-r0.3.18.dylib
## LAPACK: /usr/local/Cellar/r/4.1.2/lib/R/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] htmlwidgets_1.6.1   gapminder_0.3.0     plotly_4.10.1      
##  [4] dplyr_1.1.0         plyr_1.8.8          ggforce_0.4.1      
##  [7] mixOmics_6.18.1     MASS_7.3-58.2       sva_3.42.0         
## [10] BiocParallel_1.28.3 genefilter_1.76.0   mgcv_1.8-41        
## [13] nlme_3.1-162        ggpubr_0.6.0        ggplot2_3.4.1      
## [16] reshape2_1.4.4      plotrix_3.8-2       multcompView_0.1-8 
## [19] gridExtra_2.3       emmeans_1.8.4-1     ape_5.7            
## [22] lme4_1.1-31         Matrix_1.5-3        vegan_2.6-4        
## [25] lattice_0.20-45     permute_0.9-7       matrixStats_0.63.0 
## [28] seacarb_3.3.1       SolveSAPHE_2.1.0    oce_1.7-10         
## [31] gsw_1.1-1           lubridate_1.9.2    
## 
## loaded via a namespace (and not attached):
##   [1] backports_1.4.1        igraph_1.4.1           lazyeval_0.2.2        
##   [4] splines_4.1.2          crosstalk_1.2.0        GenomeInfoDb_1.30.1   
##   [7] digest_0.6.31          htmltools_0.5.4        fansi_1.0.4           
##  [10] magrittr_2.0.3         memoise_2.0.1          cluster_2.1.4         
##  [13] limma_3.50.3           Biostrings_2.62.0      annotate_1.72.0       
##  [16] rARPACK_0.11-0         timechange_0.2.0       colorspace_2.1-0      
##  [19] blob_1.2.3             ggrepel_0.9.3          xfun_0.37             
##  [22] crayon_1.5.2           RCurl_1.98-1.10        jsonlite_1.8.4        
##  [25] survival_3.5-3         glue_1.6.2             polyclip_1.10-4       
##  [28] gtable_0.3.1           zlibbioc_1.40.0        XVector_0.34.0        
##  [31] car_3.1-1              BiocGenerics_0.40.0    abind_1.4-5           
##  [34] scales_1.2.1           mvtnorm_1.1-3          DBI_1.1.3             
##  [37] edgeR_3.36.0           rstatix_0.7.2          Rcpp_1.0.10           
##  [40] viridisLite_0.4.1      xtable_1.8-4           bit_4.0.5             
##  [43] stats4_4.1.2           httr_1.4.5             RColorBrewer_1.1-3    
##  [46] ellipsis_0.3.2         pkgconfig_2.0.3        XML_3.99-0.13         
##  [49] farver_2.1.1           sass_0.4.5             locfit_1.5-9.7        
##  [52] utf8_1.2.3             tidyselect_1.2.0       labeling_0.4.2        
##  [55] rlang_1.0.6            AnnotationDbi_1.56.2   munsell_0.5.0         
##  [58] tools_4.1.2            cachem_1.0.7           cli_3.6.0             
##  [61] generics_0.1.3         RSQLite_2.3.0          broom_1.0.3           
##  [64] evaluate_0.20          stringr_1.5.0          fastmap_1.1.1         
##  [67] yaml_2.3.7             knitr_1.42             bit64_4.0.5           
##  [70] purrr_1.0.1            KEGGREST_1.34.0        compiler_4.1.2        
##  [73] rstudioapi_0.14        png_0.1-8              ggsignif_0.6.4        
##  [76] tibble_3.1.8           tweenr_2.0.2           bslib_0.4.2           
##  [79] stringi_1.7.12         highr_0.10             RSpectra_0.16-1       
##  [82] nloptr_2.0.3           vctrs_0.5.2            pillar_1.8.1          
##  [85] lifecycle_1.0.3        jquerylib_0.1.4        estimability_1.4.1    
##  [88] data.table_1.14.8      cowplot_1.1.1          bitops_1.0-7          
##  [91] corpcor_1.6.10         R6_2.5.1               IRanges_2.28.0        
##  [94] boot_1.3-28.1          withr_2.5.0            S4Vectors_0.32.4      
##  [97] GenomeInfoDbData_1.2.7 parallel_4.1.2         grid_4.1.2            
## [100] tidyr_1.3.0            minqa_1.2.5            rmarkdown_2.20        
## [103] carData_3.0-5          Biobase_2.54.0         ellipse_0.4.3
```
