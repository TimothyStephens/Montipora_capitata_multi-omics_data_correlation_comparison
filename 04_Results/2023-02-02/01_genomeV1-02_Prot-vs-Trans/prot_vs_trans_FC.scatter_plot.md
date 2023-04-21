---
title: "Compare FC prot-vs-trans"
author: "Timothy Stephens"
date: "27/02/2023"
output: 
  html_document:
    keep_md: yes
---



# Setup

Setup R env. Load packages and set default image export formats, size and resolution.


```r
knitr::opts_chunk$set(echo = TRUE,
                      fig.height = 8, 
                      fig.width = 8, 
                      dev = c("png", "pdf"),
                      dpi = 1000)
library(tibble)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
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
library(reshape2)
library(ggplot2)
library(ggpmisc)
```

```
## Loading required package: ggpp
```

```
## 
## Attaching package: 'ggpp'
```

```
## The following object is masked from 'package:ggplot2':
## 
##     annotate
```

```r
library(ggtext)
library(RColorBrewer)
library(ggpubr)
library(plotly)
```

```
## 
## Attaching package: 'plotly'
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
library(htmlwidgets)
options(scipen = 999) #Prevent scientific notation
```



```r
xlim.min <- -10
xlim.max <-  10
ylim.min <- -10
ylim.max <-  10

point.color <- c('DEP&DEG'='#f29790','DEG'='#9cc13f','DEP'='#5acfd2','none'='#d49cff')

scatter_plot <- function(df, title="", tag="", text.scale=1) {
  x.name <- paste(title,".Transcript.FC", sep='')
  y.name <- paste(title,".Protein.FC", sep='')
  g.name <- paste(title,".Regulation", sep='')
  
  df <- df[, c("Name", "top_hit_id", "top_hit_evalue", "top_hit_description", 
               y.name, x.name, g.name)
           ] %>% 
    arrange(desc(!!sym(g.name)))
  
  df$top_hit_evalue <- formatC(df$top_hit_evalue, format = "e", digits = 2)

  print("")
  print(paste("Plotting matrix with",nrow(df),"rows using rows",x.name,",",y.name,"and",g.name))
  
  p <- ggplot(df,aes(x=!!ensym(x.name), y=!!ensym(y.name))) +
    stat_poly_line(se=FALSE, col="red", formula = y ~ x) + # Add trend line
    stat_poly_eq(mapping = use_label(c("n", "R2", "eq")), label.y = 0.95, label.x = 0.05, size=5*text.scale) +   # Add equation
    annotate("text", x=xlim.max-1, y=ylim.max-1, label= "Q1", size=6*text.scale) +
    annotate("text", x=xlim.min+1, y=ylim.max-1, label= "Q2", size=6*text.scale) +
    annotate("text", x=xlim.min+1, y=ylim.min+1, label= "Q3", size=6*text.scale) +
    annotate("text", x=xlim.max-1, y=ylim.min+1, label= "Q4", size=6*text.scale) +
    geom_hline(yintercept = 0, color = "black", size = 0.5*text.scale) +
    geom_vline(xintercept = 0, color = "black", size = 0.5*text.scale) +
    geom_point(aes(col=!!sym(g.name), 
                   label=Name, label2=top_hit_id, label3=top_hit_evalue, label4=top_hit_description
                   ), 
               size=1*text.scale) + # Define 'col' here to stop stat_* from adding trend line for each color/group
    scale_x_continuous(limits=c(xlim.min, xlim.max), expand = c(0, 0), breaks = seq(xlim.min, xlim.max, by = 1)) +
    scale_y_continuous(limits=c(ylim.min, ylim.max), expand = c(0, 0), breaks = seq(ylim.min, ylim.max, by = 1)) +
    labs(title = title,
         x = 'Transcript log2 fold change',
         y = 'Protein log2 fold change',
         tag = tag) +
    scale_color_manual(values=point.color) +
    theme_light() +
    theme(text            = element_text(family="Helvetica"),
          plot.title      = element_text(size=24*text.scale, face="bold", hjust=0.5),
          plot.tag        = element_text(size=24*text.scale, face="bold"),
          axis.title.x    = element_text(size=16*text.scale),
          axis.title.y    = element_text(size=16*text.scale),
          axis.line       = element_line(color="black"),
          axis.line.x     = element_line(color="black"),
          axis.ticks      = element_line(color="black"),
          legend.key.size = unit(1*text.scale, 'cm'),             # Legend key size
          legend.title    = element_text(size=14*text.scale),     # Legend title font size
          legend.text     = element_text(size=12*text.scale),     # Legend text font size
          aspect.ratio=1)                                         # Make plot square
  return(p)
}
```






```r
names <- read.table("Mcapitata_V1_proteomic_data.tsv", sep="\t", header=TRUE, quote="", fill=FALSE, comment.char = '')[,1]
df <- read.table("Mcapitata_V1_multiomics_results.tsv", sep="\t", header=TRUE, quote="", fill=FALSE, comment.char = '')
dim(df)
```

```
## [1] 63227    55
```

```r
df <- df[df$Name %in% names,]
dim(df)
```

```
## [1] 4036   55
```

```r
TP1 <- scatter_plot(df, "TP1", "A", 0.6)
```

```
## [1] ""
## [1] "Plotting matrix with 4036 rows using rows TP1.Transcript.FC , TP1.Protein.FC and TP1.Regulation"
```

```
## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
## ℹ Please use `linewidth` instead.
```

```
## Warning in geom_point(aes(col = !!sym(g.name), label = Name, label2 =
## top_hit_id, : Ignoring unknown aesthetics: label, label2, label3, and label4
```

```r
TP3 <- scatter_plot(df, "TP3", "B", 0.6)
```

```
## [1] ""
## [1] "Plotting matrix with 4036 rows using rows TP3.Transcript.FC , TP3.Protein.FC and TP3.Regulation"
```

```
## Warning in geom_point(aes(col = !!sym(g.name), label = Name, label2 =
## top_hit_id, : Ignoring unknown aesthetics: label, label2, label3, and label4
```

```r
TP5 <- scatter_plot(df, "TP5", "C", 0.6)
```

```
## [1] ""
## [1] "Plotting matrix with 4036 rows using rows TP5.Transcript.FC , TP5.Protein.FC and TP5.Regulation"
```

```
## Warning in geom_point(aes(col = !!sym(g.name), label = Name, label2 =
## top_hit_id, : Ignoring unknown aesthetics: label, label2, label3, and label4
```

```r
out.dir <- "interactive_html_plots"
dir.create(out.dir, showWarnings=FALSE)
saveWidget(ggplotly(TP1), file=paste(out.dir,"prot_vs_trans.FC_TP1.scatter.html", sep='/'))
```

```
## Warning: The following aesthetics were dropped during statistical transformation:
## x_plotlyDomain, y_plotlyDomain
## ℹ This can happen when ggplot fails to infer the correct grouping structure in
##   the data.
## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
##   variable into a factor?
```

```
## Warning: The following aesthetics were dropped during statistical transformation:
## x_plotlyDomain, y_plotlyDomain
## ℹ This can happen when ggplot fails to infer the correct grouping structure in
##   the data.
## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
##   variable into a factor?
```

```
## Warning: `gather_()` was deprecated in tidyr 1.2.0.
## ℹ Please use `gather()` instead.
## ℹ The deprecated feature was likely used in the plotly package.
##   Please report the issue at <https://github.com/plotly/plotly.R/issues>.
```

```
## Warning in geom2trace.default(dots[[1L]][[1L]], dots[[2L]][[1L]], dots[[3L]][[1L]]): geom_GeomTextNpc() has yet to be implemented in plotly.
##   If you'd like to see this geom implemented,
##   Please open an issue with your example code at
##   https://github.com/ropensci/plotly/issues
```

```
## Warning: Aspect ratios aren't yet implemented, but you can manually set a
## suitable height/width

## Warning: Aspect ratios aren't yet implemented, but you can manually set a
## suitable height/width
```

```r
saveWidget(ggplotly(TP3), file=paste(out.dir,"prot_vs_trans.FC_TP3.scatter.html", sep='/'))
```

```
## Warning: The following aesthetics were dropped during statistical transformation:
## x_plotlyDomain, y_plotlyDomain
## ℹ This can happen when ggplot fails to infer the correct grouping structure in
##   the data.
## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
##   variable into a factor?
```

```
## Warning: The following aesthetics were dropped during statistical transformation:
## x_plotlyDomain, y_plotlyDomain
## ℹ This can happen when ggplot fails to infer the correct grouping structure in
##   the data.
## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
##   variable into a factor?
```

```
## Warning in geom2trace.default(dots[[1L]][[1L]], dots[[2L]][[1L]], dots[[3L]][[1L]]): geom_GeomTextNpc() has yet to be implemented in plotly.
##   If you'd like to see this geom implemented,
##   Please open an issue with your example code at
##   https://github.com/ropensci/plotly/issues
```

```
## Warning: Aspect ratios aren't yet implemented, but you can manually set a
## suitable height/width

## Warning: Aspect ratios aren't yet implemented, but you can manually set a
## suitable height/width
```

```r
saveWidget(ggplotly(TP5), file=paste(out.dir,"prot_vs_trans.FC_TP5.scatter.html", sep='/'))
```

```
## Warning: The following aesthetics were dropped during statistical transformation:
## x_plotlyDomain, y_plotlyDomain
## ℹ This can happen when ggplot fails to infer the correct grouping structure in
##   the data.
## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
##   variable into a factor?
```

```
## Warning: The following aesthetics were dropped during statistical transformation:
## x_plotlyDomain, y_plotlyDomain
## ℹ This can happen when ggplot fails to infer the correct grouping structure in
##   the data.
## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
##   variable into a factor?
```

```
## Warning in geom2trace.default(dots[[1L]][[1L]], dots[[2L]][[1L]], dots[[3L]][[1L]]): geom_GeomTextNpc() has yet to be implemented in plotly.
##   If you'd like to see this geom implemented,
##   Please open an issue with your example code at
##   https://github.com/ropensci/plotly/issues
```

```
## Warning: Aspect ratios aren't yet implemented, but you can manually set a
## suitable height/width

## Warning: Aspect ratios aren't yet implemented, but you can manually set a
## suitable height/width
```

```r
ggarrange(TP1 + rremove("legend"), 
          TP3 + rremove("legend"), 
          TP5 + rremove("legend"), 
          get_legend(TP1),
          ncol = 2, nrow = 2)
```

![](prot_vs_trans_FC.scatter_plot_files/figure-html/plot_timepoints-1.png)<!-- -->





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
##  [1] htmlwidgets_1.6.1  plotly_4.10.0      ggpubr_0.4.0       RColorBrewer_1.1-3
##  [5] ggtext_0.1.2       ggpmisc_0.5.2      ggpp_0.5.0         ggplot2_3.4.0     
##  [9] reshape2_1.4.4     dplyr_1.0.10       tibble_3.1.8      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.9         lattice_0.20-45    confintr_0.2.0     tidyr_1.2.1       
##  [5] assertthat_0.2.1   digest_0.6.29      utf8_1.2.2         R6_2.5.1          
##  [9] plyr_1.8.7         backports_1.4.1    MatrixModels_0.5-1 evaluate_0.17     
## [13] highr_0.9          httr_1.4.4         pillar_1.8.1       rlang_1.0.6       
## [17] lazyeval_0.2.2     rstudioapi_0.14    data.table_1.14.2  SparseM_1.81      
## [21] car_3.1-0          jquerylib_0.1.4    Matrix_1.5-1       rmarkdown_2.17    
## [25] splines_4.1.2      polynom_1.4-1      stringr_1.4.1      munsell_0.5.0     
## [29] gridtext_0.1.5     broom_1.0.1        compiler_4.1.2     xfun_0.33         
## [33] pkgconfig_2.0.3    htmltools_0.5.4    tidyselect_1.2.0   viridisLite_0.4.1 
## [37] fansi_1.0.3        withr_2.5.0        MASS_7.3-58.1      grid_4.1.2        
## [41] jsonlite_1.8.2     gtable_0.3.1       lifecycle_1.0.3    DBI_1.1.3         
## [45] magrittr_2.0.3     scales_1.2.1       cli_3.4.1          stringi_1.7.8     
## [49] cachem_1.0.6       carData_3.0-5      farver_2.1.1       ggsignif_0.6.3    
## [53] xml2_1.3.3         bslib_0.4.0        ellipsis_0.3.2     generics_0.1.3    
## [57] vctrs_0.5.2        cowplot_1.1.1      boot_1.3-28        tools_4.1.2       
## [61] glue_1.6.2         purrr_0.3.5        crosstalk_1.2.0    abind_1.4-5       
## [65] fastmap_1.1.0      survival_3.4-0     yaml_2.3.5         colorspace_2.0-3  
## [69] rstatix_0.7.0      knitr_1.40         sass_0.4.2         quantreg_5.94
```
