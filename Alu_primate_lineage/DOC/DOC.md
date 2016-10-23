# Alu exon evolution on primate lineage

### Example usage

Just run `bash run.sh bedfile.bed` to go over all the programs and plot figures

## Software version

Pipeline was develop and tested on:

GNU bash, version 4.3.11(1)-release (x86_64-pc-linux-gnu)

Python 2.7.6 [GCC 4.8.2] on linux2

Pycharm IDE PyCharm 2016.1.2
Build #PY-145.844, built on April 8, 2016
JRE: 1.8.0_40-release-b132 x86_64
JVM: OpenJDK 64-Bit Server VM by JetBrains s.r.o

R version 3.3.1 (2016-06-21) -- "Bug in Your Hair"
Platform: x86_64-pc-linux-gnu (64-bit)

Rstudio IDE Version 0.99.489 – © 2009-2015 RStudio, Inc.
Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/538.1 (KHTML, like Gecko) RStudio Safari/538.1 Qt/5.4.0

## Dependencies

The folowing programs must be on the path to be called by the pipeline:

1.[Bedtools v2.25.0]( http://bedtools.readthedocs.io/en/latest/)

2.[ MaxEntScan  ]( https://github.com/Congenica/maxentscan )
NOTE!! splicemodels folder  is needed to run the program. Refer to [ MaxEntScan  ]( https://github.com/Congenica/maxentscan ) for instalation and usage.


3.[ R  ] (https://www.r-project.org/) Libraries must be available for R

    library(vcd)
    library(gmodels)
    library(plyr)
    library(pheatmap)
    library(ggplot2)
    require(gplots)
    require(reshape)
    library(system)
    library(parallel)
    library(GMD)
    library(RColorBrewer)
    library(gridGraphics)
    library(ColorPalette)
    library(scales)
    library(GGally)
    library(sys)
    library(Bio) import SeqIO




