# tRNAdbImport <img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/tRNA/tRNA.png" height="200" align="right">

<!-- badges: start -->
![R-CMD-check](https://github.com/FelixErnst/tRNAdbImport/workflows/R-CMD-check/badge.svg)
[![codecov](https://codecov.io/gh/FelixErnst/tRNAdbImport/branch/master/graph/badge.svg)](https://codecov.io/gh/FelixErnst/tRNAdbImport)
<!-- badges: end -->

The tRNAdb and mttRNAdb ([Juehling et al. 2009](#Literature)) is a compilation of
tRNA sequences and tRNA genes. It is a follow up version of the database of
[Sprinzl et al. (2005)](#Literature).
Using `tRNAdbImport` the tRNAdb can be accessed as outlined on the website
[trna.bioinf.uni-leipzig.de](trna.bioinf.uni-leipzig.de) directly via R. The
results are returned as a `GRanges` object and can be further used in a
BioC context.

# Installation

The current version of the `tRNAdbImport` package is available from Bioconductor.
 
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("tRNAdbImport")
# Load and attach the package
library("tRNAdbImport")
```

## Literature

- Jühling, Frank; Mörl, Mario; Hartmann, Roland K.; Sprinzl, Mathias; Stadler,
Peter F.; Pütz, Joern (2009): "TRNAdb 2009: Compilation of tRNA Sequences and
tRNA Genes." Nucleic Acids Research 37 (suppl_1): D159–D162.
doi:[10.1093/nar/gkn772](https://doi.org/10.1093/nar/gkn772). 
- Sprinzl, Mathias; Vassilenko, Konstantin S. (2005): "Compilation of tRNA 
Sequences and Sequences of tRNA Genes." Nucleic Acids Research 33 (suppl_1): 
D139–D140. doi:[10.1093/nar/gki012](https://doi.org/10.1093/nar/gki012).
