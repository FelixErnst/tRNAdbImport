# tRNAdbImport

The tRNAdb and mttRNAdb [Juehling et al. 2009](#Literature) is a compilation of
tRNA sequences and tRNA genes. It is a follow up version of the database of
Sprinzl et al. [Sprinzl et al. 2005](#Literature).
Using `tRNAdbImport` the tRNAdb can be accessed as outlined on the website
[trna.bioinf.uni-leipzig.de](trna.bioinf.uni-leipzig.de) directly via R and the
results are returned as a `GRanges` object. The results can be further used in
BioC context.

## Installation

The current version of `tRNAscanImport` is available from GitHub.
 
```{r}
devtools::install_github("FelixErnst/tRNA")
devtools::install_github("FelixErnst/tRNAdbImport")
# Load and attach thepackage
library("tRNAdbImport")
# this also attaches 
```
A submission to Bioconductor might follow

## Literature

- Jühling, Frank, Mario Mörl, Roland K. Hartmann, Mathias Sprinzl, Peter F.
Stadler, and Joern Pütz. 2009. “TRNAdb 2009: Compilation of tRNA Sequences and
tRNA Genes.” Nucleic Acids Research 37 (suppl_1): D159–D162.
doi:[10.1093/nar/gkn772](https://doi.org/10.1093/nar/gkn772). - Sprinzl,
Mathias, and Konstantin S. Vassilenko. 2005. “Compilation of tRNA Sequences and
Sequences of tRNA Genes.” Nucleic Acids Research 33 (suppl_1): D139–D140.
doi:[10.1093/nar/gki012](https://doi.org/10.1093/nar/gki012).