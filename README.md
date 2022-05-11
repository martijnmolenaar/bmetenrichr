# bmetenrichr

## metabolite set enrichment using bootstrapping, for single-cell metabolomics data

This R-package aims to perform metabolite set enrichment analysis (MSEA) on single-cell metabolomics datasets.
In contrast to bulk-metabolomics, metabolite annotation is often more ambiguous with fully resolved molecular structures.
That means, annotations are vectors of isomeric (and/or isobaric) molecules, complicating downstream MSEA. This package uses a 
boostrapping approach by performing enrichment analyses many times with random sampling of the isomers/isobars.
