## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------

options(stringsAsFactors = FALSE, warn = -1)

## install devtools if not installed
if(!("devtools" %in% rownames(installed.packages()))){
  install.packages("devtools",  repos = c(CRAN = "http://cran.rstudio.com"))
  }

## install bmetenrichr if not installed
if(!("bmetenrichr" %in% rownames(installed.packages()))){
   devtools::install_github(repo = "martijnmolenaar/bmetenrichr", build_vignettes = TRUE)
}

library(bmetenrichr)

## -----------------------------------------------------------------------------
data("Rappez_et_al")

## the main input is a single-cell metabolomics matrix with molecules as rows and cells as columns
Rappez_et_al$sc_matrix[1:10,1:10]

## for the molecules, a vector with molecular formulas plus adduct is required
rownames(Rappez_et_al$sc_matrix)[1:10]

## a conditions vector is required to define to which condition a given cell belongs
Rappez_et_al$conditions[1:10]

## in this analysis, only specific annotations should be included as others are extracellular 
Rappez_et_al$cellular[1:10]


## -----------------------------------------------------------------------------

myTestRun <-
  initEnrichment(scmatrix = Rappez_et_al$sc_matrix,
                 annotations = rownames(Rappez_et_al$sc_matrix),
                 conditions = Rappez_et_al$conditions,
                 include = Rappez_et_al$cellular,
                 condition.x = "U",
                 condition.y = "F"                    )

## ----fig.height=4.5, fig.width=8, warning=FALSE-------------------------------
## rank metabolites, in this case by t.test statistic

myTestRun <- rankScore(myTestRun, ranking.by = 't.test')


## ----fig.height=4.5, fig.width=8, warning=FALSE-------------------------------

myTestRun <- calcEnrichment(myTestRun, n = 100)


## ----fig.height=4.5, fig.width=8, warning=FALSE-------------------------------

plotEnrichment(myTestRun, min.annotations = 5, q.value.cutoff = .1, by.statistic = "ES")

## ----fig.height=4.5, fig.width=8, warning=FALSE-------------------------------

plotEnrichment(myTestRun, min.annotations = 5, q.value.cutoff = .05, plotIDs = T, 
               by.statistic = "q.value")


## ----fig.height=4.5, fig.width=8, warning=FALSE-------------------------------

enrichmentTable(myTestRun)[1:10,]


## ----fig.height=4.5, fig.width=8, warning=FALSE-------------------------------

## now, let's test FIT vs F

myTestRun <-  setConditions(object = myTestRun, condition.x = 'F', condition.y = 'FIT')
myTestRun

## rank metabolites, in this case by t.test statistic

myTestRun <- rankScore(myTestRun, ranking.by = 't.test')

## and perform enrichment analysis

myTestRun <- calcEnrichment(myTestRun, n = 100)

plotEnrichment(myTestRun, min.annotations = 5, q.value.cutoff = .05)

## ----fig.height=4.5, fig.width=8, warning=FALSE-------------------------------

## create object

myTestRun <-
  initEnrichment(scmatrix = Rappez_et_al$sc_matrix, 
                 isobars = TRUE,                      ## to include isobars (default is FALSE)
                 mass_range_ppm = 3,                  ## mass range to define isobars
                 polarization_mode = "positive",      ## mode is important to include the right adducts
                 annotations = rownames(Rappez_et_al$sc_matrix), 
                 conditions = Rappez_et_al$conditions,
                 include = Rappez_et_al$cellular,
                 condition.x = "U",
                 condition.y = "F"                    )

## ----fig.height=4.5, fig.width=8, warning=FALSE-------------------------------
## rank metabolites, in this case by t.test statistic

myTestRun <- rankScore(myTestRun, ranking.by = 't.test')

## perform enrichment analysis with n = 100 bootstraps

myTestRun <- calcEnrichment(myTestRun, n = 100)

## example of the annotations, that now also include isobars

myTestRun$annotations[[
  sample(which(sapply(myTestRun$isobars_list, length) > 1), size = 1)]][1:10]

## plot enrichment analysis, with q.values on x-axis, and with LION IDs

plotEnrichment(myTestRun, min.annotations = 5, q.value.cutoff = .05, 
               by.statistic = "q.value")

