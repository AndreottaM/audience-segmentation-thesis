# MM-audiencesegmentation

Audience segmentation of the Australian public. Contains stimuli, data, and scripts for Chapter 4 of the Matthew Andreotta's thesis. Co-authored by Mark Hurlstone, Fabio Boschetti, Simon Farrell, Cecile Paris, and Iain Walker. Any questions can be sent to matthewandreotta@gmail.com.

**Please note the study numbers of this repository do not match the study numbers in the text.** Study 1 of Chapter 4 is Study 2 of this repository. Study 2 of Chapter 4 is Study 3 of this repository. The pilot study reported in the Supplementary Materials is Study 1 of this repository.

Studies were approved by the University of Western Australia Human Research Ethics Office (RA/4/20/5104), and reciprocated by the Commonwealth Scientific and Industrial Research Organisation Human Research Ethics Committee (026/19).

Study 2 was pre-registered using the Open Science Framework (https://osf.io/e7zhx/).

## Repository Structure

### Studies

**study-1**: Free list study to pilot the mental model scales.

**study-2**: Audience segmentation Q sort alongside surveys of psychological characteristics (including mental models).

**study-3**: Audience segmentation Q sort alongside belief-updating tasks.

### File Structure

**codebook**: For qualitative data analysis, codebooks were required to standardise participant responses. There are only codebooks for study-1 and study-4.

**data**: Rawest version of anonymised data available.

**out**: Figures and output files for examination or inclusion in a manuscript.

**paper**: Assets from chapter.

**relevant-notes**: Ethics documentation and collated notes. The latter is unlikely to be of use to anyone.

**scripts**: Analysis scripts.

**study1, study2, ..., studyX**: Stimuli and Qualtrics surveys (in .QSF file for easy import) for each study.


## Requirements

### R

R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

#### Attached packages

ggpubr (0.3.0), kableExtra (1.1.0), lme4 (1.1-17), doParallel (1.0.14), iterators (1.0.10), glmnet (2.0-18), foreach (1.4.4), Matrix (1.2-14), forcats (0.5.0), stringr (1.4.0), dplyr (0.8.5), purrr (0.3.4), readr (1.3.1), tidyr (1.0.2), tibble (3.0.1), ggplot2 (3.3.1), tidyverse (1.3.0), RGenData (1.0), psych (1.8.4), qmethod (1.5.4).

#### All packages
ggpubr (0.3.0), kableExtra (1.1.0), lme4 (1.1-17), doParallel (1.0.14), iterators (1.0.10), glmnet (2.0-18), foreach (1.4.4), Matrix (1.2-14), forcats (0.5.0), stringr (1.4.0), dplyr (0.8.5), purrr (0.3.4), readr (1.3.1), tidyr (1.0.2), tibble (3.0.1), ggplot2 (3.3.1), tidyverse (1.3.0), RGenData (1.0), psych (1.8.4), qmethod (1.5.4), httr (1.4.1), jsonlite (1.6.1), viridisLite (0.3.0), splines (3.5.0), carData (3.0-3), modelr (0.1.8), assertthat (0.2.1), blob (1.2.1), cellranger (1.1.0), yaml (2.1.19), pillar (1.4.4), backports (1.1.6), lattice (0.20-35), glue (1.4.1), digest (0.6.25), ggsignif (0.6.0), rvest (0.3.5), GPArotation (2014.11-1), minqa (1.2.4), colorspace (1.4-1), htmltools (0.4.0), pkgconfig (2.0.3), broom (0.5.6), haven (2.2.0), xtable (1.8-3), scales (1.1.1), webshot (0.5.1), openxlsx (4.1.4), rio (0.5.16), car (3.0-7), generics (0.0.2), ellipsis (0.3.1), withr (2.2.0), cli (2.0.2), mnormt (1.5-5), magrittr (1.5), crayon (1.3.4), readxl (1.3.1), evaluate (0.14), fs (1.4.1), fansi (0.4.1), nlme (3.1-137), MASS (7.3-49), rstatix (0.6.0), xml2 (1.3.2), foreign (0.8-71), data.table (1.11.6), tools (3.5.0), hms (0.5.3), lifecycle (0.2.0), munsell (0.5.0), reprex (0.3.0), zip (2.0.4), compiler (3.5.0), rlang (0.4.6), grid (3.5.0), nloptr (1.0.4), rstudioapi (0.11), rmarkdown (2.2), gtable (0.3.0), codetools (0.2-15), abind (1.4-5), curl (3.2), DBI (1.0.0), R6 (2.4.1), lubridate (1.7.4), knitr (1.28), stringi (1.1.7), Rcpp (1.0.4.6), vctrs (0.2.4), dbplyr (1.4.4), tidyselect (1.1.0), xfun (0.14), doRNG (1.8.2), rngtools (1.5)

### TeX Installation

All LaTeX documents are generated using R sweave with a local TeX installation. All documents were compiled using XeLaTex.

### SVG

All .SVG files were created with Inkscape and converted to .tex files.
