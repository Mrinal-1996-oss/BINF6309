#!/usr/bin/env bash
# renderMerge.sh
R -e "rmarkdown::render('mergeKo.Rmd', output_format='md')"
