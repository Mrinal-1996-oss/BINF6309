#!/usr/bin/env bash
# renderMerge.sh
R -e "rmarkdown::render('lungData.Rmd', output_format='all')"
