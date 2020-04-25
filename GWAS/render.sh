#!/usr/bin/env bash
#render.sh
R -e "rmarkdown::render('runPlink.Rmd', output_format='all')"
