#!/bin/bash

RMDFILE=$1
echo $RMDFILE

Rscript -e "require(knitr); require(markdown); rmarkdown::render('${RMDFILE}.Rmd', output_format = 'all')"