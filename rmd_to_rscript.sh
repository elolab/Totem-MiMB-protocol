#!/bin/bash

# Convert Rmd notebooks into R scripts
R --no-echo --no-restore --no-save -e "knitr::purl(input = 'notebooks/QuickStart_Totem.Rmd', output = 'scripts/QuickStart_Totem.R', documentation=0L)"
R --no-echo --no-restore --no-save -e "knitr::purl(input = 'notebooks/GuidedStart_Totem.Rmd', output = 'scripts/GuidedStart_Totem.R', documentation=0L)"

