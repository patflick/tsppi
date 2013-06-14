#!/usr/bin/bash

# filter APE with High or Medium Reliability
cat cerebral.cortex_neuronal.cells.csv | awk -F',' ' ($5 == "\"APE\"" && ($6 == "\"High\"" || $6 == "\"Medium\""))  {print}'
