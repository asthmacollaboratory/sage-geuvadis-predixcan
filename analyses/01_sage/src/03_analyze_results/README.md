# Analyze PrediXcan predictions in SAGE

Compile summary statistics comparing gene expression predictions versus RNA-Seq measurements in SAGE.

## Prerequisites
The analysis requires the following R packages:
- `data.table`
- `purrr`
- `broom`
- `ggplot2`
- `dplyr`
- `optparse`

Install these with any standard technique, such as
```R
install.packages(c("data.table","purrr","broom","ggplot2","dplyr","optparse")
```

## Running
```bash
./analyze_predictions.sh
```

## Notes

This script generates output organized into 3 categories:
- plots
- tables
- `Rdata`

The plots are saved as PNG format by default. Changing the default plot format, e.g. to PDF, requires editing `analyze_prediction.sh` accordingly.
The tables contain summary statistics only. They are saved as `.txt` files.
The `Rdata` object contains all generated output. Load it from within R via
```R
load("/PATH/TO/OUTPUT/FILES/sage_predixcan_allresults_allplots.Rdata")
```
