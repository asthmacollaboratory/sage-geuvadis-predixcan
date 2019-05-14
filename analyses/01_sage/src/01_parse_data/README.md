# Parse SAGE data

Parse and format genotype and expression data from SAGE for analysis with PrediXcan

## Prerequisites
The following R libraries are used:
-- `argparser`
-- `peer`
-- `data.table`
-- `optparse`
-- `preprocessCore`

`preprocessCore` must be installed from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/preprocessCore.html).
The other packages can be installed from within R via normal means, such as
```R
install.packages(c("argparser", "peer", "data.table", optparse"))
```

## Running
```bash
./parse_sage_data.sh
```

## Notes
The script `run_PEER.R` is slightly modified from code by François Aguet. The original code is [here](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/run_PEER.R).
