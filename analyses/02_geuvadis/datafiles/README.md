# Data files for analysis
This directory includes various small data files used for downloading, parsing, and analyzing GEUVADIS data.

## Files
| Name | Description | Used in |
| --- | --- | --- |
| `geuvadis.yri89.sampleids.txt` | A list of 89 YRI sample IDs used for parsing GEUVADIS files | [`download_parse_data.sh`](../src/01_download_parse_data/download_parse_data.sh) | 
| `geuvadis.eur373.sampleids.txt` | A list of 373 EUR sample IDs (CEU + GBR + TSI + FIN) used for parsing GEUVADIS files | [`download_parse_data.sh`](../src/01_download_parse_data/download_parse_data.sh) | 
| `genelist_plus_500kb` | A simple 4-column list with gene info (Ensembl ID, chr, start, end), with 500kb padding on each end | [`qsub_geuvadis_compute_weights.sh`](../src/common/qsub_geuvadis_compute_weights.sh) |
