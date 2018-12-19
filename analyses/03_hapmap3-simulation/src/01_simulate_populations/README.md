# Simulate populations

This folder contains code to simulate haplotypes from [HapMap3 data](https://mathgen.stats.ox.ac.uk/impute/impute_v1.html#Using_IMPUTE_with_the_HapMap_Data) using [HAPGEN2](http://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html).

## Analysis steps
1. Download and unpack both HAPGEN2 and the HapMap3 data.
2. Download gene positions from [Ensembl](https://uswest.ensembl.org/index.html) for _Homo sapiens_. By default, this analysis uses genes from chromosome 22, though any autosome will do.
3. Using HapMap3 haplotypes for CEU and YRI, use HAPGEN2 to forward-simulate 1000 individuals from each of CEU and YRI.
4. Construct an admixed "African American" population by sampling from the simulated haplotypes. This analysis uses admixture proportions of 80\% YRI haplotypes and 20\% CEU haplotypes.
5. Format resulting haplotypes into genotype files. This anlaysis produces files in [PLINK RAW](https://www.cog-genomics.org/plink2/formats#raw) format. 

See comments in `simulate_populations.sh`.

## Simulation 
From a Bourne shell, execute
```
./simulate_populations
```

## Notes
* HAPGEN2 does not allow users to set random seeds. As a result, simulated haplotypes are most likely **not** reproducible. 

## References
* Su, Z., Marchini, J., and Donnelly, P. (2011). HAPGEN2: simulation of multiple disease SNPs. _Bioinformatics_ **27**, 2304-2305.
* Baharian, S., Barakatt, M., Gignoux, C.R., Shringarpure, S., Errington, J., Blot, W.J., Bustamante, C.D., Kenny, E.E., Williams, S.M., Aldrich, M.C., et al. (2016). The Great Migration and African-American Genomic Diversity. _PLOS Genetics_ **12**, e1006059.
