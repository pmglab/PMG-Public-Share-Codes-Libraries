This Python script formats GWAS summary statistics to KGGSEE compatible format, GCTA format and
LDSC format. A SNP reference file contains columns of chromosome, basepair coordinate, SNP ID,
allele 1, allele 2, and allele 1 frequency will be loaded. The input GWAS summary statistics will
be mapped to the SNP reference by either SNP coordinates or SNP IDs. SNPs with any confliction in
coordinate, SNP ID or alleles between the GWAS summary statistics and the SNP reference will be
removed. Allele frequencies and effect sizes will be flipped to match the allele specified by the
`--ref-a1-col` flag. All columns of the output files will be compatible to tht SNP reference file.

The `--sum-file-out` flag specified file will contain the following columns:
```
CHR    chromosome
BP     basepair coordinate
SNP    SNP ID
A1     the effect allele
A2     the other allele
FRQ    frequency of A1
BETA   effect size of A1
SE     stderr of effect size
P      p-value
Neff   effective sample size
```

The `--gcta-out` flag specified output file will contain columns of SNP, A1, A2, FRQ, BETA, SE, P and
Neff. The `--ldsc-out` flag specified output file will contain columns of SNP, A1, A2, N and Z. These
two outputs can be masked by booleanizable-value columns (e.g., 0 for False and 1 for True) in the
SNP reference file.