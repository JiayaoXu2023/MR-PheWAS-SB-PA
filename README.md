# MR-PheWAS-SB-PA
This repository provides code for the paper "The effect of sedentary behaviour and physical activity on 1719 diseases: a Mendelian randomisation phenome-wide association study (MR-PheWAS)"

## 1. Prepare exposure data
1_exposure data.R demonstrates how we selected genetic instruments for the exposures sedentary behaviours (LST) and physical activity (MVPA). 
SNPs for LST and MVPA were selected based on a genome-wide significance threshold (p < 5 × 10⁻⁸) and clumped using a linkage disequilibrium (LD) threshold of r² < 0.001 within a 10,000 kb window, based on the 1000 Genomes European reference panel.

## 2. Prepare outcome data
2_outcome data.R demonstrates how we selected genetic instruments from [FinnGen](https://finngen.gitbook.io/documentation/r11) GWAS summary statistics and formatted the data.

## 3. Run main analyses (MR PheWAS)
3_main analysis.R demonstrates how we performed MR-PheWAS using the inverse variance weighted (IVW) method.

## 4. Sensitivity analyses
4_sensitivity analysis demonstrates how we performed Cochran’s Q test, MR weighted median (WM), MR-Egger and Steiger filtering methods.

## 5. Replication analyses
Results from the main analyses were replicated using available GWAS data from UK Biobank.

5_1_replication neale.R demonstrates how we replicated results using GWAS summary statistics from [Nealelab](https://www.nealelab.is/uk-biobank).

5_2_replication geneatlas.R demonstrates how we replicated results using GWAS summary statistics from [GeneATLAS](http://geneatlas.roslin.ed.ac.uk/downloads/).

5_3_replication opengwas.R demonstrates how we replicated results using GWAS summary statistics from [OpenGWAS](https://opengwas.io/datasets/).

