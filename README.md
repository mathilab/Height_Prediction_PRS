# Polygenic scores for height in admixed populations

[Link to manuscript](https://www.g3journal.org/content/early/2020/09/02/g3.120.401658)

## Contents

# Table of contents
1. [Introduction](#intro)
2. [Clone this repo and recreate analyses](#recreating-all-analyses-in-the-paper)
	1. [Prerequisites](#prerequisites)
	2. [Preparing input data](#prepare-genotype-input-data)
	3. [C+T and PRS](#clump/prune/calculate-prs)
	4. [Plots](#plots)
	5. [LDpred](#ldpred-analyses)
	6. [Stratified PRS by bins of recombination rate](#stratified-prs)
	7. [PCA and GWAS in the UKB indivudals with admixed African ancestry](#pca-and-gwas-in-ukb-individuals-with-admixed-African-ancestry)
	8. [PRS incorporating local ancestry](#prs-incorporating-local-ancestry)
	9. [Effect size vs frequencies across ancestries](#effect-size-differences-accross-ancestries-and-their-relationship-with-frequency-differences)
3. [Skip details. Take me to the scripts to recreate tables and figures in the paper](#make-all-figures-and-supplementary-tables-in-the-paper)

## Intro

Height is a very polygenic trait and well-studied in humans. GWAS summary statistics for height used hundreds of thousands of individuals of European ancestry. It is unclear how well polygenic risk scores (PRS) predict height in non-Europeans in comparison to Europeans. Here we explore and quantify this, and then provide improvements.

## Recreating all analyses in the paper

### Prerequisites

*R and packages listed in the scripts.
*shell
*python (for ldpred1 analyses)
*ldpred
*PLINK and PLINK2
*bcftools
*there is an assumption that you have access to the datasets used in this paper. We are not allowed to share the raw data.


Go to a direcotry in your computer and clone this repo:

```
git clone https://github.com/mathilab/Height_Prediction_PRS.git
#from now on it is assumed your project root directory is called "Height_Prediction_PRS".
```

### Prepare genotype input data

Go to the [input directory](input/README.md) and follow instructions. 

Within each dataset's directory you will find a README.md with instructions on how to prepare input data. 

*Note:* if you do get access to the data you will need to fix the path to the data accordingly. 


### Clump/prune/calculate PRS

Pruning/clumping using both the GWAS effect sizes ('gwas') and the sibling-estimated effect sizes ('sib_betas'). In both cases, p-values used for clumping come from the full UKB GWAS.

**Using GWAS summmary statistics** 

Go to [gwas](gwas/README.md) and follow instructions. 

**Using sibling pairs summmary statistics**

Go to [sib_betas](sib_betas/README.md) and follow instructions.

### Combine PRS results

Combine all PRS results per dataset:
```
for D in JHS WHI ukb_afr ukb_eur HRS_eur HRS_afr;
do
for F in gwas sib_betas unweighted_prs;
do 
rm $F/${D}/scripts/test2.txt;
rm $F/${D}/scripts/run_this_PGS.sh;
scripts/combine_Rds_PGS.sh $F $D
scripts/combine_Rds_PGS.sh gwas $D
unweighted_prs/combine_Rds_PGS.sh unweighted_prs $D #need to fix this
done
done

```

### Plots

These scripts will produce plots for each pruning/clumping strategy. Throughout the paper we show the one called "phys_100000_0.0005":

```
for J in gwas sib_betas unweighted_prs;
do
Rscript --vanilla ${J}/WHI/scripts/Plots_WHI.R
Rscript --vanilla ${J}/JHS/scripts/Plots_JHS.R
Rscript --vanilla ${J}/ukb_afr/scripts/Plots_ukb_afr.R
Rscript --vanilla ${J}/ukb_eur/scripts/Plots_ukb_eur.R
Rscript --vanilla ${J}/HRS_afr/scripts/Plots_HRS_afr.R
Rscript --vanilla ${J}/HRS_eur/scripts/Plots_HRS_eur.R
done
for J in gwas sib_betas;
do
Rscript --vanilla scripts/combine_datasets.R ${J
done
Rscript --vanilla unweighted_prs/combine_datasets.R
```


### LDpred analyses

See [README.md in the ldpred directory](ldpred/README.md)


### Stratified PRS

PRS for SNPs in bins of recombination rate.

See [README.md in the stratified prs directory](strat_prs/README.md)


### PCA and GWAS in UKB individuals with admixed African ancestry


See [README.md in the PCA and GWAS directory](PCA_and_GWAS/UKB_AFR_imputed/README.md)


### PRS incorporating local ancestry

See [README.md in the local ancestry directory](loc_anc_analyses/README.md)




## Effect size differences accross ancestries and their relationshop with frequency differences


See [README.md in the epistasis directory](epistasis/README.md)


## Make all figures and supplementary tables in the paper



```
#Figure 1
Rscript --vanilla scripts/Fig1.R
#Figure 2
Rscript --vanilla scripts/Fig2.R
#Figures 3, S11, S14
Rscript --vanilla scripts/qqplot_betadiff.R
Rscript --vanilla scripts/Fig3.R gwas AA 20000 WHI  #Fig 3 & Fig S11
Rscript --vanilla scripts/Fig3.R gwas CEU 20000 WHI  Fig S14
Rscript --vanilla scripts/Fig4.R
Rscript --vanilla scripts/Fig5.R
Rscript --vanilla scripts/make_all_tables.R
Rscript --vanilla scripts/FigS7.R
```

*Note:* I have tried my best to make this research reproducible. Please do reach out if you have any suggestions on how to improve that if you have any trouble.

