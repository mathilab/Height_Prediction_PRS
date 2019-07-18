
*copy files to my local directory*
```
cp /project/mathilab/data/JHS/data/JHS_b37_strand_prune.* .
cp /project/mathilab/data/JHS/data/JHS_phenotypes.txt .
cp /project/mathilab/data/JHS/admixture/JHS_b37_strand_prune.2.Q . ##order comes from file below
cp /project/mathilab/data/JHS/data/WHI_b37_strand_prune_include.fam .  #the order of samples
```
*convert to vcf format*

```
plink --allow-extra-chr --allow-no-sex --autosome \
--bfile JHS_b37_strand_prune \
--recode vcf \
--keep-allele-order \
--out JHS_b37_strand_prune.bas
```

```
awk '/^ *#/ { print; }' JHS_b37_strand_prune.bas.vcf > header.txt

for chr in {1..22};
do
touch chr${chr}_bas.vcf
cat header.txt >  chr${chr}_bas.vcf
awk -v i=$chr '$1==i' JHS_b37_strand_prune.bas.vcf >> chr${chr}_bas.vcf
echo 'chr'
echo $chr
echo 'done'
done
```

*compress/index*
```
for chr in {1..22};
do
bgzip -c chr${chr}_bas.vcf > chr${chr}_bas.vcf.gz &&
tabix -p vcf chr${chr}_bas.vcf.gz
done

bgzip -c JHS_b37_strand_prune.bas.vcf > JHS_b37_strand_prune.bas.vcf.gz