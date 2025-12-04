### Overview
NewbornRisk location on server: `/mnt/storage_pool/Genomics/Genome/NewbornRisk`
1. Script for family analysis: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/NewbornRisk_family.sh`
2. Script for manage the annotated variants from VEP: ``
3. Script for manage the python result: `Pedigree_Analysis20251203.R`
4. Example output files: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/examples/F4`
5. Input file for example:
   - Merged VCF file: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/examples/F4/F4_merged.vcf.gz`
   - Family information PED file: `/mnt/storage_pool/Genomics/Genome/NewbornRisk/examples/F4/F4.ped`
  
### Command for analysis
```bash
newbornrisk='/mnt/storage_pool/Genomics/Genome/NewbornRisk'
bash NewbornRisk_family.sh \
  -i F4 \
  -o $newbornrisk/examples/F4 \
  -v $newbornrisk/examples/F4/F4_merged.vcf.gz \
  --ped $newbornrisk/examples/F4/F4.ped \
  --only-pass yes --genome GRCH38 --genedb TR \
  --fork 20 --threads 20
```

### UI design
It will be nice if the following contents could be included in the webpage:
1) Couple details and variant summary.
2) Potential recessive disease that will be inherited to the offsprings.
3) PGS information.
4) PGx information.

#### 1) Couple details and variant summary.
1. The table of family details: `family_information.txt`
2. The table of variant summary: `variant_summary.txt`

#### 2) Potential recessive disease that will be inherited to the offspring
1. Table of potential recessive diseases: `carrier_screening_result_couple.txt`
* Note: Only diseases caused by ClinVar or ACMG P/LP variants will be shown as default. 

#### 3) PGS information (have not done yet)
1. Table of the family PGS scores: `cohort_PGS.txt`
2. Density plot compared with the reference populations: `PGS_Znorm1_DensityPlot.png` and `PGS_Znorm2_DensityPlot.png`
3. Report for the PGS score: `report.html`

#### 4) PGx information (have not done yet)
Results from PharmCat in HTML and/or JSON file.
