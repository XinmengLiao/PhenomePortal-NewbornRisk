Family-based anlaysis could be only be conducted for a single family by `gatk` or `gatk4` installed with conda.

####  GATK PhaseByTransmission
```bash
gatk PhaseByTransmission \
  -R /mnt/storage_pool/Genomics/VarXOmics/tools/picard/hg38.fa \
  -V F1_merged_rmmissingalt_biallelic.vcf.gz \
  --pedigree F1.ped \
  -O F1_gatk.vcf.gz \
  --default-qual 30.0
```
