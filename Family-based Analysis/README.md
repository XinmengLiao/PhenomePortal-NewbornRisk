Family based analysis will be performed simply by Python and R scripts, since the GATK and slivar always fail. 

#### Input files
1. Family merged vcf.gz file
2. Family ped file

Command for analysis
```bash
newbornrisk='/mnt/storage_pool/Genomics/Genome/NewbornRisk'
bash NewbornRisk_family.sh \
  -i F4 \
  -o $newbornrisk/examples/F4 \
  -v $newbornrisk/examples/F4/F4_merged.vcf.gz \
  --ped $newbornrisk/examples/F4/F4.ped \
  --only-pass yes --genome GRCH38 --genedb TR
```

<img width="1920" height="968" alt="image" src="https://github.com/user-attachments/assets/e5e2181f-5fee-4816-a1d9-479c980d8f26" />
<img width="1914" height="737" alt="image" src="https://github.com/user-attachments/assets/e662948e-4134-4d15-8f65-dc333222be90" />
<img width="1917" height="732" alt="image" src="https://github.com/user-attachments/assets/ff26aba1-6aff-462f-b903-83b59b456f69" />

