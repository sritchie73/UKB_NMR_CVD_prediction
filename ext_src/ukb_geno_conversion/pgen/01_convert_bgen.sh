#!/bin/bash

ref_dir=/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/COMMON/post_qc_data/imputed/HRC_UK10K
out_dir=$ref_dir/plink_format/GRCh37/pgen

chr=$SLURM_ARRAY_TASK_ID

if [ $chr = 23 ]; then 
  # X chromosome
  plink2 --bgen $ref_dir/ukb_imp_chrX_v3.bgen \
         --sample $out_dir/dummyX.sample \
         --remove $out_dir/dummyX.remove \
         --threads $SLURM_CPUS_ON_NODE \
         --memory $SLURM_MEM_PER_NODE \
         --silent \
         --make-pgen \
         --out $out_dir/ukb_imp_v3_chrX
elif [ $chr = 24 ]; then
  # XY chromosome
  plink2 --bgen $ref_dir/ukb_imp_chrXY_v3.bgen \
         --sample $out_dir/dummyXY.sample \
         --threads $SLURM_CPUS_ON_NODE \
         --memory $SLURM_MEM_PER_NODE \
         --silent \
         --make-pgen \
         --out $out_dir/ukb_imp_v3_chrXY
else 
  # Chromosomes 1-22
  plink2 --bgen $ref_dir/ukb_imp_chr${chr}_v3.bgen \
         --sample $out_dir/dummy.sample \
         --remove $out_dir/dummy.remove \
         --threads $SLURM_CPUS_ON_NODE \
         --memory $SLURM_MEM_PER_NODE \
         --silent \
         --make-pgen \
         --out $out_dir/ukb_imp_v3_chr${chr}
fi

