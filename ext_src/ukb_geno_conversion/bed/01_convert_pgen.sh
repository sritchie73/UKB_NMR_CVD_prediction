#!/bin/bash

ref_dir=/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/COMMON/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen
out_dir=/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/COMMON/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/bed

chr=$SLURM_ARRAY_TASK_ID
if [ $chr -eq 23 ]; then
  chr="X"
elif [ $chr -eq 24 ]; then
  chr="XY"
fi

# also compute hard call frequency information
plink2 --pfile $ref_dir/ukb_imp_v3_dedup_chr${chr} \
       --geno-counts \
       --make-bed \
       --out $out_dir/ukb_imp_v3_dedup_chr${chr} \
       --threads $SLURM_CPUS_ON_NODE \
       --memory $SLURM_MEM_PER_NODE \
       --silent 
