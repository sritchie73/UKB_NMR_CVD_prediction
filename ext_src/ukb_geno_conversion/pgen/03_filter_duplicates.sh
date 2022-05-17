#!/bin/bash

out_dir=/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/COMMON/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen
chr=$SLURM_ARRAY_TASK_ID

if [ $chr -eq 23 ]; then
  chr="X"
elif [ $chr -eq 24 ]; then
  chr="XY"
fi

if [ -f $out_dir/chr${chr}_duplicates.txt ]; then
  # Exclude duplicates and create new pgen files
  plink2 --pfile $out_dir/ukb_imp_v3_chr${chr} \
         --exclude $out_dir/chr${chr}_duplicates.txt \
         --threads $SLURM_CPUS_ON_NODE \
         --memory $SLURM_MEM_PER_NODE \
         --silent \
         --make-pgen \
         --out $out_dir/ukb_imp_v3_dedup_chr${chr}

  # Remove old pgen files.
  rm $out_dir/ukb_imp_v3_chr${chr}.pgen
  rm $out_dir/ukb_imp_v3_chr${chr}.pvar
  rm $out_dir/ukb_imp_v3_chr${chr}.psam

  # Remove temporary exclusion file
  rm $out_dir/chr${chr}_duplicates.txt
else
  # No duplicates to remove, just mv the files.
  mv $out_dir/ukb_imp_v3_chr${chr}.pgen $out_dir/ukb_imp_v3_dedup_chr${chr}.pgen
  mv $out_dir/ukb_imp_v3_chr${chr}.pvar $out_dir/ukb_imp_v3_dedup_chr${chr}.pvar
  mv $out_dir/ukb_imp_v3_chr${chr}.psam $out_dir/ukb_imp_v3_dedup_chr${chr}.psam
  echo "There were no duplicate variants to remove." > $out_dir/ukb_imp_v3_dedup_chr${chr}.log
fi

