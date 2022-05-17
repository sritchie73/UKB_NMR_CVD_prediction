# README
--------

Author: Scott Ritchie <sr827@medschl.cam.ac.uk>
Date: 20-04-2020

--------

These files contain the UK Biobank genotype data in plink version 2's binary 
format (https://www.cog-genomics.org/plink/2.0/input#pgen). This format supports:

    - Reliable tracking of REF vs. ALT alleles.
    - Computationally efficient compression of low-MAF variants and high-LD 
      variant pairs.
    - Phased genotypes.
    - Dosages.
    - VCF-style header information (including species-specific chromosome 
      info, so you don't have to constantly use --chr-set).
    - Multiallelic variants.
    - Multiple phenotypes.
    - Named categorical phenotypes (a phenotype string which doesn't start 
      with a number is interpreted as a category name).

Notably, it preserves probabilistic dosage information present in the BGEN files,
while allowing for rapid and memory efficient access to the data by the plink2 
program.

Some additional QC has also been performed on the imputed genotype data, 
specifically I've removed duplicate variants, keeping from each pair/set
the variant with the highest imputation INFO score. There are typically a
handful of these variants in each bgen file (<100) which can cause problems
with many plink functions crashing. 

Note multi-allelic variants have been kept and are spread across multiple 
rows in the variant information and genotype data files (this is also true 
in the BGEN files).If you are working with these, you may want to create 
a new copy of the .pvar files with the ID column modified to always include
the REF and ALT alleles so that plink does not otherwise crash when 
encountering duplicate variants (i.e. because multiple rows will share the
same rsid for multi-allelic variants).

If you require functionality present in plink version 1.9 but not yet present 
in plink version 2, we also provide files in the plink 1 binary format under
(../bed/)

