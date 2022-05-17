# README
--------

Author: Scott Ritchie <sr827@medschl.cam.ac.uk>
Date: 20-04-2020

--------

These files contain the UK Biobank genotype data in plink version 1's binary 
format (https://www.cog-genomics.org/plink/2.0/input#bed). 

Note that this format stores the data as *hard called genotype data*. Alleles
will have been set to missing based on their uncertainty. The .gcount files
for each chromosome detail the hard call genotype counts and missing rates 
for each allele (see https://www.cog-genomics.org/plink/2.0/formats#gcount)

Some additional QC has also been performed on the imputed genotype data, 
specifically I've removed duplicate variants, keeping from each pair/set
the variant with the highest imputation INFO score. There are typically a
handful of these variants in each bgen file (<100) which can cause problems
with many plink functions crashing. 

Note multi-allelic variants have been kept and are spread across multiple 
rows in the variant information and genotype data files (this is also true 
in the BGEN files). If you are working with these, you may want to create 
a new copy of the .bim files with the ID column modified to always include
the A1 and A2 alleles so that plink does not otherwise crash when 
encountering duplicate variants (i.e. because multiple rows will share the
same rsid for multi-allelic variants). 

