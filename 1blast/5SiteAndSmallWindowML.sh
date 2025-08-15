#step1.indexgenome.sh
BatMeth2 index -g ../12.ref/genome.fa
#step2.run_calmeth.sh
nohup BatMeth2 calmeth -Q 20 --remove_dup --coverage 4 -nC 5 --Regions 200 --step 50000 --genome /Climon_v1_primary-scaffolds.fa --binput /1_clean_1_bismark_bt2_pe.deduplicated.sorted.bam --methratio 1 1>1.log 2>&1 &
#step3.merge.sh
Rscript ../11.software/batmeth2merge.R -p methBins.txt -s 1 -e 3
Rscript ../11.software/batmeth2merge.R -p mCcatero.txt -s 1 -e 3
Rscript ../11.software/batmeth2merge.R -p Region.CG.txt -s 1 -e 3
Rscript ../11.software/batmeth2merge.R -p Region.CHG.txt -s 1 -e 3
Rscript ../11.software/batmeth2merge.R -p Region.CHH.txt -s 1 -e 3
cat merged.Region.*.txt >merged.Region.C.txt