#step1.prepare_bed.sh
awk '$3=="gene"' Climon_v1_primary-annotation.gtf | awk '{print $1"\t"$4"\t"$5"\t"$10"\t"$7}'| sed 's/"//g'| sed 's/;//' > genes.bed5
awk -F'\t' '$7=="LTR"||$7=="DNA"{print $1"\t"$2"\t"$3"\t"$7"\t"$6}' ../15.TEannotation/genome.fa_rm.bed >TEs.bed5
#step2.run_methyGff_genes.sh
nohup BatMeth2 methyGff --coverage 4 -nC 5 --genome Climon_v1_primary-scaffolds.fa -b5 genes.bed5 --body --promoter --GENE --distance 2000 --step 0.01 --methratio 1.methratio.txt --out 1.genes &
#step2.run_methyGff_TEs.sh
nohup BatMeth2 methyGff --coverage 4 -nC 5 --genome Climon_v1_primary-scaffolds.fa -b5 TEs.bed5 --body --promoter --GENE --distance 2000 --step 0.01 --methratio 1.methratio.txt --out NB1.TEs &
Rscript ../11.software/batmeth2merge.R -p genes.AverMethylevel.txt -s 1 -e 3
Rscript ../11.software/batmeth2merge.R -p TEs.AverMethylevel.txt -s 1 -e 3
Rscript ../11.software/batmeth2merge.R -p genes.Promoter.c.txt -s 1 -e 3
