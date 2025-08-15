awk '$3=="gene"' Climon_v1_primary-annotation.gtf | awk '{print $1"\t"$4"\t"$5"\t"$10"\t"$7}'| sed 's/"//g'| sed 's/;//' > genes.bed5
awk -F'\t' '$7=="LTR"||$7=="DNA"{print $1"\t"$2"\t"$3"\t"$7"\t"$6}' TEannotation/genome.fa_rm.bed >TEs.bed5
