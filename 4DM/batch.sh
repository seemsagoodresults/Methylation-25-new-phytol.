# base
awk '{print $1"\t"$2"\tCpG\t"$2"_loci.CG.txt"}' ../14.data/wgbs/reads_info.txt >samples_info.txt
awk '{print $1"\t"$2"\tCHG\t"$2"_loci.CHG.txt"}' ../14.data/wgbs/reads_info.txt >>samples_info.txt
awk '{print $1"\t"$2"\tCHH\t"$2"_loci.CHH.txt"}' ../14.data/wgbs/reads_info.txt >>samples_info.txt