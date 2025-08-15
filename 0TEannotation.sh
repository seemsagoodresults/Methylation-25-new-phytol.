BuildDatabase -name nm -engine ncbi /home/wangyilei2020/All_genome_data/lemon_genome/Climon_v1_primary-scaffolds.fa
nohup RepeatModeler -database nm -engine ncbi -pa 20 &> nm.out&
nohup RepeatMasker -lib ../02repeatModeler/nm-families.fa -e ncbi -pa 10 -dir . /home/wangyilei2020/All_genome_data/lemon_genome/Climon_v1_primary-scaffolds.fa&
