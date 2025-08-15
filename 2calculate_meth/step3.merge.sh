Rscript ../11.software/batmeth2merge.R -p methBins.txt -s 1 -e 3
Rscript ../11.software/batmeth2merge.R -p mCcatero.txt -s 1 -e 3
Rscript ../11.software/batmeth2merge.R -p Region.CG.txt -s 1 -e 3
Rscript ../11.software/batmeth2merge.R -p Region.CHG.txt -s 1 -e 3
Rscript ../11.software/batmeth2merge.R -p Region.CHH.txt -s 1 -e 3
cat merged.Region.*.txt >merged.Region.C.txt
