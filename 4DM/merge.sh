head -1 CpG/MA_vs_NB.dmr >merged.dmr
cat  C*/*.dmr | grep -v "start" >> merged.dmr
