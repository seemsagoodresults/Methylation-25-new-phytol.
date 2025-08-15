#!/usr/bin/env Rscript

###############################################
# parse parameter 
###############################################
library(argparser, quietly=TRUE)

p <- arg_parser("merge BatMath2 result")
p <- add_argument(p, "--indir", help="input dir", type="character", default = "./")
p <- add_argument(p, "--postfix", help="input file postfix", type="character")
p <- add_argument(p, "--sample_name_start", short = '-s', help="start pos for sample name", type="numeric")
p <- add_argument(p, "--sample_name_end", short = '-e', help="end pos for sample name", type="numeric")
p <- add_argument(p, "--header", short = '-c', help="input file contain header ?", type="logical", flag = T)
argv <- parse_args(p)

library(dplyr, quietly = T)
library(stringr, quietly = T)

###############################################
# merge files 
###############################################
infiles = list.files(path = argv$indir, pattern = argv$postfix, full.names = T)

i = 1
f = infiles[1]
s = str_sub(basename(f), argv$sample_name_start, argv$sample_name_end)
merged_tbl <- read.table(file = f, sep = "\t", header = argv$header)
merged_tbl$sample = s

for (i in 2:length(infiles)){
  f = infiles[i]
  s = str_sub(basename(f), argv$sample_name_start, argv$sample_name_end)
  t <- read.table(file = f, sep = "\t", header = argv$header)
  t$sample = s
  merged_tbl <- bind_rows(merged_tbl, t)
}

output = str_c(argv$indir, 'merged.', argv$postfix)
output
write.table(merged_tbl, file = output, quote = F, col.names = argv$header, row.names = F, sep = "\t")

