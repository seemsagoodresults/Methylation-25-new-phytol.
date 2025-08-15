#!/usr/bin/env Rscript

###############################################
# parse parameter 
###############################################
library(argparser, quietly=TRUE)

p <- arg_parser("run dm analysis use methylKit")
p <- add_argument(p, "--samples_file", help="tab-delimited text file indicating biological replicate relationships.", type="character")
p <- add_argument(p, "--contrasts", help="file (tab-delimited) containing the pairs of sample comparisons to perform.", type="character")
p <- add_argument(p, "--fmt", help="input file format, calmeth | methylKit", type="character", default = "methylKit")

p <- add_argument(p, "--dm", help="dmc|dmr|dmp|all", type="character", default = "dmr")
p <- add_argument(p, "--gff", help="gff or gtf, required for dmp", type="character")
p <- add_argument(p, "--feature", help="gene|transcript|mRNA", type="character", default = "gene")
p <- add_argument(p, "--upstream", help="upstream TSS", type="numeric", default = 2000)
p <- add_argument(p, "--downstream", help="downstream TSS", type="numeric", default = 200)

p <- add_argument(p, "--context", help="CpG | CHG | CHH | all", type="character", default = "CpG")
p <- add_argument(p, "--window", help="window size", type="numeric", default = 1000)
p <- add_argument(p, "--step", help="step size", type="numeric", default = 200)
p <- add_argument(p, "--mincov", help="min cov", type="numeric", default = 4)
p <- add_argument(p, "--threads", help="number of threads", type="numeric", default = 4)
p <- add_argument(p, "--outdir", help="output dir", type="character", default = "./")

argv <- parse_args(p)
stopifnot(argv$context %in% c('CpG', 'CHH', 'CHG', 'all'))

# ###############################################
# # load libraries
# ###############################################
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(methylKit))
suppressMessages(library(genomation))

dir.create(argv$outdir, recursive = TRUE)

# ###############################################
# read samples file
# ###############################################

samples_info = read.table(file = argv$samples_file)
contrasts = read.table(file = argv$contrasts)

samples_info = dplyr::filter(samples_info, V1 %in% unique(c(contrasts$V1, contrasts$V2)))

context = argv$context
if (context == 'all') {
  context = c('CpG', 'CHG', 'CHH')
}

dm = argv$dm

if (dm == 'all') {
  dm = c('dmc', 'dmr', 'dmp')
}

samples_info = dplyr::filter(samples_info, V3 %in% context)

# ###############################################
# transform batmeth2 result to methlyKit
# ###############################################
for (i in 1:nrow(samples_info)) {
  fmt = argv$fmt
  
  if (fmt == 'calmeth') {# BatMeth2  calmeth result transform
    f = samples_info[i,4]
    f1 = str_c(argv$outdir, '/', basename(f))
    t = read.table(file = f, header = F)
    
    t1 = mutate(t, chrBase = str_c(V1, V2, sep = "."), chr = V1, base = V2,
                strand = if_else(V3 == '+', 'F', 'R'),
                coverage = V6, freqC = 100*V5/V6, freqT = 100-(100*V5/V6)) %>%
      dplyr::select(chrBase, chr, base, strand, coverage, freqC, freqT)
    write.table(t1, file = f1, quote = F, row.names = F, sep = '\t')
    samples_info[i,5] = f1
  }else if (fmt == 'methlyKit') {
    samples_info[i,5] = samples_info[i,4]
  }else{
    stop("fmt must be one of calmeth, methyGff, methlyKit")
  }
}

###############################################
# DM analysis use methylKit
###############################################
for (i in 1:nrow(contrasts)){
  for (c in context) {
    log = str_c("run de analysis for", contrasts[i,1], "vs", contrasts[i,2], ", context:", c, sep = " ")
    print(log)
    samples_info_i = dplyr::filter(samples_info, V1 %in% contrasts[i,], V3 == c) %>%
      mutate(treatment = if_else(V1 == contrasts[i, 1], 1, 0))
    location = samples_info_i$V5
    sample.id = samples_info_i$V2
    treatment = samples_info_i$treatment
    
    mrl = methRead(
      location = as.list(location),
      sample.id = as.list(sample.id),
      assembly = "mygenome",
      context = c,
      treatment = treatment,
      mincov = argv$mincov)
    
    # filter
    mrl_base=filterByCoverage(mrl,lo.count=NULL,lo.perc=NULL,
                         hi.count=NULL,hi.perc=99.9)
    # DMC
    if ('dmc' %in% dm) {
      mb_base=methylKit::unite(mrl_base)
      dm_base = calculateDiffMeth(mb_base, mc.cores = argv$threads)
      dm_base_tbl = as_tibble(dm_base) 
      dm_base_tbl = mutate(dm_base_tbl, sampleA = contrasts[i,1], sampleB = contrasts[i,2], context = c)
      write.table(dm_base_tbl,
                  file = str_c(argv$outdir, '/', contrasts[i,1], '_vs_', contrasts[i,2], '.dmc'),
                  sep = "\t", quote = F, row.names = F)    
    }
    
    # DMR
    if ('dmr' %in% dm) {
      mrl_region = tileMethylCounts(mrl_base,win.size=argv$window,step.size=argv$step,cov.bases = argv$mincov)
      mb_region=methylKit::unite(mrl_region)
      dm_region = calculateDiffMeth(mb_region, mc.cores = argv$threads)
      dm_region_tbl = as_tibble(dm_region) 
      dm_region_tbl = mutate(dm_region_tbl, sampleA = contrasts[i,1], sampleB = contrasts[i,2], context = c)
      write.table(dm_region_tbl,
                  file = str_c(argv$outdir, '/', contrasts[i,1], '_vs_', contrasts[i,2], '.dmr'),
                  sep = "\t", quote = F, row.names = F)      
    }
    
    # DMP
    if ('dmp' %in% dm) {
      feature_gr = gffToGRanges(argv$gff, filter = argv$feature)
      promoter_gr = promoters(feature_gr, upstream = argv$upstream, downstream = argv$downstream, use.names = TRUE)
      
      mb_base=methylKit::unite(mrl_base)
      mb_promoter = regionCounts(mb_base, regions = promoter_gr, strand.aware = FALSE, cov.bases = argv$mincov)
      
      dm_promoter = calculateDiffMeth(mb_promoter, mc.cores = argv$threads)
      dm_promoter_tbl = as_tibble(dm_promoter) 
      dm_promoter_tbl = mutate(dm_promoter_tbl, sampleA = contrasts[i,1], sampleB = contrasts[i,2], context = c)
      
      dm_promoter_tbl = left_join(dm_promoter_tbl,
        as_tibble(promoter_gr) %>%
          dplyr::select(chr = seqnames, start, ID),
        by = c("chr", "start")
        )
      
      write.table(dm_promoter_tbl,
                  file = str_c(argv$outdir, '/', contrasts[i,1], '_vs_', contrasts[i,2], '.dmp'),
                  sep = "\t", quote = F, row.names = F)   
    }
    
    
    
  }
}
