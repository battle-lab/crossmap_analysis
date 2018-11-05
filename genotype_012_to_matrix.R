library(argparser)
library(data.table)

args <- arg_parser('program')

args <- add_argument(args, '-geno012',
                     help='genotype 012 file',
                     default='/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/genotype_process/chr22_maf_0.05_repeat_masked.012')
args <- add_argument(args, '-vcf',
                     help='vcf file where 012 file was generated from',
                     default='/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/genotype_process/chr22_maf_0.05_repeat_masked.vcf')
args <- add_argument(args, '-o',
                     help='output file',
                     default='results/geno.txt')

argv = parse_args(args)
genotype_012_fn = argv$geno012
vcf_fn = argv$vcf
out_fn = argv$o

### read util
read_df <- function(fn, sep = '\t',  header = T, quote = "", row.names=T, stringsAsFactors = F, check.names = F, lessColInHeader=F, skip = 0){
  if(header==T && lessColInHeader==T){
    header_line = readLines(fn, n = 1)
    headers = strsplit(header_line, split = sep)[[1]]
    skip = skip + 1
  }
  
  data_df = fread(fn, 
                  sep = sep,
                  header = header & !lessColInHeader, 
                  skip = skip,
                  quote = quote,    # not available in old package
                  stringsAsFactors = stringsAsFactors, 
                  check.names = check.names, 
                  data.table = FALSE)
  if (row.names == TRUE){
    rownames(data_df) = data_df[,1]
    data_df = data_df[,-1, drop=F]
  }
  if(header==T && lessColInHeader==T){
    colnames(data_df) = headers
  }
  return(data_df)
}

write_df <- function(x, file, sep = "\t", quote = F, row.names = T, col.names = NA){
  write.table(x = x, file = file, sep = sep, quote = quote, row.names = row.names, col.names = col.names)
}

read_vcf_variant_info <- function(vcf_fn){
  cmd = sprintf("grep -v '^#' '%s' | cut -f1-3", vcf_fn)
  lines = system(cmd, intern=T)
  infos = t(sapply(strsplit(lines, split = '\t'), function(x) x))
  colnames(infos) = c('Chr','Pos','VariantID')
  return(infos)
}


### read expression
geno_df = read_df(fn = genotype_012_fn, sep = '\t', header = F, quote = "", row.names = T, stringsAsFactors = F, lessColInHeader = F)
indiv_df = read_df(paste0(genotype_012_fn, '.indv'), header = F, row.names = F, stringsAsFactors = F)
pos_df = read_df(paste0(genotype_012_fn, '.pos'), header = F, row.names = F, stringsAsFactors = F)
vcf_infos = read_vcf_variant_info(vcf_fn)

stopifnot(nrow(geno_df)==nrow(indiv_df))
stopifnot(ncol(geno_df)==nrow(pos_df))
stopifnot(nrow(vcf_infos)==nrow(pos_df))

### get ids from variant info, and double check unique mapping
pos_id_df = cbind(pos_df, vcf_infos)
stopifnot(all(pos_id_df$V1 == pos_id_df$Chr))
stopifnot(all(pos_id_df$V2 == pos_id_df$Pos))

### set genotype id and indiv id
rownames(geno_df) = indiv_df[,1]
colnames(geno_df) = pos_id_df[,'VariantID']

### save 
geno_df = t(geno_df)
write_df(geno_df, file = out_fn)
