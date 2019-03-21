library(argparser)
library(data.table)
source('io_util.R')

args <- arg_parser('program')
args <- add_argument(args, "-eqtl", 
                     help="eqtl file", 
                     default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/gtex_v7/eqtl_crossmap/trans_eqtl_cross_chr/combined_cross_chr_trans_eqtl_fdr_0.05.all.unique.crossmap.txt")
args <- add_argument(args, '-segdup',
                     help='segdup annotation file',
                     default='/work-zfs/abattle4/lab_data/hg19_tracks/segdup/GRCh37GenomicSuperDup.tab')
args <- add_argument(args, '-o',
                     help='output file -- eqtls where the snp is in segdup',
                     default='results/segdup_eqtl.txt')


argv <- parse_args(args)
eqtl_fn = argv$eqtl
segdup_annot_fn = argv$segdup
out_fn <- argv$o

### read files
eqtl_df = read_df(eqtl_fn, header = T, row.names = F)
segdup_df = read_df(segdup_annot_fn, header = T, row.names = F)

get_eqtls_in_segdup <- function(eqtl_df, segdup_df){
  # split segdups by chr
  chromosomes = unique(segdup_df$chrom)
  segdup_by_chr = lapply(chromosomes, function(chr) return(segdup_df[segdup_df$chrom == chr, c('chrom', 'chromStart', 'chromEnd') ,drop=F]))
  names(segdup_by_chr) = chromosomes
  
  is_snp_in_segdup <- function(snp_chr, snp_pos){
    if(!snp_chr %in% chromosomes)
      return(FALSE)
    chr_segdup = segdup_by_chr[[snp_chr]]
    is_in_segdup = sum(chr_segdup$chromStart <= snp_pos & chr_segdup$chromEnd >= snp_pos) > 0
    return(is_in_segdup)
  }
  
  is_eqtl_in_segdup = mapply(is_snp_in_segdup, snp_chr = eqtl_df$snps_chr, snp_pos = eqtl_df$snps_pos)
  return(eqtl_df[is_eqtl_in_segdup, ,drop=F])
}

eqtl_in_segdup_df = get_eqtls_in_segdup(eqtl_df, segdup_df)
write_df(eqtl_in_segdup_df, file = out_fn, row.names = F, col.names = T)

