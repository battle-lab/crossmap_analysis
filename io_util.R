library(data.table)

read_df <- function(fn, sep = '\t',  header = T, quote = "", row.names=T, stringsAsFactors = F, check.names = F, lessColInHeader=F, skip = 0){
  if(header==T && lessColInHeader==T){
    header_line = readLines(fn, n = skip+1)[skip+1]
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

read_gsea_db <- function(fn){
  gsea_lines = readLines(fn)
  parssed_lines = strsplit(gsea_lines, split = '\t')
  gs_titles = sapply(parssed_lines, function(x) x[1])
  genesets = lapply(parssed_lines, function(x) x[-(1:2)])
  names(genesets) = gs_titles
  return(genesets)
}
