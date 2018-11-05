parse_delimitted_param <- function(param_str, delim=',', rm.empty=T){
  stopifnot(class(param_str)=='character' && length(param_str)==1)
  parts = strsplit(param_str, split = delim)[[1]]
  if(rm.empty==T){
    is_valid_parts = sapply(parts, function(s) nchar(s)>0)
    parts = parts[is_valid_parts]
  }
  return(parts)
}
