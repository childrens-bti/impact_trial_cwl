get_seq <- function(chr, start, end, strand) {
  if (is.na(start) | is.na(end)) return(DNAString(""))
  if (!chr %in% names(genome)) return(DNAString(""))
  start <- as.integer(start); end <- as.integer(end)
  chr_len <- length(genome[[chr]])
  if (is.na(start) | is.na(end) | start < 1 | end > chr_len | start > end) return(DNAString(""))
  seq <- subseq(genome[[chr]], start = start, end = end)
  if (strand == "-") seq <- reverseComplement(seq)
  return(seq)
}