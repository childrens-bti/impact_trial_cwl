translate_all_phases <- function(dna_seq, coord_str, gene) {
  lapply(0:2, function(f) {
    s <- if (f > 0) subseq(dna_seq, start = f + 1) else dna_seq
    aa_seq <- AAStringSet(sub("\\*.*", "", as.character(
      translate(s, if.fuzzy.codon = "X", no.init.codon = TRUE))))
    names(aa_seq) <- glue::glue("{coord_str}_{gene}_phase{f}")
    return(aa_seq)
  })
}