# STRING's own rule for combining evidence channels into one score (the rule behind `combined_score` in
# protein.links.full): per channel, remove the prior (0.041), merge the direct and transferred parts of a channel as
# independent evidence, then combine all channels the same way and add the prior back:
#   s' = (max(s, prior) - prior) / (1 - prior);  combined = (1 - prod(1 - s'_channel)) * (1 - prior) + prior
# Source: STRING FAQ "How are the scores computed?" (https://string-db.org/help/faq/) and STRING's published script
# combine_subscores.v2.py, "algorithm updated for versions v12 and up"
# (https://stringdb-downloads.org/download/combine_subscores.v2.py). Same arithmetic step for step (no homology
# correction; truncated to 1/1000), but for a chosen list of channels: leaving a channel out equals setting it to 0
# in that script.
# Applied to ALL channels it reproduces STRING's combined_score in 9606.protein.links.full.v12.5 to within 3/1000
# (74.5% exact, 99.55% within 1/1000; all 10.8 million rows, 2026-11-11).
STRING_PRIOR <- 0.041

GENOMIC_CHANNELS <- c("neighborhood", "fusion", "cooccurence")   # genomic context; homology is not a channel

#' Combine the chosen evidence channels of a protein.links.full table with STRING's rule.
#' @param fullint protein.links.full table (scores 0-1000)
#' @param channels channel names (e.g. "experiments"); a channel's "<name>_transferred" column is included only
#'   when `transferred = TRUE` and the column exists
#' @return combined score, 0-1000 integer scale like STRING's combined_score
combine_string_channels <- function(fullint, channels, transferred = FALSE) {
  away <- function(s) { s <- pmax(s / 1000, STRING_PRIOR); (s - STRING_PRIOR) / (1 - STRING_PRIOR) }
  one_minus <- rep(1, nrow(fullint))
  for (ch in channels) {
    p <- away(fullint[[ch]])
    tcol <- paste0(ch, "_transferred")
    if (transferred && tcol %in% names(fullint)) p <- 1 - (1 - p) * (1 - away(fullint[[tcol]]))
    one_minus <- one_minus * (1 - p)
  }
  floor(((1 - one_minus) * (1 - STRING_PRIOR) + STRING_PRIOR) * 1000)
}

#' Build the PPI network from STRING's protein.links.full table.
#'
#' The edge score is recomputed with STRING's own combination rule from the enabled channels only, so evidence
#' that is switched off (transferred from other species, genomic context, any channel not in `best_evidence`)
#' contributes nothing to the score, not only to the keep/drop decision. (Before 2026-11-11 the switched-off
#' channels were removed only from the keep test, and the cut-off was applied to STRING's full combined_score,
#' which still contained them.)
#'
#' @param links unused (kept for the call in generate_string_ppi.rmd)
#' @param info protein.info table, or its path/URL (first two columns: STRING id, preferred name)
#' @param fullint protein.links.full table
#' @param best_evidence direct channels to use
#' @param transferred include each channel's evidence transferred from other species
#' @param genomic include the genomic-context channels (neighborhood, fusion, cooccurence)
#' @param to_kinograte keep edges with score >= `min_score` and return head/tail/cost
#' @param min_score cut-off on the recomputed score, 0-1 scale
generate_string_ppi <- function(links, info, fullint, best_evidence = c('experiments', 'textmining', 'database', 'coexpression'),
                                transferred = F, genomic = F, to_kinograte = T, min_score = 0.5) {
  names <- if (is.character(info)) readr::read_tsv(info, show_col_types = FALSE) else info
  names <- names[, 1:2]; colnames(names) <- c("string_protein_id", "preferred_name")

  channels <- c(best_evidence, if (genomic) GENOMIC_CHANNELS)
  fullint$combined_score <- combine_string_channels(fullint, channels, transferred = transferred)

  network <- fullint %>%
    filter(combined_score > floor(STRING_PRIOR * 1000)) %>%   # at least one enabled channel has evidence
    left_join(names, by = c("protein1" = "string_protein_id")) %>%
    rename(protein1_symbol = preferred_name) %>%
    left_join(names, by = c("protein2" = "string_protein_id")) %>%
    rename(protein2_symbol = preferred_name) %>%
    select(protein1_symbol, protein2_symbol, combined_score) %>%
    dplyr::rename(head = protein1_symbol, tail = protein2_symbol)

  if (to_kinograte) {
    # keeps edges with score >= min_score, then converts to cost
    network <- network %>%
      mutate(combined_score = combined_score / 1000) %>%
      filter(combined_score >= min_score) %>%
      mutate(cost = pmax(0.01, 1 - combined_score)) %>%
      select(head, tail, cost)
  }

  return(network)
}
