generate_string_ppi <- function(links, info, fullint, best_evidence = c('experiments', 'textmining', 'database', 'coexpression'), 
                                transferred = F, genomic = F, to_kinograte = T ) {
  
  # neighborhood, fusion, cooccurence, homology evidence types are summarized in genomic type.
  
  # remove transferred evidence
  if (!transferred) {
    transferred_cols <- grep("_transferred$", names(fullint), value = TRUE)
    fullint[transferred_cols] <- 0
  }
  
  if (genomic) {
    # some evidence types are to be summarized in the genomic column:
    fullint <- fullint %>% 
      mutate(genomic = max(neighborhood, fusion, cooccurence, homology, na.rm = T))
    
    best_evidence <- c(best_evidence, "genomic")
  }
  
  
  # Filter out rows that have only untrusted evidence. 
  # In other words: “Keep interaction if it has support from at least one enabled evidence type.”

  fullint_best <- fullint %>%
    filter(rowSums(across(all_of(best_evidence)) != 0) > 0)
  # parse
  
  network <- fullint_best %>%
    left_join(names, by = c("protein1" = "#string_protein_id")) %>%
    rename(protein1_symbol = preferred_name) %>%
    left_join(names, by = c("protein2" = "#string_protein_id")) %>%
    rename(protein2_symbol = preferred_name) %>% 
    select(protein1_symbol, protein2_symbol, combined_score) %>%
    dplyr::rename(head = protein1_symbol, tail = protein2_symbol)

  if (to_kinograte){
    # keeps top 50% of combined score, then converts to cost
    network <- network %>%
      mutate(combined_score = combined_score / 1000) %>%
      filter(combined_score >= .5) %>%
      mutate(cost = pmax(0.01, 1-combined_score)) %>%
      select(head, tail, cost) 
  } 
  
  
  return(network)
}
  
  
  