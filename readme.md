# Download STRING

Download the latest version of STRING, filter for the evidence types of your choice.

The edge score is recomputed from the enabled channels only, with STRING's own combination rule (`R/helper.R`,
`combine_string_channels()`). Switched-off evidence (transferred from other species, genomic context, channels
not in `best_evidence`) therefore does not count in the score, and the `min_score` cut-off (default 0.5) applies to
that score. The rule is STRING's own: FAQ [How are the scores computed?](https://string-db.org/help/faq/) and
STRING's script [combine_subscores.v2.py](https://stringdb-downloads.org/download/combine_subscores.v2.py) (v12 and
up). The R function does the same arithmetic for a chosen subset of channels. With all channels enabled it
reproduces STRING's `combined_score` to within 3/1000.
Builds before 2026-11-11 (e.g. `ppi_network_human_filtered_v12.5.rda`, networkGen's `string_12.5`) filtered on
STRING's full `combined_score`, which still contained transferred and genomic evidence.

| Evidence Channel | Type | Reliability | Key Strength |
| --- | --- | --- | --- |
| Experiments | Direct | Very High | laboratory experiments that were conducted with the expressed goal to uncover protein–protein association evidence. BioGRID, IMEx consortium. |
| Curated Databases | Direct | Very High | well-established knowledge about protein complexes, pathways and other functional connections from dedicated knowledge resources: KEGG, Reactome, MetaCyc, EBI Complex Portal, and Gene Ontology Complexes |
| Co-expression | Functional | Indirect but useful | By comparing expression profiles (both transcript and protein abundances across various conditions), it identifies gene pairs with similar expression patterns, suggesting functional linkage. |
| Text mining | Literature | Indirect but useful | co-mentions of protein names in scientific literature |
| Gene Fusion | Genomic context | High (Prokaryote) / Low (Eukaryote) |  |
| Gene Neighborhood | Genomic context | High (Prok) / Low (Euk) | Evidence for operons and co-regulation (mainly in bacteria). |
| Gene Co-occurrence | Genomic context | High (Prok) / Low (Euk) | Evidence for participation in the same pathway/system. |

source: \[STRING publication\](<https://academic.oup.com/nar/article/53/D1/D730/7903368>)