library(tidyverse)
library(topGO)


#-------------- Combining Annotation Files --------------------------------

# Make a list of all of the files that will be combined 
file_list <- list.files(path = choose.dir(), # Selected Hap2_results file
                        full.names = TRUE)

# Combine all of the files (only include IDs and GO terms) and clean #query names to match given table
combined <- map_dfr(file_list, ~read_tsv(.x, comment = "##") %>% 
                      dplyr::select(`#query`, GOs) %>%
                      mutate(`#query` = gsub("^(MSTRG\\.\\d+).*$", "\\1", `#query`),  # trims MSTRG
                             `#query` = gsub("\\.t.*$", "", `#query`)))    

#------------- Filter DEGs -------------------------------------------------

# Load whole genome GO file
whole_genome <- read_csv(choose.files()) # Selected the "combined" file from above

# Load DEG list
deg_list <- read_csv(choose.files()) # Selected the table of DEGs sent to me 

# Filter whole genome GO file to just DEGs
deg_GO <- whole_genome %>%
  filter(`#query` %in% deg_list$id)

# Do some number checking 
# How many have GO terms?
deg_GO %>% filter(GOs != "-") %>% nrow()

# Which DEGs are NOT in the whole genome GO file at all
missing <- deg_list$id[!deg_list$id %in% whole_genome$`#query`]
length(missing)

# Make the final list of DEGs filtered to be only the 72 present 
deg_GO_annotated <- deg_GO %>% 
  filter(GOs != "-")

write_csv(deg_GO_annotated, "genesdesc_GO_DEG")

#------------------- Running topGO -------------------------------------------

# Create a named list of GO terms for the whole genome (background)
# Split the comma-separated GO terms into a list
gene2GO <- whole_genome %>%
  filter(GOs != "-") %>%
  mutate(GOs = strsplit(GOs, ",")) %>%
  deframe()

# Create a named vector of your DEGs (1 = DEG, 0 = not DEG)
gene_list <- as.integer(whole_genome$`#query` %in% deg_list$id)
names(gene_list) <- whole_genome$`#query`

# Create topGO object
GOdata <- new("topGOdata",
              ontology = "BP",  # Use Biological Process
              allGenes = gene_list,
              geneSel = function(x) x == 1, # Select genes where value = 1 (DEGs)
              gene2GO = gene2GO,
              annot = annFUN.gene2GO)

# ----------------------- weight01 Version -------------------------------

# Run the hypergeometric test with weight01 
results_weight01 <- runTest(GOdata, algorithm = "weight01",
                   statistic = "fisher") # Fisher's = hypergeometric test


# Get ALL GO terms that were tested 
results_weight01_full <- GenTable(GOdata, 
                         weight01_p = results_weight01,
                         topNodes = length(score(results_weight01)),
                         numChar = 1000) 

# Clean up and filter to only statistically significant results
sig_results_weight01 <- results_weight01_full %>% 
  mutate(weight01_p = as.numeric(gsub("< ", "", weight01_p))) %>% 
  filter(weight01_p < 0.05) 

# Export the topology-filtered data frame
write_csv(sig_results_weight01, "enriched_GO_terms_weight01.csv")

# ----------------------- Genes per GO Term -------------------------------
sig_genes <- names(which(gene_list == 1))

terms_of_interest <- c("GO:0032989", "GO:0008219", "GO:0030168", "GO:0009860",
                       "GO:0009741", "GO:0009911", "GO:0010483", "GO:0009740",
                       "GO:0048574", "GO:0006468", "GO:0010371", "GO:0010321", 
                       "GO:0007178")

term_gene_table <- lapply(terms_of_interest, function(term) {
  tibble(
    GO_ID   = term,
    gene_id = intersect(genesInTerm(GOdata, term)[[1]], sig_genes)
  )
}) %>%
  bind_rows() %>%
  left_join(sig_results_weight01 %>% dplyr::select(GO.ID, Term, weight01_p),
            by = c("GO_ID" = "GO.ID"))

write_csv(term_gene_table, "genes_per_GO_term.csv")

# ------------------- Combining DEG and New Table with Genes ----------------

# Import the file with GO terms and their genes 
GO_gene_table <- read_csv(file.choose()) # genes_per_GO_term

GO_gene_deg <- GO_gene_table %>%
  left_join(deg_list, by = c("gene_id" = "id")) %>%
  mutate(bias = ifelse(fc > 1, "male-biased", "female-biased"))

write_csv(GO_gene_deg, "GO_gene_deg.csv")

# Look into the biases found between significant genes
# Which genes are more expressed in females/males
bias_summary <- GO_gene_deg %>%
  group_by(GO_ID, desc, bias) %>%
  summarise(n_genes = n(), .groups = "drop")

#fc value
GO_gene_deg %>%
  arrange(desc(fc)) %>%
  dplyr::select(GO_ID, desc, gene_id, fc, qval, bias)


#------------- Make a Version that Extracts ALL GENES for ALL GO TERMS--------

sig_genes <- names(which(gene_list == 1))

# Use ALL significant terms instead of just the 13
all_sig_terms <- sig_results_weight01$GO.ID

term_gene_table_full <- lapply(all_sig_terms, function(term) {
  tibble(
    GO_ID   = term,
    gene_id = intersect(genesInTerm(GOdata, term)[[1]], sig_genes)
  )
}) %>%
  bind_rows() %>%
  left_join(sig_results_weight01 %>% dplyr::select(GO.ID, Term, weight01_p),
            by = c("GO_ID" = "GO.ID"))

GO_gene_deg_full <- term_gene_table_full %>%
  left_join(deg_list, by = c("gene_id" = "id")) %>%
  mutate(
    bias = ifelse(fc > 1, "male-biased", "female-biased"),
    desc = Term
  )

write_csv(GO_gene_deg_full, "All_terms_GO_gene_deg.csv")
