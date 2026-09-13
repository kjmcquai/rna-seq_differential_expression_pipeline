# Functional Genomics & GO Enrichment Pipeline (`topGO`)

An end-to-end R workflow for transcriptomic annotation, Gene Ontology (GO) over-representation analysis, and fold-change expression mapping. This pipeline parses multi-file annotation sets, builds custom Biological Process (BP) ontologies using `topGO`, applies graph-topology weighted Fisher's tests, and extracts gene-level directional expression bias.

## Computational Workflow

1. **Multi-File Annotation Parsing & Cleaning:**
   * Consolidates raw annotation files across directory listings using `map_dfr`.
   * Standardizes transcript identifiers via regex string operations (`MSTRG` trimming and suffix stripping).

2. **Foreground/Background Target Filtering:**
   * Constructs whole-genome background GO maps (`gene2GO`) and maps Differentially Expressed Gene (DEG) targets.
   * Evaluates dataset missingness and isolates annotated gene sets for statistical testing.

3. **Topology-Weighted Over-Representation Analysis (`topGO`):**
   * Initializes custom `topGOdata` objects focused on Biological Process (BP) terms.
   * Utilizes the `weight01` graph algorithm to account for GO DAG hierarchy dependencies and minimize parent-child redundancy.
   * Executes hypergeometric/Fisher’s Exact Tests to extract statistically significant nodes ($p < 0.05$).

4. **Gene-to-Term Mapping & Phenotypic Bias Layering:**
   * Queries graph structures (`genesInTerm`) to extract specific gene IDs assigned to enriched terms.
   * Intersects output tables with differential fold-change values to classify directional expression bias (e.g., male- vs. female-biased).
   * Exports comprehensive CSV summary reports across both targeted biological pathways and all significant ontology nodes.

## Tech Stack & Packages
* **Language:** R
* **Core Libraries:** `topGO`, `tidyverse` (`dplyr`, `stringr`, `purrr`, `readr`, `tibble`)

---
*Developed as part of research in the Guerrero Computational Genetics Lab at North Carolina State University.*
