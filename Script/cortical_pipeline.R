# Cortical Development Gene Mining Pipeline in R -  Version 3
# Includes "all(?)" mammalian species + visualizations
# By: Samira Salimiyan-- CBC class (wk11)
# Date: 10-28-2025

setwd("C:/Users/szs0394/Downloads/wk11")

# Install packages
install.packages(c("rentrez", "dplyr", "tidyr", "stringr", "readr"))



library(rentrez)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(viridis)
library(scales)

# ============================================================================
# SETUP
# ============================================================================

# Define search queries for gyrencephalic cortical development
queries <- c(
  '("cortical development" OR "cortical folding" OR "gyrification") AND (ferret OR primate OR gyrencephalic OR carnivore)',
  '"radial glia" AND (OSVZ OR "outer subventricular zone") AND mammal',
  '"basal progenitor" AND cortex AND (ferret OR primate OR cat OR dog)',
  'gyrencephalic AND "gene expression" AND cortex',
  '("outer radial glia" OR oRG) AND cortex'
)

# Known cortical development genes (expanded)
known_genes <- c(
  'PAX6', 'TBR2', 'EOMES', 'TBR1', 'CTNNB1', 'TRNP1', 'ARHGAP11B',
  'FAM107A', 'NEUROG2', 'NEUROD1', 'SOX2', 'NES', 'NESTIN',
  'ASCL1', 'HES1', 'HES5', 'NOTCH1', 'NOTCH2', 'DLL1', 'DLL3',
  'FGF2', 'FGFR1', 'FGFR2', 'FGFR3', 'EGF', 'EGFR',
  'WNT3A', 'WNT7A', 'GSK3B', 'APC', 'LEF1', 'TCF7L2',
  'FOXG1', 'EMX1', 'EMX2', 'LHX2', 'NR2F1', 'NR2F2',
  'HOPX', 'PDGFD', 'CDC42', 'CENPE', 'CENPF',
  'TP53', 'P53', 'CDK4', 'CCND1', 'CCND2',
  'PTEN', 'AKT1', 'MTOR', 'TSC1', 'TSC2',
  'ASPM', 'MCPH1', 'CDK5RAP2', 'CENPJ', 'CEP152',
  'REELIN', 'DAB1', 'LIS1', 'DCX', 'TUBA1A',
  'MAP2', 'TUBB3', 'MAPT', 'NCAM1',
  'SATB2', 'CUX1', 'CUX2', 'BCL11B', 'FEZF2',
  'TLE4', 'INSM1', 'HEY1', 'HEY2', 'ID2', 'ID4'
)

# Gene aliases
gene_aliases <- list(
  'BETA-CATENIN' = 'CTNNB1',
  'β-CATENIN' = 'CTNNB1',
  'P53' = 'TP53',
  'NESTIN' = 'NES',
  'COUP-TFI' = 'NR2F1',
  'COUP-TFII' = 'NR2F2'
)

# All mammalian species studied in cortical development
species_keywords <- list(
  # Gyrencephalic species
  ferret = c('ferret', 'mustela putorius'),
  cat = c('cat', 'feline', 'felis catus', 'felid'),
  dog = c('dog', 'canine', 'canis familiaris'),
  sheep = c('sheep', 'ovine', 'ovis aries'),
  pig = c('pig', 'porcine', 'sus scrofa', 'swine'),
  
  # Primates (gyrencephalic)
  macaque = c('macaque', 'rhesus', 'macaca mulatta', 'cynomolgus'),
  marmoset = c('marmoset', 'callithrix'),
  human = c('human', 'homo sapiens'),
  chimpanzee = c('chimpanzee', 'pan troglodytes'),
  
  # Lissencephalic species (for comparison)
  mouse = c('mouse', 'mice', 'murine', 'mus musculus'),
  rat = c('rat', 'rattus norvegicus'),
  
  # Other mammals
  rabbit = c('rabbit', 'oryctolagus'),
  hamster = c('hamster', 'mesocricetus')
)

# Create gyrencephalic groups
gyrencephalic_species <- c('ferret', 'cat', 'dog', 'sheep', 'pig', 
                           'macaque', 'marmoset', 'human', 'chimpanzee')
lissencephalic_species <- c('mouse', 'rat')

# ============================================================================
# STEP 1: FETCH PUBMED PAPERS
# ============================================================================

cat("[STEP 1] Fetching PubMed papers...\n")
cat("--------------------------------------------------------------------------------\n")

fetch_abstracts <- function(query, max_results = 200) {
  cat(sprintf("\nSearching: %s...\n", substr(query, 1, 60)))
  
  tryCatch({
    search_result <- entrez_search(
      db = "pubmed",
      term = query,
      retmax = max_results,
      use_history = TRUE
    )
    
    id_list <- search_result$ids
    cat(sprintf("  Found %d papers\n", length(id_list)))
    
    if (length(id_list) == 0) {
      return(data.frame())
    }
    
    papers <- data.frame()
    batch_size <- 50
    
    for (i in seq(1, length(id_list), batch_size)) {
      batch_ids <- id_list[i:min(i + batch_size - 1, length(id_list))]
      
      tryCatch({
        summaries <- entrez_summary(db = "pubmed", id = batch_ids)
        abstracts <- entrez_fetch(
          db = "pubmed",
          id = batch_ids,
          rettype = "abstract",
          retmode = "text"
        )
        
        for (j in seq_along(batch_ids)) {
          pmid <- batch_ids[j]
          
          tryCatch({
            summary <- summaries[[as.character(pmid)]]
            title <- summary$title
            year <- summary$pubdate
            
            abstract_text <- ""
            if (is.character(abstracts)) {
              abstract_text <- abstracts
            }
            
            papers <- rbind(papers, data.frame(
              pmid = pmid,
              title = title,
              abstract = abstract_text,
              year = substr(year, 1, 4),
              query = substr(query, 1, 50),
              stringsAsFactors = FALSE
            ))
          }, error = function(e) { })
        }
        
        Sys.sleep(0.34)
        
      }, error = function(e) {
        cat(sprintf("  Error in batch: %s\n", e$message))
      })
    }
    
    return(papers)
    
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
    return(data.frame())
  })
}

# Fetch all papers
all_papers <- data.frame()
for (query in queries) {
  papers <- fetch_abstracts(query, max_results = 200)
  all_papers <- rbind(all_papers, papers)
}

all_papers <- all_papers %>%
  distinct(pmid, .keep_all = TRUE)

cat(sprintf("\n✓ Total unique papers: %d\n", nrow(all_papers)))
write_csv(all_papers, "cortical_papers.csv")
cat("✓ Saved to: cortical_papers.csv\n")

# ============================================================================
# STEP 2: EXTRACT GENE MENTIONS
# ============================================================================

cat("\n[STEP 2] Extracting gene mentions...\n")
cat("--------------------------------------------------------------------------------\n")

extract_genes <- function(text) {
  if (is.na(text) || text == "") {
    return(character(0))
  }
  
  text_upper <- toupper(text)
  found_genes <- character(0)
  
  for (gene in known_genes) {
    pattern <- paste0("\\b", gene, "\\b")
    if (str_detect(text_upper, pattern)) {
      found_genes <- c(found_genes, gene)
    }
  }
  
  for (alias in names(gene_aliases)) {
    pattern <- paste0("\\b", alias, "\\b")
    if (str_detect(text_upper, pattern)) {
      canonical <- gene_aliases[[alias]]
      if (!(canonical %in% found_genes)) {
        found_genes <- c(found_genes, canonical)
      }
    }
  }
  
  return(unique(found_genes))
}

cat("Extracting genes from abstracts...\n")

gene_mentions <- data.frame()

for (i in 1:nrow(all_papers)) {
  row <- all_papers[i, ]
  full_text <- paste(row$title, row$abstract)
  
  genes <- extract_genes(full_text)
  
  if (length(genes) > 0) {
    for (gene in genes) {
      gene_mentions <- rbind(gene_mentions, data.frame(
        pmid = row$pmid,
        gene = gene,
        year = row$year,
        title = substr(row$title, 1, 100),
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat(sprintf("\n✓ Total gene mentions: %d\n", nrow(gene_mentions)))
cat(sprintf("✓ Unique genes found: %d\n", length(unique(gene_mentions$gene))))
write_csv(gene_mentions, "gene_mentions.csv")
cat("✓ Saved to: gene_mentions.csv\n")

# ============================================================================
# STEP 3: RANK GENES BY FREQUENCY
# ============================================================================

cat("\n[STEP 3] Ranking genes...\n")
cat("--------------------------------------------------------------------------------\n")

gene_counts <- gene_mentions %>%
  group_by(gene) %>%
  summarise(paper_count = n_distinct(pmid)) %>%
  arrange(desc(paper_count))

top_genes <- head(gene_counts, 30)

cat("\n================================================================================\n")
cat("TOP 30 GENES IN CORTICAL DEVELOPMENT LITERATURE\n")
cat("================================================================================\n\n")
cat(sprintf("%-6s %-15s %-10s %-12s\n", "Rank", "Gene", "Papers", "% of Total"))
cat("--------------------------------------------------------------------------------\n")

total_papers <- nrow(all_papers)
for (i in 1:nrow(top_genes)) {
  gene <- top_genes$gene[i]
  count <- top_genes$paper_count[i]
  percentage <- (count / total_papers) * 100
  cat(sprintf("%-6d %-15s %-10d %6.1f%%\n", i, gene, count, percentage))
}

write_csv(top_genes, "top_cortical_genes.csv")
cat("\n✓ Saved to: top_cortical_genes.csv\n")

# ============================================================================
# STEP 4: SPECIES-SPECIFIC ANALYSIS
# ============================================================================

cat("\n[STEP 4] Species-specific analysis...\n")
cat("--------------------------------------------------------------------------------\n")

detect_species <- function(text) {
  if (is.na(text) || text == "") {
    return(character(0))
  }
  
  text_lower <- tolower(text)
  found_species <- character(0)
  
  for (species in names(species_keywords)) {
    keywords <- species_keywords[[species]]
    for (keyword in keywords) {
      if (str_detect(text_lower, keyword)) {
        found_species <- c(found_species, species)
        break
      }
    }
  }
  
  return(unique(found_species))
}

cat("Detecting species mentions...\n")

species_mentions <- data.frame()

for (i in 1:nrow(all_papers)) {
  row <- all_papers[i, ]
  full_text <- paste(row$title, row$abstract)
  species <- detect_species(full_text)
  
  if (length(species) > 0) {
    for (sp in species) {
      species_mentions <- rbind(species_mentions, data.frame(
        pmid = row$pmid,
        species = sp,
        stringsAsFactors = FALSE
      ))
    }
  }
}

species_counts <- species_mentions %>%
  count(species, sort = TRUE)

cat("\n================================================================================\n")
cat("PAPERS BY SPECIES\n")
cat("================================================================================\n\n")

# Categorize species
species_counts <- species_counts %>%
  mutate(
    category = case_when(
      species %in% gyrencephalic_species ~ "Gyrencephalic",
      species %in% lissencephalic_species ~ "Lissencephalic",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(desc(n))

for (i in 1:nrow(species_counts)) {
  species <- species_counts$species[i]
  count <- species_counts$n[i]
  category <- species_counts$category[i]
  percentage <- (count / nrow(all_papers)) * 100
  cat(sprintf("%15s (%s): %4d papers (%.1f%%)\n", 
              species, category, count, percentage))
}

# ============================================================================
# STEP 5: GENE x SPECIES MATRIX
# ============================================================================

cat("\n[STEP 5] Creating gene-species associations...\n")
cat("--------------------------------------------------------------------------------\n")

# Create gene-species matrix
gene_species_matrix <- data.frame()

for (i in 1:min(20, nrow(top_genes))) {
  gene <- top_genes$gene[i]
  
  gene_pmids <- gene_mentions %>%
    filter(gene == !!gene) %>%
    pull(pmid) %>%
    unique()
  
  for (sp in names(species_keywords)) {
    count <- species_mentions %>%
      filter(pmid %in% gene_pmids, species == sp) %>%
      nrow()
    
    gene_species_matrix <- rbind(gene_species_matrix, data.frame(
      gene = gene,
      species = sp,
      count = count,
      stringsAsFactors = FALSE
    ))
  }
}

# Calculate gyrencephalic enrichment score
gyren_genes <- data.frame()

for (i in 1:min(20, nrow(top_genes))) {
  gene <- top_genes$gene[i]
  
  gene_data <- gene_species_matrix %>%
    filter(gene == !!gene)
  
  gyren_count <- gene_data %>%
    filter(species %in% gyrencephalic_species) %>%
    summarise(total = sum(count)) %>%
    pull(total)
  
  lissen_count <- gene_data %>%
    filter(species %in% lissencephalic_species) %>%
    summarise(total = sum(count)) %>%
    pull(total)
  
  gyren_score <- gyren_count - lissen_count
  
  if (gyren_score > 0) {
    gyren_genes <- rbind(gyren_genes, data.frame(
      gene = gene,
      gyrencephalic = gyren_count,
      lissencephalic = lissen_count,
      gyren_score = gyren_score,
      stringsAsFactors = FALSE
    ))
  }
}

gyren_genes <- gyren_genes %>%
  arrange(desc(gyren_score))

cat("\n================================================================================\n")
cat("GYRENCEPHALIC-ENRICHED GENES\n")
cat("================================================================================\n\n")

cat(sprintf("%-15s %-15s %-15s %-10s\n", 
            "Gene", "Gyrencephalic", "Lissencephalic", "Score"))
cat("--------------------------------------------------------------------------------\n")

for (i in 1:min(15, nrow(gyren_genes))) {
  row <- gyren_genes[i, ]
  cat(sprintf("%-15s %-15d %-15d %-10d\n",
              row$gene, row$gyrencephalic, row$lissencephalic, row$gyren_score))
}

write_csv(gyren_genes, "gyrencephalic_genes.csv")
cat("\n✓ Saved to: gyrencephalic_genes.csv\n")

# ============================================================================
# STEP 6: CREATE VISUALIZATIONS
# ============================================================================

# Create plots directory
dir.create("plots", showWarnings = FALSE)

# -------------------------
# PLOT 1: Lollipop Chart of Top 20 Genes
# -------------------------

cat("Creating lollipop chart...\n")

top_20_genes <- head(gene_counts, 20)

p1 <- ggplot(top_20_genes, aes(x = reorder(gene, paper_count), y = paper_count)) +
  geom_segment(aes(x = reorder(gene, paper_count), 
                   xend = reorder(gene, paper_count),
                   y = 0, yend = paper_count),
               color = "steelblue", size = 1.5) +
  geom_point(color = "steelblue", size = 4) +
  geom_text(aes(label = paper_count), hjust = -0.5, size = 3) +
  coord_flip() +
  labs(title = "Top 20 Genes in Cortical Development Literature",
       subtitle = paste0("Based on ", nrow(all_papers), " PubMed papers"),
       x = "Gene",
       y = "Number of Papers") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ylim(0, max(top_20_genes$paper_count) * 1.1)

ggsave("plots/01_lollipop_top_genes.png", p1, width = 10, height = 8, dpi = 300)
cat("  ✓ Saved: plots/01_lollipop_top_genes.png\n")

# -------------------------
# PLOT 2: Species Distribution Pie Chart
# -------------------------

cat("Creating species pie chart...\n")

species_for_pie <- species_counts %>%
  head(10) %>%
  mutate(percentage = n / sum(n) * 100)

p2 <- ggplot(species_for_pie, aes(x = "", y = n, fill = species)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(species, "\n", round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 3) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(title = "Distribution of Papers by Species",
       subtitle = "Top 10 Most Studied Species") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "none"
  )

ggsave("plots/02_pie_species_distribution.png", p2, width = 10, height = 8, dpi = 300)
cat("  ✓ Saved: plots/02_pie_species_distribution.png\n")

# -------------------------
# PLOT 3: Bubble Plot - Gene x Species
# -------------------------

# Filter for top 15 genes and species with data
top_15_genes <- head(top_genes$gene, 15)
bubble_data <- gene_species_matrix %>%
  filter(gene %in% top_15_genes, count > 0) %>%
  left_join(species_counts %>% select(species, category), by = "species")

p3 <- ggplot(bubble_data, aes(x = species, y = gene, size = count, color = category)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(2, 15), name = "Paper Count") +
  scale_color_manual(values = c("Gyrencephalic" = "#E74C3C", 
                                "Lissencephalic" = "#3498DB",
                                "Other" = "#95A5A6"),
                     name = "Brain Type") +
  labs(title = "Gene-Species Association Matrix",
       subtitle = "Top 15 genes across all species (bubble size = paper count)",
       x = "Species",
       y = "Gene") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    legend.position = "right"
  )

ggsave("plots/03_bubble_gene_species.png", p3, width = 14, height = 10, dpi = 300)
cat("  ✓ Saved: plots/03_bubble_gene_species.png\n")

# -------------------------
# PLOT 4: Heatmap - Gene x Species
# -------------------------

heatmap_data <- gene_species_matrix %>%
  filter(gene %in% top_15_genes) %>%
  pivot_wider(names_from = species, values_from = count, values_fill = 0) %>%
  pivot_longer(-gene, names_to = "species", values_to = "count")

p4 <- ggplot(heatmap_data, aes(x = species, y = gene, fill = count)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = ifelse(count > 0, count, "")), 
            color = "white", size = 3, fontface = "bold") +
  scale_fill_viridis(name = "Papers", option = "plasma") +
  labs(title = "Gene Expression Studies Across Species",
       subtitle = "Heatmap showing number of papers for each gene-species combination",
       x = "Species",
       y = "Gene") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave("plots/04_heatmap_gene_species.png", p4, width = 14, height = 10, dpi = 300)
cat("  ✓ Saved: plots/04_heatmap_gene_species.png\n")

# -------------------------
# PLOT 5: Gyrencephalic Enrichment Bar Chart
# -------------------------

gyren_top15 <- head(gyren_genes, 15)

# Reshape for stacked bar
gyren_long <- gyren_top15 %>%
  select(gene, gyrencephalic, lissencephalic) %>%
  pivot_longer(-gene, names_to = "brain_type", values_to = "count")

p5 <- ggplot(gyren_long, aes(x = reorder(gene, count), y = count, fill = brain_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("gyrencephalic" = "#E74C3C", 
                               "lissencephalic" = "#3498DB"),
                    labels = c("Gyrencephalic", "Lissencephalic"),
                    name = "Brain Type") +
  coord_flip() +
  labs(title = "Gyrencephalic-Enriched Genes",
       subtitle = "Comparison of mentions in gyrencephalic vs lissencephalic species",
       x = "Gene",
       y = "Number of Papers") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

ggsave("plots/05_gyren_enrichment.png", p5, width = 10, height = 8, dpi = 300)
cat("  ✓ Saved: plots/05_gyren_enrichment.png\n")

# -------------------------
# PLOT 6: Species Comparison - Gyrencephalic vs Lissencephalic
# -------------------------

species_comparison <- species_counts %>%
  mutate(
    brain_type = ifelse(category == "Gyrencephalic", "Gyrencephalic", "Lissencephalic")
  ) %>%
  filter(brain_type %in% c("Gyrencephalic", "Lissencephalic"))

p6 <- ggplot(species_comparison, aes(x = reorder(species, n), y = n, fill = brain_type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  scale_fill_manual(values = c("Gyrencephalic" = "#E74C3C", 
                               "Lissencephalic" = "#3498DB"),
                    name = "Brain Type") +
  coord_flip() +
  labs(title = "Cortical Development Studies by Species",
       subtitle = "Comparing gyrencephalic and lissencephalic mammals",
       x = "Species",
       y = "Number of Papers") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  ) +
  ylim(0, max(species_comparison$n) * 1.1)

ggsave("plots/06_species_comparison.png", p6, width = 10, height = 8, dpi = 300)
cat("  ✓ Saved: plots/06_species_comparison.png\n")



cat("PIPELINE COMPLETED SUCCESSFULLY\n")
