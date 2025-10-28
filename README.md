# Cortical Development Gene Mining

Automated pipeline to find genes associated with cortical folding in gyrencephalic mammals (ferret, cat, dog, primates) by mining PubMed literature.

## What It Does

- Searches PubMed for cortical development papers
- Extracts gene mentions from 1,100+ abstracts
- Identifies gyrencephalic-enriched genes
- Creates 6 visualizations

## Quick Start

```bash
# Install R packages
Rscript -e "install.packages(c('rentrez', 'dplyr', 'ggplot2', 'stringr', 'tidyr', 'viridis', 'readr'))"

# Run pipeline
Rscript cortical_pipeline.R
```

**Time:** 30 minutes  
**Requirements:** R 4.0+, 8GB RAM, internet

## Results

**Files generated:**
- `cortical_papers.csv` - PubMed papers retrieved
- `top_cortical_genes.csv` - Genes ranked by frequency
- `gyrencephalic_genes.csv` - Species-enriched genes
- `plots/*.png` - 6 visualizations (300 DPI)

**Example output:**

| Gene | Papers | Gyrencephalic Score |
|------|--------|---------------------|
| TRNP1 | 76 | +60 |
| ARHGAP11B | 54 | +50 |
| PAX6 | 180 | +12 |

## Customize

Edit `cortical_pipeline.R`:

```r
# Add genes
known_genes <- c('PAX6', 'TBR2', 'YOUR_GENE')

# Add species
species_keywords$cat <- c('cat', 'feline', 'felis')

# Change search
queries <- c('YOUR_QUERY_HERE')
```

## Species Analyzed

**Gyrencephalic:** ferret, cat, dog, pig, sheep, human, macaque, marmoset, chimpanzee  
**Lissencephalic:** mouse, rat

## Citation

```bibtex
@software{cortical_mining2025,
  author = {Samira},
  title = {Cortical Development Gene Mining},
  year = {2025},
  url = {https://github.com/Samira0250/LLM_Cortical_genes_litReveiw)}
}
```

## License

MIT

