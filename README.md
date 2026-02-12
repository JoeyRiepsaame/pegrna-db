---
title: pegRNA-DB
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: docker
pinned: false
license: mit
---

# pegRNA Database

An automated database of experimentally validated prime editing guide RNAs (pegRNAs and epegRNAs) extracted from scientific literature.

## Features

- **Search & Browse**: Filter by target gene, edit type, cell type, prime editor version, and efficiency
- **Statistics Dashboard**: Visualize database contents with interactive Plotly charts
- **Export**: Download data in CSV or JSON format
- **Detailed Entry View**: See complete sequence information, experimental conditions, and source papers

## What is Prime Editing?

Prime editing is a CRISPR-based genome editing technology that uses a prime editing guide RNA (pegRNA) to direct precise edits without requiring double-strand breaks or donor DNA templates.

## Data Source

The database extracts data from peer-reviewed publications indexed in PubMed using a hybrid pipeline:
- Rule-based parsing for structured supplementary tables
- LLM-powered extraction for unstructured text
- Validation of sequences and experimental parameters

## Using the Space

### Browse the Database
- Navigate to "Search & Browse" to explore all pegRNA entries
- Use filters to find entries matching your criteria
- Click on entries to view full details including sequences and experimental conditions

### View Statistics
- See aggregate statistics about the database
- Explore distributions of edit types, efficiency, and target genes
- Compare pegRNA vs epegRNA performance

### Export Data
- Download filtered or complete datasets in CSV or JSON format
- Use exported data in your own analysis pipelines

### Submit Papers (Local Only)
- The "Submit Paper" feature requires API keys and is disabled on HF Spaces
- To use this feature, clone the repository and run locally with proper API credentials

## Local Development

```bash
git clone https://github.com/JoeyRiepsaame/pegrna-db.git
cd pegrna-db
pip install -r requirements.txt

# Create .env file with API keys (optional, only for paper submission)
# ANTHROPIC_API_KEY=your_key_here
# NCBI_EMAIL=your_email@example.com

streamlit run app.py
```

## Data Fields

Each pegRNA entry includes:
- **Sequences**: spacer, PBS, RTT, 3' extension, full sequence, nicking sgRNA
- **Target Information**: gene, genomic locus, organism
- **Edit Details**: type (substitution/insertion/deletion), description, intended mutation
- **Experimental Conditions**: prime editor version, cell type, delivery method
- **Results**: editing efficiency, product purity, indel frequency
- **Metadata**: source paper (PMID), confidence score, validation status

## Confidence Scoring

- **0.8+**: Rule-based extraction from structured tables (highest confidence)
- **0.5-0.7**: LLM-based extraction from paper text (moderate confidence)
- **Validated**: Human-reviewed entries

## Citation

If you use this database in your research, please cite:

```
pegRNA Database (2026)
https://huggingface.co/spaces/JoeyRiepsaame/pegRNA-DB
GitHub: https://github.com/JoeyRiepsaame/pegrna-db
```

## Contributing

Contributions are welcome! Please submit issues or pull requests on [GitHub](https://github.com/JoeyRiepsaame/pegrna-db).

## License

MIT License - see LICENSE file for details

## Contact

For questions or feedback, please open an issue on GitHub.
