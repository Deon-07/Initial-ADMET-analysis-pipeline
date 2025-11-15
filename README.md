# ADMET Analysis Pipeline

A high-performance, memory-efficient computational pipeline for analyzing ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity) properties of chemical compounds. Developed by **Dip Kumar Ghosh** for drug discovery and compound prioritization.

## 🚀 Features

- **Comprehensive ADMET Scoring**: Multi-parameter evaluation including solubility, lipophilicity, absorption, bioavailability, CYP metabolism, and toxicity
- **Memory-Efficient Processing**: Chunked processing and streaming for large datasets
- **Parallel Computing**: Multi-processor support for accelerated scoring
- **Scientific Validation**: Built-in data quality checks and validation against experimental data
- **Flexible Configuration**: Customizable scoring weights and thresholds via JSON
- **Multiple Output Formats**: CSV, Parquet, and professional visualizations
- **Automatic Column Detection**: Smart pattern matching for diverse input formats

## 📋 Requirements

```bash
pip install pandas numpy matplotlib scikit-learn seaborn
# Optional: for Parquet output
pip install pyarrow
```

## 🛠 Installation

```bash
git clone https://github.com/DipKumarGhosh/admet-analysis-pipeline.git
cd admet-analysis-pipeline
```

## 💻 Usage

### Basic Usage
```bash
python admet_analysis.py compounds.csv
```

### Advanced Usage
```bash
python admet_analysis.py compounds.csv \
  --output ./results \
  --chunk_size 10000 \
  --workers 4 \
  --top_n 50 \
  --parquet
```

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `input_file` | Input CSV/TSV file with compound data | Required |
| `-o, --output` | Output directory | `./admet_results` |
| `-n, --top_n` | Number of top compounds to highlight | 20 |
| `--sep` | Input file separator | `\t` (tab) |
| `--config` | Custom JSON configuration file | None |
| `--chunk_size` | Process data in chunks (0 = single shot) | 0 |
| `--workers` | Number of parallel workers (0 = disabled) | 0 |
| `--no_stream` | Disable disk streaming (keep in memory) | False |
| `--parquet` | Write Parquet output (requires pyarrow) | False |

## 📊 Input Format

The pipeline expects a CSV/TSV file with compound data. Required columns:
- `cid`: Compound identifier
- `smiles`: SMILES representation

### Supported ADMET Properties (Auto-detected)
- **Solubility**: LogS, WaterSol
- **Lipophilicity**: LogP, ALogP, MLogP  
- **Absorption**: HIA, HumanIntestinalAbsorption
- **Bioavailability**: F20%, F30%, Bioavailability
- **Permeability**: Caco2, MDCK, PAMPA
- **Toxicity**: hERG, DILI, Ames, Carcinogenicity
- **CYP Metabolism**: CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4
- **Structural Alerts**: PAINS, NR-AR, SR-ARE, etc.

## ⚙️ Configuration

Create a JSON configuration file to customize scoring:

```json
{
  "adme_weights": {
    "solubility": 2.0,
    "lipophilicity": 1.5,
    "absorption": 2.0,
    "bioavailability": 2.0,
    "permeability": 1.5,
    "cyp_clean": 3.0,
    "bbb": 1.0,
    "rule_based": 1.0
  },
  "toxicity_weights": {
    "severe": 5.0,
    "moderate": 3.0,
    "mild": 1.0,
    "pains": 2.0,
    "environmental": 1.0
  },
  "thresholds": {
    "logS_pass": -4.0,
    "logS_warning": -6.0,
    "logP_optimal_min": 1.0,
    "logP_optimal_max": 5.0
  },
  "cyp_weights": {
    "CYP3A4": 3.0,
    "CYP2D6": 2.0,
    "CYP2C19": 1.5,
    "CYP2C9": 1.0,
    "CYP1A2": 1.0
  }
}
```

## 📈 Output

### Generated Files
- `admet_top_compounds_summary.csv`: Ranked list with scores
- `admet_comprehensive_results.csv`: Full analysis with all properties
- `admet_top20_detailed.csv`: Detailed analysis of top compounds
- `admet_scored_stream.csv`: Intermediate streaming output
- `admet_comprehensive_results.parquet`: Parquet format (if enabled)
- `adme_vs_toxicity.png`: Scatter plot visualization
- `final_score_distribution.png`: Score distribution histogram

### Score Components
- **ADME Score**: Drug-likeness and pharmacokinetic properties (0-100)
- **Toxicity Penalty**: Safety concerns and structural alerts
- **Final Score**: ADME Score - Toxicity Penalty
- **Rank**: Overall compound prioritization

## 🧪 Scientific Methodology

### ADME Scoring
- **Solubility**: LogS thresholds with Excellent/Acceptable/Poor grading
- **Lipophilicity**: Optimal LogP range (1-5) with penalty for extremes
- **Absorption**: Human Intestinal Absorption (HIA) classification
- **Bioavailability**: Oral bioavailability flags (F20%, F30%)
- **CYP Profile**: Weighted inhibition scoring (CYP3A4 > CYP2D6 > others)
- **Rule-based**: Compliance with Lipinski, Veber, Pfizer, GSK rules

### Toxicity Assessment
- **Severe Toxicities**: hERG, DILI, Ames, Carcinogenicity (hard failures)
- **Moderate Concerns**: Skin sensitization, eye irritation
- **Mild Alerts**: Nuclear receptor and stress response pathways
- **Structural Alerts**: PAINS filters and toxicophores

## 🔬 Validation

The pipeline includes validation against experimental data:

```python
# Example validation usage
experimental_data = pd.read_csv('experimental_results.csv')
metrics = analyzer.validate_admet_predictions(
    experimental_data,
    mapping={'LogP_pred': 'LogP_exp', 'LogS_pred': 'LogS_exp'}
)
```

## 🚀 Performance

### Memory Optimization
- **Chunked Processing**: Handle datasets larger than memory
- **Streaming to Disk**: Avoid memory exhaustion
- **Category Data Types**: Reduce memory footprint
- **Efficient Garbage Collection**: Manual memory management

### Parallel Processing
- **Multi-core Scoring**: Scale with available CPUs
- **Process-based Parallelism**: Avoid GIL limitations
- **Configurable Workers**: Balance performance and memory

## 🐛 Troubleshooting

### Common Issues
1. **Memory Errors**: Use `--chunk_size` and `--no_stream` for large datasets
2. **Column Detection**: Ensure column names follow supported patterns
3. **Separator Issues**: Specify correct `--sep` parameter
4. **Performance**: Enable `--workers` for multi-core systems

### Logging
Detailed logs are available in the console and can be redirected to files for debugging.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for:
- New ADMET property integrations
- Performance optimizations
- Additional validation metrics
- Enhanced visualizations

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 👨‍💻 Author

**Dip Kumar Ghosh**  
Computational Chemist & Drug Discovery Scientist

## 🙏 Acknowledgments

- Built for the computational drug discovery community
- Inspired by established ADMET prediction tools and literature
- Thanks to contributors and testers

---

**Note**: This tool is designed for research purposes. Always validate computational predictions with experimental data in drug discovery workflows.
