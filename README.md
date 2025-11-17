# Phasing Pipeline


<img src="imgs/logo.png" alt="Phasing Pipeline" width="50%"/>


A Nextflow pipeline for phasing unphased genotype data using Beagle with reference panels from the 1000 Genomes Project.

## Overview

This pipeline performs haplotype phasing of unphased VCF files using Beagle, a powerful tool for phasing and imputation. The pipeline includes preprocessing steps to prepare both the unphased data and reference panels, performs chromosome-by-chromosome phasing, and postprocesses the results.

## Features

- **Preprocessing**: Quality control and preparation of unphased VCF files
- **Reference Panel Processing**: Indexing and preparation of 1000 Genomes reference panels
- **Phasing**: Chromosome-by-chromosome phasing using Beagle
- **Postprocessing**: Concatenation and final processing of phased results
- **Containerized**: Uses Singularity containers for reproducibility
- **Scalable**: Supports SLURM cluster execution

## Requirements

- Nextflow (>= 22.04.0)
- Singularity (for containerized execution)
- SLURM (for cluster execution, optional)

## Quick Start

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd Phasing
   ```

2. **Configure parameters**:
   Edit `params/params_beagle.yml` with your input files:
   ```yaml
   vcf_unphased: ./data/your_unphased.vcf.gz
   refcsv: ./params/files_beagle.csv
   outdir: ./results/
   output_prefix: phased_output
   ```

3. **Run the pipeline**:
   ```bash
   # Local execution
   nextflow run main.nf -params-file params/params_beagle.yml -profile local

   # SLURM cluster execution
   nextflow run main.nf -params-file params/params_beagle.yml -profile kutral
   ```

## Input Files

### Required

- **Unphased VCF**: A VCF/BCF file containing unphased genotypes (must be indexed)
- **Reference CSV**: A CSV file with columns:
  - `chr`: Chromosome identifier
  - `ref_vcf`: Path to reference VCF file for that chromosome
  - `ref_vcf_index`: Path to reference VCF index file
  - `gmap`: Path to genetic map file

### Example Reference CSV (`files_beagle.csv`):
```csv
chr,ref_vcf,ref_vcf_index,gmap
1,/path/to/chr1_ref.vcf.gz,/path/to/chr1_ref.vcf.gz.csi,/path/to/chr1.gmap
2,/path/to/chr2_ref.vcf.gz,/path/to/chr2_ref.vcf.gz.csi,/path/to/chr2.gmap
...
```

## Pipeline Workflow

1. **Preprocessing**:
   - Index VCF files
   - Fill AC (allele count) annotations
   - Remove duplicate variants
   - Remove missing genotypes
   - Prepare reference panels

2. **Phasing**:
   - Extract chromosome-specific regions
   - Run Beagle phasing with reference panels
   - Index phased output

3. **Postprocessing**:
   - Concatenate chromosome-specific results
   - Generate final phased VCF

## Output

The pipeline generates:
- Phased VCF files per chromosome: `phased_<chr>.vcf.gz`
- Concatenated phased VCF: `<output_prefix>.vcf.gz`
- Pipeline execution reports in `pipeline_info/`

## Configuration

### Profiles

- **local**: For local execution with Singularity
- **kutral**: For SLURM cluster execution on the `ngen-ko` queue

### Resource Requirements

Default resource allocation:
- Memory: 120GB per process
- CPUs: 16 per process

Adjust in `nextflow.config` if needed.

## Parameters

Key parameters (set in `params/params_beagle.yml`):

| Parameter | Description | Default |
|----------|-------------|---------|
| `vcf_unphased` | Path to unphased VCF file | - |
| `refcsv` | Path to reference CSV file | - |
| `outdir` | Output directory | `./results/` |
| `output_prefix` | Prefix for output files | - |

## Tools Used

- **Beagle**: Haplotype phasing and imputation
- **bcftools**: VCF/BCF manipulation and indexing
- **Nextflow**: Workflow orchestration

## Citation

If you use this pipeline, please cite:
- Beagle: Browning, B. L., & Browning, S. R. (2016). Genotype imputation with millions of reference samples. *The American Journal of Human Genetics*, 98(1), 116-126.

## License

See [LICENSE](LICENSE) file for details.

## Author

Gabriel Cabas

## Support

For issues and questions, please open an issue on the repository.
