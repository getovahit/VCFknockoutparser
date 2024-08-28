# VCFknockoutparser
Script to parse a single sample VCF file and determine gene knockouts 
# VCF Gene Knockout Parser

## Table of Contents
1. [Introduction](#introduction)
2. [Features](#features)
3. [Prerequisites](#prerequisites)
4. [Installation](#installation)
5. [Usage](#usage)
6. [How It Works](#how-it-works)
7. [Output](#output)
8. [Running Tests](#running-tests)
9. [Troubleshooting](#troubleshooting)
10. [Contributing](#contributing)
11. [License](#license)

## Introduction

The VCF Gene Knockout Parser is a Python script designed to analyze Variant Call Format (VCF) files from whole genome sequencing data. It identifies potential gene knockouts by examining genetic variants and their predicted effects on gene function. This tool is particularly useful for researchers and clinicians working in genetics and genomics who need to quickly identify potentially significant mutations in large genomic datasets.

## Features

- Processes VCF files to identify potential gene knockouts
- Utilizes the Ensembl Variant Effect Predictor (VEP) for comprehensive variant annotation
- Handles multiple alternative alleles
- Distinguishes between heterozygous and homozygous knockouts
- Processes VCF files in batches by chromosome to optimize performance and memory usage
- Includes basic unit tests for core functions

## Prerequisites

Before you begin, ensure you have the following installed on your system:

- Python 3.6 or higher
- pip (Python package installer)
- git (for cloning the repository)
- Perl (required for VEP installation)

## Installation

1. Clone the repository:
   ```
   git clone https://github.com/getovahit/VCFknockoutparser.git
   cd vcf-gene-knockout-parser
   ```

2. Create a virtual environment (optional but recommended):
   ```
   python -m venv venv
   source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
   ```

3. Install required Python packages:
   ```
   pip install pyvcf requests
   ```

4. Install Ensembl VEP:
   - Follow the official VEP installation guide: [VEP Installation](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html)
   - Basic steps (may vary based on your system):
     ```
     git clone https://github.com/Ensembl/ensembl-vep.git
     cd ensembl-vep
     perl INSTALL.pl
     ```
   - During the VEP installation, you'll be prompted to install cache files. Install the cache for the human genome (and any other species you're interested in).
   - Ensure that the `vep` command is accessible from your system's PATH.

5. Download necessary VEP cache files:
   - The cache files contain pre-computed variant effect predictions and are essential for VEP to run efficiently.
   - Download the cache for the human genome (or other relevant species) as guided by the VEP installation process.

## Usage

To run the VCF Gene Knockout Parser:

```
python vcf_gene_knockout_parser.py path/to/your/input.vcf
```

Replace `path/to/your/input.vcf` with the path to your VCF file.

## How It Works

1. The script first splits the input VCF file into separate files by chromosome.
2. For each chromosome file:
   a. It runs VEP to annotate the variants.
   b. It processes the VEP output to identify potential knockouts based on variant consequences.
   c. It determines zygosity (heterozygous vs. homozygous) for each potential knockout.
3. Results from all chromosomes are combined and summarized.
4. The script outputs a list of all potential gene knockouts and a separate list of homozygous knockouts.

## Output

The script provides two main outputs:

1. A list of all potential gene knockouts, including both heterozygous and homozygous variants.
2. A list of homozygous gene knockouts, which are more likely to result in complete loss of gene function.

For each gene, the output includes the gene symbol and the associated variant(s) that may cause a knockout.

## Running Tests

To run the unit tests:

```
python vcf_gene_knockout_parser.py test
```

This will execute the basic unit tests for the core functions of the script.

## Troubleshooting

- If you encounter a "command not found" error for VEP, ensure that it's correctly installed and added to your system's PATH.
- For memory issues with large VCF files, try processing chromosomes individually by modifying the script to accept a chromosome number as an input parameter.
- If you're having issues with VEP annotations, check that you've downloaded the appropriate cache files for your genome build.

## Contributing

Contributions to improve the VCF Gene Knockout Parser are welcome. Please feel free to submit pull requests or open issues to discuss potential improvements.

## License

[Specify your chosen license here, e.g., MIT, GPL, etc.]
