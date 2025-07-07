# Canonical Exon Annotation

This script extracts exon numbers from canonical transcripts using MANE-Select data and UCSC hg38 annotation to provide high-confidence exon coordinates.

## Overview

The script combines two authoritative genomic data sources:
- **MANE-Select**: NCBI's Matched Annotation from NCBI and EMBL-EBI for canonical transcript identification
- **UCSC hg38**: Human genome reference annotation for exon coordinate mapping

## Workflow

1. **Download MANE-Select data**: Retrieves the official canonical transcript set from NCBI
2. **Process transcript identifiers**: Removes version numbers from RefSeq IDs and standardizes chromosome prefixes
3. **Create genomic ranges**: Converts data to GRanges objects for efficient genomic operations
4. **Extract exon information**: Uses UCSC hg38 annotation to obtain precise exon coordinates
5. **Filter to canonical transcripts**: Retains only high-confidence MANE-Select transcripts
6. **Export results**: Saves processed data as CSV and TXT files

## Output

The script generates output files in the working directory:
- `exons_mane.csv`: Exon coordinates for canonical transcripts in CSV format
- `exons_mane.txt`: Same data in TXT format for broader compatibility

## Requirements

- R environment with genomic annotation packages
- Internet connection for downloading MANE-Select data
- UCSC hg38 genome annotation access

## Usage

Run the script to automatically download, process, and export canonical exon annotations to the working directory.

