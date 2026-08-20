# signalp-nextflow

Nextflow pipeline that predicts signal peptides in protein sequences using SignalP versions 4, 5, and 6, and merges the results into a single indexed GFF3 file.

## Overview

This pipeline runs [SignalP](https://services.healthtech.dtu.dk/services/SignalP-6.0/) to identify signal peptide cleavage sites in predicted proteins, as part of VEuPathDB's genome annotation workflows. Signal peptide evidence is used downstream to support subcellular localization and functional annotation of gene products.

Because SignalP 6 is significantly more compute-intensive than earlier versions, the pipeline first runs SignalP 5 over the full input set and uses its scores to filter down to a smaller candidate set of proteins before running SignalP 4 and SignalP 6. The filtering step (`filterProteinsByScore.pl`) keeps proteins scoring above a configurable SignalP score cutoff, while guaranteeing that at least a configurable minimum percentage of the input proteins are retained even if their score falls below the cutoff. SignalP 6 is additionally skipped for a given batch if the filtered protein count exceeds a configured limit, since it does not scale to very large batches. Results from all three predictors are normalized to GFF3 (`fixAndCombineGff.pl`), combined, sorted, `bgzip`-compressed, and `tabix`-indexed.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- Docker or Singularity, depending on the execution profile in the runtime configuration

## Usage

```
nextflow run VEuPathDB/signalp-nextflow \
  -r main \
  -resume \
  --inputFilePath /path/to/proteins.fa \
  --fastaSubsetSize 1000 \
  --outputFileName signalP.gff3 \
  --outputDir /path/to/output \
  -C <config>
```

The pipeline has a single, unnamed entry point (no `-entry` flag is needed).

## Key parameters

| Parameter | Description |
| --- | --- |
| `inputFilePath` | Path to the input protein FASTA file. Required. |
| `fastaSubsetSize` | Number of sequences per FASTA subset chunk processed in each parallel SignalP 5 job (default `1000`). Required. Filtered output is subsequently re-chunked into batches of 200 for SignalP 4 and 6. |
| `outputFileName` | Name of the merged, sorted GFF3 output file (default `signalP.gff3`). |
| `outputDir` | Directory the final compressed/indexed output is published to (default `output` under the launch directory). |

Set via the runtime configuration (process-level `ext` settings in `nextflow.config`):

| Setting | Description |
| --- | --- |
| `ext.org` | Organism class passed to SignalP (`euk` for eukaryotes, or `gram+`/`gram-` for bacteria). |
| `ext.filter_score_cutoff` | Minimum SignalP 5 SP score a protein must have to be kept for SignalP 4/6 processing. |
| `ext.filter_min_protein_percent_cutoff` | Minimum percentage of input proteins to retain for SignalP 4/6 regardless of score. |
| `ext.protein_count_limit` | Maximum number of filtered proteins in a batch for which SignalP 6 will still be run. |

Container images and process-level resources (executor, queue, memory) are supplied via the runtime configuration passed with `-C`, allowing the same pipeline definition to run under Docker, Singularity, or an LSF cluster profile (see `conf/docker.config`, `conf/singularity.config`, and `conf/lsf.config`).

## Output

For each run, the pipeline publishes to `outputDir`:

- `<outputFileName>.gz` — a sorted GFF3 file combining signal peptide predictions from SignalP 4, 5, and 6 for the input proteins.
- `<outputFileName>.gz.tbi` — a Tabix index for the compressed GFF3 file.
