# Genetical Mapping

Version **1.1**.

Genetical Mapping is a gene-family profiling tool for comparing samples through the sequences of a chosen family's members. It pools overlapping fragments from the supplied member sequences, counts recurring fragments, and builds sample-level profiles for hierarchical clustering and heatmap visualization.

This command-line edition reads nucleotide sequences in FASTA format and exports the profiles, cluster assignments, and fragment counts.

[中文使用说明](README_zh-CN.md)

## Features

- Analyze a supplied gene-family member set, with multiple member sequences per sample.
- Build a shared k-mer pool and count fragment occurrences across the member set.
- Summarize these counts into sample-by-position profiles or sample-level k-mer frequency profiles.
- Select fragments by their pooled occurrence counts to compare different frequency ranges.
- Apply row normalization or combine reverse-complement k-mers.
- Export cluster assignments, matrices, clustering diagnostics, and an interactive heatmap for position-based profiles.
- Record parameters, input checksum, record counts, and dependency versions for each run.

## Version 1.1

This revision provides a standalone command-line edition, explicit representation and record-selection options, an offline interactive heatmap, input and output checks, and run manifests. The release includes a synthetic example and regression tests.

Run `python main.py --version` to display the version. See [CHANGELOG.md](CHANGELOG.md) for the release summary.

## Quick start

Use Python 3.12. In Windows PowerShell, open the extracted `Genetical-Mapping` directory and run:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe main.py --fasta examples/synthetic_demo.fasta --out outputs/demo --k 20 --clusters 2
```

Open `outputs/demo/interactive_heatmap.html` in a browser to view the result. The HTML file embeds the plotting library and can be viewed offline.

The example contains eight artificial sequences in four demonstration groups. Run `python generate_demo.py` to regenerate it with the fixed seed 20261027.

## Input format

For each analysis, supply the member sequences of one previously identified gene family. A sample may contribute multiple members. The default `nonempty` selection retains every nonempty member record, including repeated sequences; the other `--variant` settings provide optional selection rules.

Supply a FASTA file containing A, C, G, and T. Letter case is normalized to uppercase. Each sequence needs a unique identifier; the part before its first underscore is the group name:

```text
>DemoA_sequence1 type:complete
ACGTTGCAACGTTGCA
>DemoB_sequence1 type:complete
ACGTCGCAACGTCGCA
```

This excerpt illustrates the header format. Use the complete bundled example for the quick-start command.

Choose a k-mer length from 1 to 64. Each group needs eligible sequence windows after filtering. Use at least three groups and choose a cluster count between 2 and the number of groups minus 1. The output directory should be new or empty.

## Options

| Option | Purpose | Default |
| --- | --- | --- |
| `--fasta` | Input FASTA file | Required |
| `--out` | Output directory | Required |
| `--k` | K-mer length, 1–64 | `20` |
| `--clusters` | Number of clusters | `2` |
| `--lower` | Minimum pooled occurrence count, inclusive | `1` |
| `--upper` | Maximum pooled occurrence count, inclusive | All counts |
| `--representation` | `position_raw`, `position_l1`, or `kmer_frequency` | `position_raw` |
| `--canonical` | Combine each k-mer with its reverse complement | Off |
| `--variant` | Sequence selection rule, described below | `nonempty` |
| `--legacy-end-drop` | Set the last available position in each group to zero in positional profiles | Off |

`position_raw` sums pooled occurrence counts for the distinct k-mers present at each start position in each group. `position_l1` divides each row by its sum. `kmer_frequency` divides each group's k-mer counts by its total valid windows before occurrence filtering; distances use square-root-transformed frequencies. Positional modes use Euclidean distances. All modes use complete-linkage clustering.

Positions refer to zero-based starts in the supplied, unaligned sequences. Position-based profiles compare these input positions directly.

Sequence selection rules:

- `nonempty`: keep all nonempty sequences.
- `exact_unique`: keep one copy of each sequence within each group.
- `complete_only`: keep sequences with `type:complete` in the FASTA header.
- `complete_unique`: combine the preceding completeness and duplicate filters.
- `gene_representative`: group identifiers after removing a terminal `_i<number>` or `_i<number>.p<number>` suffix; prefer a complete sequence, then the longest sequence, then the identifier in ascending order.

For example, to use frequency profiles and reverse-complement grouping:

```powershell
.\.venv\Scripts\python.exe main.py --fasta examples/synthetic_demo.fasta --out outputs/frequency --k 20 --clusters 2 --representation kmer_frequency --canonical
```

Run `.\.venv\Scripts\python.exe main.py --help` to view the command-line options.

## Output files

| File | Contents |
| --- | --- |
| `cluster_assignments.tsv` | Group names and cluster labels |
| `distance_matrix.tsv` | Pairwise group distances |
| `feature_matrix.tsv` | Positional feature matrix |
| `feature_matrix.npz` | Sparse frequency matrix for `kmer_frequency` |
| `linkage.npy` | Hierarchical clustering linkage matrix |
| `retained_kmers.tsv` | Retained k-mers and pooled occurrence counts |
| `position_support.tsv` | Number of sequences supporting each start position |
| `run_manifest.json` | Parameters, counts, checksum, clustering diagnostics, and environment versions |
| `interactive_heatmap.html` | Offline interactive heatmap and dendrogram for positional modes |

## Run the tests

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s . -p "test_*.py" -v
```

## License

MIT. See [LICENSE](LICENSE). Dependencies retain their respective licenses.
