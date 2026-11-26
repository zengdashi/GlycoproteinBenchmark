# Benchmark inputs and glycan connectivity

Each case is one full protein chain with the covalently attached glycan fragments observed in the experimental structure. Set A was assembled from entries first released from 1 June 2023 through 25 October 2026. Set B uses the later window beginning 1 April 2026 and ending on the same day. The two windows were filtered separately. All B cases also occur in A. A prediction for a case shared by both sets contributes to both analyses.

The sequence screen excluded a PDB entry if a glycosylated protein entity had at least 30% sequence identity and at least 80% coverage of both sequences against the corresponding pre-window PDB sequences or the RFAA candidate training sequences. Redundant entries were grouped using the same thresholds. Full protein sequence, protein modifications, glycosylation sites, monosaccharide identities and covalent linkages determine the distinct chemical inputs. A sequence can therefore occur in several cases with different glycans.

| Set | Selection | Cases | Distinct protein sequences | Glycosylation sites | Sugar residues | Sugar–sugar bonds |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| A | All curated | 348 | 269 | 1,130 | 2,821 | 1,691 |
| B | All curated | 59 | 49 | 212 | 1,106 | 894 |
| A | Protein length <2,000 residues | 337 | 261 | 1,045 | 2,717 | 1,672 |
| B | Protein length <2,000 residues | 58 | 48 | 203 | 1,088 | 885 |
| A | All six methods completed | 327 | 256 | 884 | 1,831 | 947 |
| B | All six methods completed | 57 | 47 | 171 | 476 | 305 |

## Files and identifiers

| File | Contents |
| --- | --- |
| [cases.tsv](cases.tsv) | One row per curated case. Use `case_id` to join the model inputs and per-case scores. Use `source_case_id` to find the corresponding record with the original protein modifications. |
| [sets/A.txt](sets/A.txt), [sets/B.txt](sets/B.txt) | De-identified model case IDs for each full curated set. In the index, intersect `in_A` or `in_B` with `common_six` to reproduce the paper's six-method selection. |
| [protein_sequences.fasta](protein_sequences.fasta) | The 269 full, unique protein sequences. FASTA headers contain only the `GPSEQ_` identifier used by `cases.tsv`, the chemical records and the MSA query. |
| [glycan_connections.json](glycan_connections.json) | Model-ready chemical annotations, including the protein sequence, protein–root attachment, all sugar–sugar bonds and per-sugar CCD identities. These are the standardized inputs used for prediction. |
| [glycan_connections_original.json.gz](glycan_connections_original.json.gz) | Complete de-identified chemical annotations before HYP and OCS were restored to their canonical parents. Eight cases change chemical ID. Glycan nodes, attachment and internal bonds, and protein sequence positions remain the same. Both annotation files use the same [JSON schema](glycan_connections.schema.json). |
| [chemical_components.cif](chemical_components.cif) | A 48 KiB excerpt of the [wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd) for the 22 sugar and modified-protein CCDs in this collection. Contains ideal component geometry. |
| [validate.py](validate.py), [glycan_graph.py](glycan_graph.py) | Standard-library checks for the graph, atom names, sequence and case fingerprints, FASTA joins, normalization, membership and denominators in the table above. |

The prefix `GLYCO_` identifies an entire chemical input; `GPSEQ_` identifies an exact amino-acid sequence. The identifier suffix comprises the first 16 hexadecimal characters of the corresponding SHA-256 digest. The complete hashes are in the JSON. `group_id` identifies the related-input groups used for statistical resampling on the common six-method sets.

The six `*_status` columns record per-method coverage. `complete` denotes five candidates; `native_token_limit` denotes the Chai-1 native token ceiling; `memory_limit` denotes a case omitted after the high-memory retry. `outside_length_scope` marks proteins of 2,000 or more residues excluded at the MSA stage. A case enters `common_six` when all six methods completed prediction. Inclusion is independent of the candidates' chemical or structural scores.

## Reading a glycan graph

Use the normalized JSON for model inputs and score joins. Each `case.protein.sequence` is numbered from 1 to the full query length. A chain is identified only as `P1` within its case. Each entry in `case.glycans` is one attached glycan. Its `attachment` gives the protein residue position, residue/atom identity, root-sugar node/atom and bond order. Each internal edge in `bonds` gives `node_1/atom_1`, `node_2/atom_2` and bond order. Repeated CCD types represent distinct nodes. Read the named edges to recover branches.

To reconstruct the reacted heavy-atom graph, take the within-sugar atoms and bonds from `ccd_components`, remove the leaving atoms recorded on each sugar node and join the explicitly named endpoints. The `expected_reacted_heavy_atom_names` record the complete sugar heavy-atom inventory after attachment. The `observed_atom_names` describe visibility in the selected experimental reference and serve as evaluation metadata. Input chemistry uses the full reacted inventory. Glycan identities and stereochemistry follow the CCD codes, including five-character codes. Discard dictionary hydrogens and assign hydrogens after assembly when a model requires them.

For example, `GLYCO_757018a5959c6a04` contains a branched pentasaccharide attached to Asn105 ND2. The root NAG uses C1 to attach to the protein, and the BMA node branches at O3 and O6. The following example prints all named bonds from the distributed record. Run it from the repository root.

```python
import json
from pathlib import Path

doc = json.loads(Path("data/glycan_connections.json").read_text())
case = next(c for c in doc["cases"] if c["case_id"] == "GLYCO_757018a5959c6a04")
for glycan in case["glycans"]:
    print(glycan["attachment"])
    for bond in glycan["bonds"]:
        print(bond)
```

Run the data checks with `python3 -B data/validate.py`. The annotation files are frozen snapshots taken during input preparation. Their historical `summary.model_inference_performed` flags describe the snapshot; the actual prediction coverage is recorded in `cases.tsv`. These input files contain no target coordinates, PDB accessions or author chain numbering.
