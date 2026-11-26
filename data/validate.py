"""Verify the released case index, chemistry, FASTA joins and published analysis subsets."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

from glycan_graph import require, resolve_ccd_companion, validate as validate_graph

DATA = Path(__file__).resolve().parent
METHODS = ("af3", "boltz2", "chai1", "protenix", "rfaa", "sweetfold")
# Counts in order: chemical inputs, distinct protein sequences, glycosylation sites,
# monosaccharide nodes, inter-sugar bonds.
EXPECTED = {
    ("A", "all"): (348, 269, 1130, 2821, 1691),
    ("B", "all"): (59, 49, 212, 1106, 894),
    ("A", "length_eligible"): (337, 261, 1045, 2717, 1672),
    ("B", "length_eligible"): (58, 48, 203, 1088, 885),
    ("A", "common_six"): (327, 256, 884, 1831, 947),
    ("B", "common_six"): (57, 47, 171, 476, 305),
}


def read_sequences(path):
    sequences = {}
    sequence_id = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            sequence_id = line[1:]
            require(sequence_id not in sequences, f"Duplicate FASTA ID {sequence_id}")
            sequences[sequence_id] = ""
        else:
            require(sequence_id is not None, "FASTA sequence without header")
            sequences[sequence_id] += line
    for sequence_id, sequence in sequences.items():
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        require(sequence_id == "GPSEQ_" + digest[:16], f"FASTA ID mismatch at {sequence_id}")
    return sequences


def keyed(rows, key):
    index = {row[key]: row for row in rows}
    require(len(index) == len(rows), f"Repeated key {key}")
    return index


def validate_release(folder=DATA):
    model_path = folder / "glycan_connections.json"
    model = json.loads(model_path.read_text())
    with gzip.open(folder / "glycan_connections_original.json.gz", "rt") as stream:
        original = json.load(stream)
    for document in (original, model):
        validate_graph(document)
        resolve_ccd_companion(document, model_path)
    require(original["ccd_components"] == model["ccd_components"], "Chemical dictionary changed")
    current = keyed(model["cases"], "case_id")
    archived = keyed(original["cases"], "case_id")
    sequences = read_sequences(folder / "protein_sequences.fasta")
    with (folder / "cases.tsv").open(newline="") as stream:
        cases = list(csv.DictReader(stream, delimiter="\t"))
    require(keyed(cases, "case_id").keys() == current.keys(), "Model case index mismatch")
    require(keyed(cases, "source_case_id").keys() == archived.keys(), "Source case index mismatch")
    require({r["sequence_id"] for r in cases} == sequences.keys(), "FASTA coverage mismatch")

    for row in cases:
        cid = row["case_id"]
        case = current[cid]
        source = archived[row["source_case_id"]]
        for field in ("datasets", "msa_query", "glycans", "additional_covalent_bonds", "reference_quality"):
            require(case[field] == source[field], f"{cid}: changed {field}")
        seq = sequences[row["sequence_id"]]
        require(case["protein"]["sequence"] == source["protein"]["sequence"] == seq, f"{cid}: sequence differs")
        require(case["msa_query"]["sequence_id"] == row["sequence_id"], f"{cid}: wrong sequence join")
        require(case["msa_query"]["query_length"] == int(row["protein_length"]) == len(seq), f"{cid}: wrong length")
        require(case["protein"]["modifications"] == [], f"{cid}: model contains a non-glycan modification")
        mods = source["protein"]["modifications"]
        for mod in mods:
            parent = {"HYP": "P", "OCS": "C"}[mod["ccd_id"]]
            require(seq[mod["position"] - 1] == parent, f"{cid}: wrong canonical parent")
        require(int(row["n_source_protein_modifications"]) == len(mods), f"{cid}: wrong modification count")
        require(bool(mods) == (cid != source["case_id"]), f"{cid}: invalid source-model ID join")
        for dataset in ("A", "B"):
            require(int(row["in_" + dataset]) == (dataset in case["datasets"]), f"{cid}: incorrect set {dataset}")
        glycans = case["glycans"]
        require(int(row["n_sites"]) == len(glycans), f"{cid}: incorrect site count")
        require(int(row["n_sugars"]) == sum(len(g["nodes"]) for g in glycans), f"{cid}: incorrect sugar count")
        require(int(row["n_internal_glycan_bonds"]) == sum(len(g["bonds"]) for g in glycans),
                f"{cid}: incorrect internal bond count")
        eligible = len(seq) < 2000
        require(int(row["length_eligible"]) == eligible, f"{cid}: incorrect length eligibility")
        states = [row[m + "_status"] for m in METHODS]
        complete = states.count("complete")
        for method, state in zip(METHODS, states):
            allowed = {"complete", "native_token_limit" if method == "chai1" else "memory_limit"}
            require(state in (allowed if eligible else {"outside_length_scope"}), f"{cid}: invalid {method} status")
        require(complete == int(row["n_methods_complete"]), f"{cid}: incorrect method count")
        require(int(row["common_six"]) == (complete == 6), f"{cid}: incorrect common-set flag")
        require(bool(row["group_id"]) == (complete == 6), f"{cid}: common-set group mismatch")

    for dataset in ("A", "B"):
        selected = [row for row in cases if row["in_" + dataset] == "1"]
        ids = [row["case_id"] for row in selected]
        require(ids == (folder / "sets" / (dataset + ".txt")).read_text().splitlines(),
                f"Set {dataset}: ID list disagrees with cases.tsv")
        for condition in ("all", "length_eligible", "common_six"):
            subset = [r for r in selected if condition == "all" or r[condition] == "1"]
            counts = (len(subset), len({r["sequence_id"] for r in subset}),
                      *(sum(int(r[k]) for r in subset) for k in ("n_sites", "n_sugars", "n_internal_glycan_bonds")))
            require(counts == EXPECTED[dataset, condition], f"Set {dataset}, {condition}: {counts}")
            print(f"Set {dataset:1}  {condition:16} {counts[0]:3} cases   {counts[1]:3} sequences   {counts[2]:4} sites")
    print("Verified all released sequences, source/model chemical graphs, atom-level connections and cohort joins.")


if __name__ == "__main__":
    validate_release()
