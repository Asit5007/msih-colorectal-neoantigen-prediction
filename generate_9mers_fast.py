import json
import os
import time
import urllib.request
from Bio.Seq import Seq
import pandas as pd

os.makedirs("data", exist_ok=True)
os.makedirs("results", exist_ok=True)

recurrent = pd.read_csv("results/recurrent_fs_indels.csv").head(15)
msi_samples = pd.read_csv("data/msi_h_samples.tsv", sep="\t")["sampleId"].unique()
mutations = pd.read_csv("data/mutations.tsv", sep="\t", comment="#", low_memory=False)
msi_fs = mutations[
    mutations["Tumor_Sample_Barcode"].isin(msi_samples)
    & mutations["Variant_Classification"].isin(["Frame_Shift_Del", "Frame_Shift_Ins"])
]

cds_cache = {}


def fetch_canonical_cds(gene_symbol: str):
    if gene_symbol in cds_cache:
        return cds_cache[gene_symbol]

    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=10) as resp:
        data = json.loads(resp.read().decode())

    for tx in data.get("Transcript", []):
        if tx.get("is_canonical") == 1:
            seq_url = f"https://rest.ensembl.org/sequence/id/{tx['id']}?type=cds"
            seq_req = urllib.request.Request(seq_url, headers={"Content-Type": "text/plain"})
            with urllib.request.urlopen(seq_req, timeout=10) as sresp:
                result = (tx["id"], sresp.read().decode().strip())
                cds_cache[gene_symbol] = result
                time.sleep(0.2)
                return result

    cds_cache[gene_symbol] = (None, None)
    return None, None


cds_records = []
nine_mers = []

for idx_row, row in recurrent.iterrows():
    gene = row["Hugo_Symbol"]
    hgvsp = row["HGVSp_Short"]
    print(f"[{idx_row + 1}/{len(recurrent)}] {gene} ({hgvsp})...", flush=True)
    sample_mut = msi_fs[(msi_fs["Hugo_Symbol"] == gene) & (msi_fs["HGVSp_Short"] == hgvsp)].iloc[0]
    hgvsc = str(sample_mut.get("HGVSc", ""))

    try:
        tx_id, cds_str = fetch_canonical_cds(gene)
        if not cds_str:
            print(f"  Skipping {gene} ({hgvsp}): no canonical CDS", flush=True)
            continue
        cds_records.append(f">{gene}|{tx_id}\n{cds_str}\n")

        pos_digits = "".join(ch for ch in hgvsc.split(".")[1] if ch.isdigit()) if "." in hgvsc else ""
        if not pos_digits:
            continue
        c_pos = int(pos_digits) - 1

        if "del" in hgvsc:
            del_bases = hgvsc.split("del")[1]
            del_len = len(del_bases) if del_bases.isalpha() else 1
            mut_dna = cds_str[:c_pos] + cds_str[c_pos + del_len :]
        elif "dup" in hgvsc:
            dup_bases = hgvsc.split("dup")[1]
            ins_seq = dup_bases if dup_bases.isalpha() else cds_str[c_pos]
            mut_dna = cds_str[: c_pos + 1] + ins_seq + cds_str[c_pos + 1 :]
        elif "ins" in hgvsc:
            ins_seq = hgvsc.split("ins")[1]
            mut_dna = cds_str[: c_pos + 1] + ins_seq + cds_str[c_pos + 1 :]
        else:
            continue

        wt_prot = str(Seq(cds_str[: len(cds_str) - len(cds_str) % 3]).translate(to_stop=True))
        mut_prot = str(Seq(mut_dna[: len(mut_dna) - len(mut_dna) % 3]).translate(to_stop=True))

        diff_idx = next(
            (i for i, (a, b) in enumerate(zip(wt_prot, mut_prot)) if a != b),
            min(len(wt_prot), len(mut_prot)),
        )
        start_idx = max(0, diff_idx - 8)
        novel_region = mut_prot[start_idx:]

        for i in range(len(novel_region) - 8):
            pep = novel_region[i : i + 9]
            if pep not in wt_prot:
                nine_mers.append((gene, hgvsp, int(row["Patient_Count"]), pep))
    except Exception as e:
        print(f"  Skipping {gene} ({hgvsp}): {e}", flush=True)

with open("data/wildtype_cds.fasta", "w") as f:
    f.writelines(cds_records)

with open("results/candidate_9mers.fasta", "w") as f:
    for idx, (gene, hgvsp, count, pep) in enumerate(nine_mers, 1):
        f.write(f">{gene}_{hgvsp}_pts{count}_pep{idx}\n{pep}\n")

print(f"Saved {len(cds_records)} wild-type CDS sequences to data/wildtype_cds.fasta")
print(f"Generated {len(nine_mers)} novel 9-mer peptides in results/candidate_9mers.fasta")
