import pandas as pd
from Bio import SeqIO
from mhcflurry import Class1AffinityPredictor


def main() -> None:
    fasta_file = "results/candidate_9mers.fasta"
    alleles = ["HLA-A*02:01", "HLA-A*11:01", "HLA-A*24:02", "HLA-B*07:02"]

    records = []
    for rec in SeqIO.parse(fasta_file, "fasta"):
        pep = str(rec.seq).strip().upper()
        parts = rec.id.split("_")
        gene = parts[0] if len(parts) > 0 else rec.id
        mutation = parts[1] if len(parts) > 1 else ""
        pts = parts[2].replace("pts", "") if len(parts) > 2 else ""

        for allele in alleles:
            records.append(
                {
                    "header": rec.id,
                    "gene": gene,
                    "mutation": mutation,
                    "patient_count": pts,
                    "allele": allele,
                    "peptide_9mer": pep,
                }
            )

    if not records:
        print("No peptides found in results/candidate_9mers.fasta")
        return

    df = pd.DataFrame(records)
    predictor = Class1AffinityPredictor.load()
    preds = predictor.predict(
        peptides=df["peptide_9mer"].tolist(),
        alleles=df["allele"].tolist(),
    )

    df["mhcflurry_affinity"] = preds
    df.to_csv("results/all_mhc_predictions.csv", index=False)

    top_df = df[df["mhcflurry_affinity"] <= 500.0].sort_values(
        by="mhcflurry_affinity", ascending=True
    )
    top_df.to_csv("results/top_neoantigen_candidates.csv", index=False)

    print(f"Saved {len(df)} predictions to results/all_mhc_predictions.csv")
    print(f"Saved {len(top_df)} binders (<= 500 nM) to results/top_neoantigen_candidates.csv")


if __name__ == "__main__":
    main()
