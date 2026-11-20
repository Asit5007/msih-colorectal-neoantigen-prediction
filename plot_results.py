import os
import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    csv_path = "results/top_neoantigen_candidates.csv"
    out_path = "results/neoantigen_summary.png"

    df = pd.read_csv(csv_path)
    if df.empty:
        print("No rows found in results/top_neoantigen_candidates.csv")
        return

    df_sorted = df.sort_values(by="mhcflurry_affinity", ascending=True)
    top15 = df_sorted.head(15).copy()
    top15["label"] = (
        top15["gene"].astype(str)
        + " ("
        + top15["peptide_9mer"].astype(str)
        + " | "
        + top15["allele"].astype(str)
        + ")"
    )

    gene_counts = df["gene"].value_counts().sort_values(ascending=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Left: Top 15 strongest binders (lowest IC50 nM)
    ax1.barh(top15["label"][::-1], top15["mhcflurry_affinity"][::-1], color="#2b6cb0")
    ax1.set_xlabel("Predicted Affinity IC50 (nM) - Lower = Stronger")
    ax1.set_title("Top 15 Strongest Candidate 9-mers")
    ax1.invert_yaxis()
    ax1.grid(axis="x", linestyle="--", alpha=0.4)

    # Right: Binder counts per gene
    ax2.barh(gene_counts.index, gene_counts.values, color="#2f855a")
    ax2.set_xlabel("Number of Binders (<= 500 nM)")
    ax2.set_title("Candidate Binders per Gene")
    ax2.grid(axis="x", linestyle="--", alpha=0.4)

    plt.tight_layout()
    os.makedirs("results", exist_ok=True)
    plt.savefig(out_path, dpi=300)
    print(f"Saved chart to {out_path}")


if __name__ == "__main__":
    main()
