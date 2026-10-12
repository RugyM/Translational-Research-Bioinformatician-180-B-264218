
import pandas as pd
import re
from collections import Counter
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(
    description="Summarise Sambamba exon coverage by gene."
)

parser.add_argument("-i", "--input",     required=True, help="Input Sambamba file")
parser.add_argument("-g", "--genelist", required=True, help="Panel gene list (one gene per line)")
parser.add_argument("-t", "--threshold", required=False, type=float, default=100.0, help="Minimum percentage of bases covered at 30x (default: 100)")
parser.add_argument("-o", "--output",    required=False, default="sample.output.txt",        help="Output report filename")
parser.add_argument("-m", "--missing",   required=False, default="sample_missing_genes.txt", help="Missing genes output filename (default: sample_missing_genes.txt)")
args = parser.parse_args()


input_file = args.input

header = ["chromosome", "StartPosition", "EndPosition", 
          "FullPosition", "NotUsed", "NotUsed2", 
          "GeneSymbol;Accession", "Size", "readCount", 
          "meanCoverage", "percentage30", "sampleName"]
df = pd.read_csv(input_file,sep="\t",skiprows=1,header=None,names=header)


for line, position in zip(df["GeneSymbol;Accession"], df["FullPosition"]):
    if "GENE_UNKNOWN" in line:
        print(f"Warning: GENE_UNKNOWN at {position} - Skipping this entry.")
        continue
    gene_symbol = line.split(";")[0].strip()
    if gene_symbol in ["."]:
        print(f"Warning: Gene symbol missing at {position} - Skipping this entry.")
        continue
    df.loc[df["GeneSymbol;Accession"] == line, "GeneSymbol"] = gene_symbol
df = df.dropna(subset=["GeneSymbol"])



# Convert coverage to numeric
df["percentage30"] = pd.to_numeric(df["percentage30"],errors="coerce")

invalid_rows = df[df["percentage30"].isna()]
for _, row in invalid_rows.iterrows():
    print(f"Warning: Invalid coverage at {row['FullPosition']} - Skipping entry.")

df = df.dropna(subset=["percentage30"])

# Identify failing exons
MIN_PERCENTAGE30 = args.threshold
df["Failing"] = df["percentage30"] < MIN_PERCENTAGE30


'''
From README:
Note that some genes are represented by more than one transcript 
accession in the panel — these should be treated as a single gene in your report.
'''
df = df.drop_duplicates(subset=["sampleName", "GeneSymbol", "chromosome","StartPosition", "EndPosition"])


gene_df = (
    df.groupby(["sampleName", "GeneSymbol"], as_index=False)
      .agg(
          TotalExons=("StartPosition", "count"),
          FailingExons=("Failing", "sum"),
          MinCoverage=("percentage30", "min")
      )
)

df["ExonRegion"] = (df["chromosome"].astype(str) + ":"+df["StartPosition"].astype(str) + "-"+df["EndPosition"].astype(str))

worst_indices = (df.groupby(["sampleName", "GeneSymbol"])["percentage30"].idxmin())
worst_exons = df.loc[worst_indices,["sampleName", "GeneSymbol", "ExonRegion"]]

gene_df = gene_df.merge(
    worst_exons,
    on=["sampleName", "GeneSymbol"],
    how="left"
)

gene_df["Status"] = "PASS"
gene_df.loc[gene_df["FailingExons"] > 0,"Status"] = "FAIL"


failing_genes = gene_df[gene_df["Status"] == "FAIL"].copy()


failing_genes = failing_genes.rename(columns={
    "MinCoverage": "MinCoverage%",
    "ExonRegion": "WorstExon(chr:start-end)"
})

report_columns = [
    "GeneSymbol",
    "TotalExons",
    "FailingExons",
    "MinCoverage%",
    "WorstExon(chr:start-end)",
    "Status"
]

with open(args.genelist) as f:
    panel_genes = {
        line.strip()
        for line in f
        if line.strip()
    }

for sample, sample_genes in gene_df.groupby("sampleName"):

    # make output directory for sample - avoid overwriting if multiple samples
    sample_dir = Path(f"sample{sample}")
    sample_dir.mkdir(parents=True, exist_ok=True)

    failing_exons = df[(df["sampleName"] == sample) & (df["Failing"])]

    # write failing exon positions
    with open(sample_dir / "failing_positions.txt", "w") as f:
        for gene, exons in failing_exons.groupby("GeneSymbol"):
            f.write(f"{gene}\n")
            
            for region in exons["ExonRegion"]:
                f.write(f"{region}\n")
            f.write("\n")


    sample_failures = failing_genes[failing_genes["sampleName"] == sample]

    # write gene coverage report
    sample_failures[report_columns].to_csv(
        sample_dir / Path(args.output).name,
        sep="\t",
        index=False
    )

    observed_genes = set(sample_genes["GeneSymbol"].dropna())
    missing_genes = sorted(panel_genes - observed_genes)

    # write missing genes
    with open(sample_dir / Path(args.missing).name, "w") as f:
        for gene in missing_genes:
            f.write(gene + "\n")

    # summary
    print(
        f"Sample {sample}: "
        f"{len(sample_failures)}/{len(panel_genes)} "
        "genes below 30x threshold"
    )
    
"""
Part B: Q9
I use minimum aggregation here - any exon below the coverage threshold results 
in the gene failing the coverage requirement. Weighted mean can allow small exons
with < threshold coverage to be masked by larger exons with good coverage, resulting 
in a gene passing the coverage requirement.

Diagnostically, clinically significant variants may be in small exons with coverage below
threshold, and these would be missed if a gene is reported as passing coverage requirements. 
"""

