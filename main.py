from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

bact_names = []
gc_content = []
for record in SeqIO.parse("16Ssequences.fasta", "fasta"):
    name = record.description.split()[1:3]
    joined_name = " ".join(name)
    seq = record.seq
    gc_count = seq.count("G") + seq.count("C")
    gc_perc = round(gc_count / len(seq) * 100, 2)
    bact_names.append(joined_name)
    gc_content.append(gc_perc)

df = pd.DataFrame({"Bacteria": bact_names, "GC_percent": gc_content})
fig, (ax_table, ax_bar) = plt.subplots(
    1, 2, figsize=(12, 6), gridspec_kw={"width_ratios": [1, 2]})
ax_table.axis("off")
table = ax_table.table(
    cellText=df.values,
    colLabels=df.columns,
    loc="center",
    cellLoc="center",
)

ax_bar.bar(
    df["Bacteria"],
    df["GC_percent"],
    color="skyblue",
    edgecolor="black",
    linewidth=1.0,
    width=0.6,
    alpha=0.9,
)

ax_bar.set_xticks(range(len(df["Bacteria"])))
ax_bar.set_xticklabels(df["Bacteria"], rotation=90, fontsize=8)
ax_bar.set_ylabel("GC content (%)", fontsize=11)
ax_bar.set_title("GC content of 16S rRNA genes", fontsize=13)
ax_bar.tick_params(axis="y", labelsize=10)

plt.tight_layout()
plt.show()
