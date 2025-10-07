from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import os

# Change working directory to the script's directory
os.chdir(os.path.dirname(__file__))

# Parse FASTA and calculate GC content
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


# Create DataFrame
df = pd.DataFrame({"Bacteria": bact_names, "GC content(%)": gc_content})

# Plotting
fig, (ax_table, ax_bar) = plt.subplots(
    1, 2, figsize=(12, 6))
ax_table.axis("off")

table = ax_table.table(
    cellText=df.values,
    colLabels=df.columns,
    colColours=["#60D1C9"] * len(df.columns),
    loc="center",
    cellLoc="center"
)
table.scale(1, 2.5)
table.auto_set_font_size(False)
table.set_fontsize(11)
for (row, col), cell in table.get_celld().items():
    if col == 0 and row != 0:
        cell.get_text().set_fontstyle("italic")

ax_bar.bar(
    df["Bacteria"],
    df["GC content(%)"],
    color="#60D1C9",
    edgecolor="black",
    linewidth=1.0,
    width=0.6,
    alpha=0.9)
ax_bar.set_xticks(range(len(df)))
ax_bar.set_xticklabels(df["Bacteria"], rotation=90,
                       fontsize=12, fontstyle='italic')
ax_bar.set_ylabel("GC content (%)", fontsize=11)
for i, value in enumerate(df["GC content(%)"]):
    ax_bar.text(
        i,  # x-coordinate (index of the bar)
        value + 1,  # y-coordinate (slightly above the bar)
        f"{round(value, 2)}%",
        ha="center",
        va="bottom",
        fontsize=10
    )
ax_bar.set_ylim(0, max(df["GC content(%)"]) + 10)
ax_bar.set_title("GC content of 16S rRNA genes", fontsize=13)
ax_bar.tick_params(axis="y", labelsize=10)

plt.tight_layout()
plt.show()
