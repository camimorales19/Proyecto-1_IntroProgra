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
    names_list = (bact_names.append(joined_name))
    gc_list = gc_content.append(gc_perc)

df = pd.DataFrame({"Bacteria": bact_names, "GC_percent": gc_content})
fig, ax = plt.subplots()
ax.axis("off")

tabla = ax.table(
    cellText=df.values,
    colLabels=df.columns,
    loc="center"
)

plt.show()
