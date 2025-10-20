# app.py
import io
from pathlib import Path
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from Bio import SeqIO

st.set_page_config(page_title="GC 16S rRNA", layout="centered")
st.title("GC content – 16S rRNA")

st.markdown(
    "Sube uno o varios archivos **FASTA** con secuencias de **16S rRNA**. "
    "Calcularemos el **%GC** por especie y mostraremos una **tabla** y un **gráfico**."
)

def gc_percent(seq: str) -> float:
    s = seq.upper()
    g = s.count("G")
    c = s.count("C")
    atgc = sum(s.count(x) for x in "ATGC")
    return round(100 * (g + c) / atgc, 2) if atgc else 0.0

def name_from_record(r) -> str:
    parts = r.description.split()
    return f"{parts[0]} {parts[1]}" if len(parts) >= 2 else r.id

with st.sidebar:
    st.header("Parámetros")
    uploaded = st.file_uploader(
        "Sube archivos FASTA (.fa, .fasta, .fna)",
        type=["fa", "fasta", "fna"],
        accept_multiple_files=True
    )
    sort_order = st.radio("Ordenar por %GC", ["Descendente", "Ascendente"], index=0)
    show_labels = st.checkbox("Mostrar etiquetas en las barras", value=True)

if uploaded:
    rows = []
    for f in uploaded:
        try:
            data = f.read().decode("utf-8", errors="ignore")
            for r in SeqIO.parse(io.StringIO(data), "fasta"):
                rows.append((name_from_record(r), gc_percent(str(r.seq))))
        except Exception as e:
            st.error(f"No se pudo leer {f.name}: {e}")

    if not rows:
        st.warning("No se encontraron secuencias válidas en los archivos subidos.")
        st.stop()

    df = (
        pd.DataFrame(rows, columns=["Bacteria", "GC content (%)"])
        .groupby("Bacteria", as_index=False)
        .mean()
        .round(2)
    )
    df = df.sort_values(
        "GC content (%)",
        ascending=(sort_order == "Ascendente"),
        ignore_index=True
    )

    st.subheader("Tabla comparativa")
    st.dataframe(df, use_container_width=True)

    st.subheader("Gráfico comparativo")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(df["Bacteria"], df["GC content (%)"])
    ax.set_ylabel("GC content (%)")
    ax.set_title("GC content of 16S rRNA genes")
    ax.tick_params(axis="x", rotation=75)

    if show_labels:
        for i, v in enumerate(df["GC content (%)"]):
            ax.text(i, v + 0.5, f"{v:.2f}%", ha="center", va="bottom", fontsize=8)

    st.pyplot(fig, clear_figure=True)

    # Descargas
    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        "⬇️ Descargar tabla (CSV)",
        data=csv,
        file_name="gc_table.csv",
        mime="text/csv",
    )

    buf = io.BytesIO()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(df["Bacteria"], df["GC content (%)"])
    ax.set_ylabel("GC content (%)")
    ax.set_title("GC content of 16S rRNA genes")
    ax.tick_params(axis="x", rotation=75)
    if show_labels:
        for i, v in enumerate(df["GC content (%)"]):
            ax.text(i, v + 0.5, f"{v:.2f}%", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(buf, format="png", dpi=300)
    st.download_button(
        "⬇️ Descargar gráfico (PNG)",
        data=buf.getvalue(),
        file_name="gc_plot.png",
        mime="image/png",
    )

else:
    st.info("Sube los archivos FASTA en la barra lateral para calcular el %GC.")