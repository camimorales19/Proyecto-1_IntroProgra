# app.py
import io
from pathlib import Path
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from Bio import SeqIO

# AI imports
from google import genai
from dotenv import load_dotenv
import os

# Load API key from .env
load_dotenv()
# client connects to google AI servers and accesses API key
client = genai.Client(api_key=os.getenv("GEMINI_API_KEY"))

# Streamlit config
st.set_page_config(page_title="GC 16S rRNA", layout="centered")
st.title("GC content – 16S rRNA")

st.markdown(
    "Sube uno o varios archivos **FASTA** con secuencias de **16S rRNA**. "
    "Calcularemos el **%GC** por especie y mostraremos una **tabla**, un **gráfico**, "
    "y un **análisis generado por IA** para identificar la bacteria con mayor y menor %GC."
)

# -----------------------------------------------------------
# FUNCTIONS


def gc_percent(seq: str) -> float:  # function returns a float value
    s = seq.upper()  # avoids case sensitivity
    g = s.count("G")
    c = s.count("C")
    # total count of A, T, G, C (ignores N and other letters)
    atgc = sum(s.count(x) for x in "ATGC")
    # returns 0.0 if atgc is 0 (to avoid division by zero)
    return round(100 * (g + c) / atgc, 2) if atgc else 0.0

# TODO: manage sequences with ambiguous letters (should they be ignored?)


def name_from_record(r) -> str:
    parts = r.description.split()
    return f"{parts[1]} {parts[2]}" if len(parts) >= 3 else r.id
# captures first two words of description (bacteria name),else, uses record ID


def ask_gemini(df):
    # sends dataframe as text to Gemini AI for analysis
    table_as_text = df.to_string(index=False)

    prompt = f"""
    Tengo esta tabla con porcentajes de GC de secuencias 16S rRNA:

    {table_as_text}

    Necesito que me digas:
    1. Cuál bacteria tiene el mayor %GC.
    2. Cuál bacteria tiene el menor %GC.
    3. Devuélvelo en máximo 4 líneas, muy claro y conciso.
    """

    response = client.models.generate_content(
        model="gemini-2.5-flash",
        contents=prompt
    )
    # generates content based on the prompt using Gemini model
    return response.text.strip()


# -----------------------------------------------------------
# SIDEBAR SECTION

with st.sidebar:
    st.header("Parámetros")
    uploaded = st.file_uploader(
        "Sube archivos FASTA (.fa, .fasta, .fna)",
        type=["fa", "fasta", "fna"],
        accept_multiple_files=True
    )
    sort_order = st.radio("Ordenar por %GC", [
                          "Descendente", "Ascendente"], index=0)
    show_labels = st.checkbox("Mostrar etiquetas en las barras", value=True)

# -----------------------------------------------------------
# MAIN LOGIC

if uploaded:
    rows = []
    for f in uploaded:
        try:
            data = f.read().decode("utf-8", errors="ignore")
            # iterate over sequences in the fasta file
            for r in SeqIO.parse(io.StringIO(data), "fasta"):
                rows.append((name_from_record(r), gc_percent(str(r.seq))))
        except Exception as e:
            st.error(f"No se pudo leer {f.name}: {e}")

    if not rows:
        st.warning(
            "No se encontraron secuencias válidas en los archivos subidos.")
        st.stop()

    # Convert to DataFrame
    df = (
        pd.DataFrame(data=rows, columns=["Bacteria", "GC content (%)"])
        .groupby("Bacteria", as_index=False)
        .mean()
        .round(2)
    )

    df = df.sort_values(
        "GC content (%)",
        # If sort_order is "Ascendente", ascending is True, else False ("Descendente")
        ascending=(sort_order == "Ascendente"),
        ignore_index=True
    )

    # -----------------------------------------------------------
    #  TABLE SECTION
    st.subheader("Tabla comparativa")
    st.dataframe(df, use_container_width=True)

    # -----------------------------------------------------------
    # PLOT SECTION

    st.subheader("Gráfico comparativo")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(df["Bacteria"], df["GC content (%)"])
    ax.set_ylabel("GC content (%)")
    ax.set_title("GC content of 16S rRNA genes")
    ax.tick_params(axis="x", rotation=75)

    if show_labels:
        for i, v in enumerate(df["GC content (%)"]):
            ax.text(i, v + 0.5, f"{v:.2f}%",
                    ha="center", va="bottom", fontsize=8)

    st.pyplot(fig, clear_figure=True)

    # -----------------------------------------------------------
    # AI ANALYSIS SECTION

    st.subheader("Análisis por IA (Gemini)")

    try:
        ai_result = ask_gemini(df)
        st.info(ai_result)
    except Exception as e:
        st.error(f"Error al consultar a Gemini: {e}")

    # -----------------------------------------------------------
    #  DOWNLOADS SECTION

    # CSV Download
    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        "Descargar tabla (CSV)",
        data=csv,
        file_name="gc_table.csv",
        mime="text/csv",
    )

    # PNG Download
    buf = io.BytesIO()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(df["Bacteria"], df["GC content (%)"])
    ax.set_ylabel("GC content (%)")
    ax.set_title("GC content of 16S rRNA genes")
    ax.tick_params(axis="x", rotation=75)
    if show_labels:
        for i, v in enumerate(df["GC content (%)"]):
            ax.text(i, v + 0.5, f"{v:.2f}%",
                    ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(buf, format="png", dpi=300)

    st.download_button(
        "Descargar gráfico (PNG)",
        data=buf.getvalue(),
        file_name="gc_plot.png",
        mime="image/png",
    )

else:
    st.info("Sube los archivos FASTA en la barra lateral para calcular el %GC.")
