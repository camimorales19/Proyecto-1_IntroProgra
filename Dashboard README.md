# Dashboard README

## Integrantes
- Camila Morales
- Belén Angulo
- Luis Eduardo Muñoz

## Descripción
Este README documenta **el dashboard en Streamlit** para calcular el **%GC** de **genes 16S rRNA** a partir de archivos **FASTA**.  
Permite:
- Subir uno o varios archivos `.fasta` / `.fa` / `.fna`.
- Calcular el %GC por especie (promedia si hay múltiples registros).
- Visualizar **tabla** y **gráfico de barras**.
- Descargar **CSV** (tabla) y **PNG** (gráfico).

## Requisitos
- Python 3.9+ (recomendado)
- Dependencias (instalar con `requirements.txt`):  
  `streamlit`, `biopython`, `pandas`, `matplotlib`

## Instalación rápida (Mac/Linux/Windows PowerShell)
```bash
python3 -m venv .venv
source .venv/bin/activate          # En Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Ejecutar el dashboard
```bash
streamlit run app.py
```
1. En la **barra lateral**, sube los FASTA.
2. Revisa la **tabla** con el %GC por especie.
3. Observa el **gráfico** (opción de etiquetas y orden).
4. Usa los botones para **descargar** CSV/PNG.

## Estructura sugerida del repositorio
```
.
├─ app.py                 # Dashboard Streamlit
├─ requirements.txt       # Dependencias
├─ Dashboard README.md    # Este archivo
├─ .gitignore             # Ignora .venv/ y archivos generados
└─ data/ (opcional)       # FASTA de prueba
```

## Notas técnicas
- El cálculo de GC cuenta únicamente caracteres **A/T/G/C** (ignora ambiguos en el denominador).
- El nombre de especie se infiere como “Género especie” desde la descripción del registro; si no, se usa `id`.
- El gráfico se genera con **matplotlib** (sin estilos adicionales, como requiere la rúbrica).

## Entrega sugerida
- Enviar `app.py` + `requirements.txt` + este `Dashboard README.md` (o link al repo).  
- Si se pide notebook, exportar una versión `.ipynb` con una celda que llame a las mismas funciones.
