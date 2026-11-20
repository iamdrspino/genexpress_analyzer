# Módulo de importación
import pandas as pd
# función para cargar los datos de exprsión desde un .csv o un .xls
def cargar_archivo(archivo):
    nombre = archivo.name.lower()
    try:
        if nombre.endswith(".csv"):
            try:
                df = pd.read_csv(
                    archivo,
                    sep=None,
                    engine="python",
                    encoding="utf-8"
                )
            except UnicodeDecodeError:
                archivo.seek(0)
                df = pd.read_csv(
                    archivo,
                    sep=None,
                    engine="python",
                    encoding="latin-1"
                )

        elif nombre.endswith((".xlsx", ".xls")):
            df = pd.read_excel(archivo)

        else:
            raise ValueError(
                "Formato no compatible."
            )

        return df

    except Exception as error:
        raise ValueError(f"No se pudo importar el archivo: {error}")



def validar_datos(df):

    valores_faltantes = int(df.isna().sum().sum())

    resultados = {
        "filas": int(df.shape[0]),
        "numero_columnas": int(df.shape[1]),
        "valores_faltantes": valores_faltantes,
        "columnas_duplicadas": int(
            df.columns.duplicated().sum()
        ),
        "filas_vacias": int(df.isna().all(axis=1).sum()),
        "columnas_vacias": int(df.isna().all(axis=0).sum()),
        "nombre_columnas": df.columns.tolist()
    }

    return resultados

