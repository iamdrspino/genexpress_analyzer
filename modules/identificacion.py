# Módulo de identificación de columnas

import pandas as pd

CATEGORIAS = ("id","caracteristica","control","expresion")
CARACTERISTICA = ("categorica","numerica")
CONTROL_TIPO = ("endogeno","exogeno")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def obtener_columnas(df):

    if not isinstance(df, pd.DataFrame):
        raise TypeError("Se esperaba un DataFrame.")
      
    if df.columns.has_duplicates:
        raise ValueError(
            "Existen encabezados duplicados. "
            "Deben corregirse antes de clasificar.")
            
    return df.columns.tolist()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def validar_clasificacion(df, clasificacion):

    columnas = obtener_columnas(df)

    if set(columnas) != set(clasificacion):
        raise ValueError(
            "Todas las columnas deben clasificarse."
        )

    for columna, datos in clasificacion.items():

        categoria = datos.get("categoria")
        subtipo = datos.get("subtipo")

        if categoria not in CATEGORIAS:
            raise ValueError(
                f"Categoría inválida: {columna}"
            )

        if categoria == "caracteristica":
            if subtipo not in CARACTERISTICA:
                raise ValueError(
                    f"Selecciona el tipo de caracteristica: {columna}"
                )

        if categoria == "control":
            if subtipo not in CONTROL_TIPO:
                raise ValueError(
                    f"Selecciona el tipo de Control: {columna}"
                )

    if not any(
        datos["categoria"] == "id"
        for datos in clasificacion.values()
    ):
        raise ValueError(
            "Debe existir al menos una columna de identificación / ID."
        )

    if not any(
        datos["categoria"] == "expresion"
        for datos in clasificacion.values()
    ):
        raise ValueError(
            "Debe existir al menos una columna de expresión."
        )

    return True


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

def identificar_variables(df, clasificacion):

    validar_clasificacion(df, clasificacion)

    resultado = {
        "id": [],
        "caracteristica": [],
        "control": [],
        "expresion": [],
        "caracteristica_categorica": [],
        "caracteristica_numerica": [],
        "control_endogeno": [],
        "control_exogeno": []
    }

    for columna, datos in clasificacion.items():

        categoria = datos["categoria"]
        subtipo = datos.get("subtipo")

        resultado[categoria].append(columna)

        if categoria == "caracteristica":

            if subtipo == "categorica":
                resultado["caracteristica_categorica"].append(columna)

            elif subtipo == "numerica":
                resultado["caracteristica_numerica"].append(columna)

        elif categoria == "control":

            if subtipo == "endogeno":
                resultado["control_endogeno"].append(columna)

            elif subtipo == "exogeno":
                resultado["control_exogeno"].append(columna)

    return resultado
