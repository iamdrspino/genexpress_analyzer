#Módulo de normalizacion

import pandas as pd
import numpy as np


def preparar_ct(serie):
#  Convierte valores CT a números, Los valores no numéricos se convierten en NaN.

    valores = pd.to_numeric(serie, errors="coerce")
    return valores.replace([np.inf, -np.inf], np.nan)


def normalizar_qpcr(
    df,
    variables,
    columna_referencia=None,
    mirnas_seleccionados=None
):

    columnas_mirna = variables.get("expresion", [])
    columnas_control = variables.get("control", [])

    if mirnas_seleccionados is not None:
        if any(
            mirna not in columnas_mirna
            for mirna in mirnas_seleccionados
        ):
            raise ValueError(
                "Se seleccionaron columnas que no "
                "están clasificadas como Expresión."
            )
        columnas_mirna = list(mirnas_seleccionados)

    if not columnas_mirna:
        raise ValueError(
            "No hay microRNAs seleccionados."
        )

    if not columnas_control:
        raise ValueError(
            "No se identificaron columnas Control."
        )

    if columna_referencia is None:
        if len(columnas_control) != 1:
            raise ValueError(
                "Selecciona un control de referencia."
            )
        columna_referencia = columnas_control[0]

    if columna_referencia not in columnas_control:
        raise ValueError(
            "La referencia debe estar clasificada "
            "como Control."
        )

    columnas_requeridas = (
        columnas_mirna + [columna_referencia]
    )

    for columna in columnas_requeridas:
        if columna not in df.columns:
            raise ValueError(
                f"No existe la columna: {columna}"
            )


    df_normalizado = df.copy()

    datos_numericos = {
        columna: preparar_ct(df[columna])
        for columna in columnas_requeridas
    }

    ct_referencia = datos_numericos[
        columna_referencia
    ]

    if ct_referencia.notna().sum() == 0:
        raise ValueError(
            "El control seleccionado no contiene "
            "valores CT numéricos válidos."
        )


    resultados = {}

    for mirna in columnas_mirna:

        ct_mirna = datos_numericos[mirna]

        delta_ct = ct_mirna - ct_referencia

        expresion_normalizada = np.power(
            2.0, -delta_ct
        )

        nombre_dct = f"dCt_{mirna}"
        nombre_expresion = f"2^-ΔCt_{mirna}"

        # Evitar sobrescribir columnas originales
        for nombre in (nombre_dct, nombre_expresion):
            if nombre in df.columns:
                raise ValueError(
                    f"La columna '{nombre}' ya existe."
                )

        df_normalizado[nombre_dct] = delta_ct
        df_normalizado[
            nombre_expresion
        ] = expresion_normalizada

        n_validos = int(delta_ct.notna().sum())

        resultados[mirna] = {
            "columna_dct": nombre_dct,
            "columna_expresion": nombre_expresion,
            "muestras_validas": n_validos,
            "muestras_no_normalizadas": (
                len(df) - n_validos
            )
        }


    resumen = {
        "metodo": "Delta Ct",
        "referencia": columna_referencia,
        "numero_muestras": int(len(df)),
        "numero_mirnas": len(columnas_mirna),
        "resultados": resultados
    }

    return df_normalizado, resumen


def evaluar_referencia(df, columna_referencia):

    if columna_referencia not in df.columns:
        raise ValueError(
            "La columna de referencia no existe."
        )

    ct = preparar_ct(
        df[columna_referencia]
    ).dropna()

    if ct.empty:
        raise ValueError(
            "No hay valores CT válidos para evaluar."
        )

    q1 = ct.quantile(0.25)
    q3 = ct.quantile(0.75)

    resumen = {
        "referencia": columna_referencia,
        "n": int(ct.count()),
        "faltantes": int(len(df) - len(ct)),
        "media": float(ct.mean()),
        "mediana": float(ct.median()),
        "desviacion_estandar": (
            float(ct.std()) if len(ct) > 1
            else np.nan
        ),
        "minimo": float(ct.min()),
        "maximo": float(ct.max()),
        "q1": float(q1),
        "q3": float(q3),
        "iqr": float(q3 - q1)
    }

    return resumen
