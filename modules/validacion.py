# Módulo de validación 

import pandas as pd
import numpy as np

NO_DETECTADO = {
    "undetermined",
    "no ct",
    "noct",
    "not detected",
    "nd"
}

CATEGORIAS = (
    "id",
    "caracteristica",
    "control",
    "expresion"
)


def obtener_categoria(columna, variables):
    for categoria in CATEGORIAS:
        if columna in variables.get(categoria, []):
            return categoria
    return "sin_clasificar"


def valores_ausentes(serie):
    texto = serie.astype("string").str.strip()
    return serie.isna() | texto.eq("").fillna(False)


def validar_integridad(
    df,
    variables,
    ct_min=0,
    ct_max=40
):

    if ct_min >= ct_max:
        raise ValueError(
            "El CT mínimo debe ser menor que el máximo."
        )

    reporte = {
        "identificadores": {},
        "caracteristicas": {},
        "expresion": {},
        "valores_faltantes": {},
        "advertencias": [],
        "total_muestras": int(len(df)),
        "total_columnas": int(len(df.columns)),
        "total_faltantes": 0
    }


    for columna in variables.get("id", []):

        serie = df[columna]
        faltantes = valores_ausentes(serie)

        duplicados = int(
            serie[~faltantes]
            .astype(str)
            .str.strip()
            .duplicated()
            .sum()
        )

        reporte["identificadores"][columna] = {
            "faltantes": int(faltantes.sum()),
            "duplicados": duplicados
        }

        if faltantes.any():
            reporte["advertencias"].append(
                f"ID '{columna}': existen valores faltantes."
            )

        if duplicados > 0:
            reporte["advertencias"].append(
                f"ID '{columna}': {duplicados} duplicados."
            )


    categoricas = variables.get(
        "caracteristica_categorica", []
    )

    numericas = variables.get(
        "caracteristica_numerica", []
    )

    for columna in variables.get("caracteristica", []):

        serie = df[columna]
        faltantes = valores_ausentes(serie)

        if columna in categoricas:

            valores = serie.mask(faltantes)
            frecuencias = (
                valores.value_counts(dropna=True)
                .to_dict()
            )

            reporte["caracteristicas"][columna] = {
                "tipo": "categorica",
                "faltantes": int(faltantes.sum()),
                "numero_categorias": len(frecuencias),
                "frecuencias": frecuencias
            }

        elif columna in numericas:

            numeros = pd.to_numeric(
                serie.mask(faltantes),
                errors="coerce"
            )

            invalidos = int(
                ((~faltantes) & numeros.isna()).sum()
            )

            validos = numeros.dropna()
            validos = validos[
                np.isfinite(validos)
            ]

            reporte["caracteristicas"][columna] = {
                "Tipo": "numerica",
                "Faltantes": int(faltantes.sum()),
                "Invalidos": invalidos,
                "N": int(len(validos)),
                "Promedio": validos.mean() if len(validos) else np.nan,
                "Mediana": validos.median() if len(validos) else np.nan,
                "Desviación Estándar": validos.std() if len(validos) > 1 else np.nan,
                "Valor mínimo": validos.min() if len(validos) else np.nan,
                "Valor máximo": validos.max() if len(validos) else np.nan
            }

            if invalidos:
                reporte["advertencias"].append(
                    f"Feature '{columna}': "
                    f"{invalidos} valores no numéricos."
                )

        else:
            reporte["advertencias"].append(
                f"Feature '{columna}' sin subtipo válido."
            )

    columnas_ct = (
        variables.get("control", []) +
        variables.get("expresion", [])
    )

    for columna in columnas_ct:

        serie = df[columna]
        faltantes = valores_ausentes(serie)

        texto = (
            serie.astype("string")
            .str.strip()
            .str.lower()
        )

        no_detectados = (
            texto.isin(NO_DETECTADO).fillna(False)
        )

        numeros = pd.to_numeric(
            serie.mask(faltantes | no_detectados),
            errors="coerce"
        )

        finitos = pd.Series(
            np.isfinite(numeros),
            index=serie.index
        ).fillna(False)

        invalidos = (
            (~faltantes) &
            (~no_detectados) &
            (~finitos)
        )

        fuera_rango = (
            finitos &
            ((numeros < ct_min) | (numeros > ct_max))
        )

        reporte["expresion"][columna] = {
            "categoria": obtener_categoria(columna, variables),
            "valores_numéricos": int(finitos.sum()),
            "valores_invalidos": int(invalidos.sum()),
            "valores_faltantes": int(faltantes.sum()),
            "no_detectados": int(no_detectados.sum()),
            "fuera_rango": int(fuera_rango.sum())
        }

        if invalidos.any():
            reporte["advertencias"].append(
                f"'{columna}': valores Ct inválidos."
            )

        if fuera_rango.any():
            reporte["advertencias"].append(
                f"'{columna}': Ct fuera del rango "
                f"{ct_min}–{ct_max}."
            )


    for columna in df.columns:
        n = int(valores_ausentes(df[columna]).sum())

        if n:
            reporte["valores_faltantes"][columna] = n

    reporte["total_faltantes"] = sum(
        reporte["valores_faltantes"].values()
    )

    return reporte


def generar_reporte_calidad(df, variables):

    resultados = []
    mascara_faltantes = pd.DataFrame(index=df.index)

    for columna in df.columns:

        faltantes = valores_ausentes(df[columna])
        mascara_faltantes[columna] = faltantes

        total = len(df)
        n_faltantes = int(faltantes.sum())

        resultados.append({
            "columna": columna,
            "categoria": obtener_categoria(
                columna, variables
            ),
            "total_muestras": total,
            "valores_validos": total - n_faltantes,
            "valores_faltantes": n_faltantes,
            "porcentaje_faltantes": round(
                100 * n_faltantes / total, 2
            ) if total else 0
        })

    reporte = pd.DataFrame(resultados)

    muestras_con_faltantes = int(
        mascara_faltantes.any(axis=1).sum()
    )

    resumen = {
        "total_muestras": int(len(df)),
        "muestras_con_faltantes": muestras_con_faltantes,
        "muestras_sin_faltantes": (
            int(len(df)) - muestras_con_faltantes
        ),
        "total_faltantes": int(
            mascara_faltantes.sum().sum()
        )
    }

    return reporte, resumen
