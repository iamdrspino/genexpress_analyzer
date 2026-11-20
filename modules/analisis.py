#Módulo de análisis

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from itertools import combinations


def _numeros(df, columna):
    if columna not in df.columns:
        raise ValueError(f"No existe la columna '{columna}'.")
    return pd.to_numeric(df[columna], errors="coerce").replace([np.inf, -np.inf], np.nan)


def _validar_grupo(df, columna_grupo):
    if columna_grupo not in df.columns:
        raise ValueError(f"No existe la Feature '{columna_grupo}'.")


def calcular_expresion_relativa(df, columna_grupo, columnas_mirna):
    """Recibe nombres de columnas ΔCt y resume 2^-ΔCt por grupo."""
    _validar_grupo(df, columna_grupo)
    datos = df.copy()
    resultados = []
    for columna in columnas_mirna:
        delta = _numeros(datos, columna)
        nombre = f"ER_{columna}"
        datos[nombre] = np.power(2.0, -delta)
        for grupo, sub in datos.groupby(columna_grupo, dropna=True, observed=True):
            v = _numeros(sub, nombre).dropna()
            resultados.append({"Grupo": grupo, "miRNA": columna.removeprefix("dCt_"),
                               "n": len(v), "Media": v.mean(), "Mediana": v.median(),
                               "DE": v.std(), "Minimo": v.min(), "Maximo": v.max()})
    return datos, pd.DataFrame(resultados)


def calcular_fold_change(df, columna_grupo, columnas_mirna, grupo_calibrador):
    """ΔΔCt usando media de ΔCt por grupo; columnas_mirna son nombres dCt_."""
    _validar_grupo(df, columna_grupo)
    if grupo_calibrador not in df[columna_grupo].dropna().values:
        raise ValueError("El grupo calibrador no existe.")
    resultados = []
    for columna in columnas_mirna:
        datos = df[[columna_grupo]].copy()
        datos["delta"] = _numeros(df, columna)
        calibrador = datos.loc[datos[columna_grupo] == grupo_calibrador, "delta"].mean()
        if not np.isfinite(calibrador):
            raise ValueError(f"No hay ΔCt válidos en el calibrador para {columna}.")
        for grupo, sub in datos.groupby(columna_grupo, dropna=True, observed=True):
            media = sub["delta"].mean()
            ddct = media - calibrador
            resultados.append({"Grupo": grupo, "Calibrador": grupo_calibrador,
                               "miRNA": columna.removeprefix("dCt_"),
                               "n": int(sub['delta'].notna().sum()),
                               "Media_DeltaCt": media, "Media_DeltaCt_calibrador": calibrador,
                               "DeltaDeltaCt": ddct, "Fold_Change": 2.0 ** (-ddct),
                               "Log2_FC": -ddct})
    return pd.DataFrame(resultados)


def analizar_estadistica(df, columna_grupo, columnas_mirna, metodo="no_parametrico"):
    """Comparación de grupos independientes en escala ΔCt. Holm para pares, BH entre miRNAs."""
    _validar_grupo(df, columna_grupo)
    if metodo not in ("no_parametrico", "parametrico"):
        raise ValueError("Método no reconocido.")
    globales, pares_resultados = [], []
    for columna in columnas_mirna:
        datos = pd.DataFrame({"grupo": df[columna_grupo], "ct": _numeros(df, columna)}).dropna()
        grupos = [(nombre, sub["ct"].to_numpy()) for nombre, sub in datos.groupby("grupo", observed=True)]
        grupos = [(nombre, vals) for nombre, vals in grupos if len(vals) > 0]
        if len(grupos) < 2:
            continue
        nombres, valores = zip(*grupos)
        if any(len(v) < 2 for v in valores):
            continue
        try:
            if len(valores) == 2:
                if metodo == "parametrico":
                    prueba, p = "Welch", stats.ttest_ind(*valores, equal_var=False).pvalue
                else:
                    prueba, p = "Mann-Whitney U", stats.mannwhitneyu(*valores, alternative="two-sided").pvalue
            else:
                if metodo == "parametrico":
                    prueba, p = "ANOVA de Welch", stats.f_oneway(*valores, equal_var=False).pvalue
                else:
                    prueba, p = "Kruskal-Wallis", stats.kruskal(*valores).pvalue
        except (ValueError, FloatingPointError):
            continue
        if not np.isfinite(p):
            continue
        globales.append({"miRNA": columna.removeprefix("dCt_"), "Columna_DeltaCt": columna,
                         "Prueba": prueba, "p_global": float(p), "n_grupos": len(valores)})


        temp = []
        for i, j in combinations(range(len(valores)), 2):
            try:
                if metodo == "parametrico":
                    pp = stats.ttest_ind(valores[i], valores[j], equal_var=False).pvalue
                    nombre_prueba = "Welch"
                else:
                    pp = stats.mannwhitneyu(valores[i], valores[j], alternative="two-sided").pvalue
                    nombre_prueba = "Mann-Whitney U"
            except ValueError:
                pp = np.nan
            temp.append({"miRNA": columna.removeprefix("dCt_"), "Grupo_1": nombres[i],
                         "Grupo_2": nombres[j], "Prueba": nombre_prueba,
                         "p_sin_ajustar": float(pp)})
        validos = [i for i, r in enumerate(temp) if np.isfinite(r["p_sin_ajustar"])]
        if validos:
            ajustados = multipletests([temp[i]["p_sin_ajustar"] for i in validos], method="holm")[1]
            for i, ajustado in zip(validos, ajustados):
                temp[i]["p_ajustado_holm"] = float(ajustado)
        pares_resultados.extend(temp)
    if globales:
        global_p = [r["p_global"] for r in globales]
        for r, q in zip(globales, multipletests(global_p, method="fdr_bh")[1]):
            r["q_BH"] = float(q)
    return pd.DataFrame(globales), pd.DataFrame(pares_resultados)


def correlacion_spearman(df, columna_feature, columna_expresion):
    datos = pd.DataFrame({"feature": _numeros(df, columna_feature),
                          "expresion": _numeros(df, columna_expresion)}).dropna()
    if len(datos) < 3 or datos.nunique().min() < 2:
        raise ValueError("Se requieren al menos 3 pares válidos y variables no constantes.")
    resultado = stats.spearmanr(datos["feature"], datos["expresion"])
    return {"Feature": columna_feature, "Expresion": columna_expresion,
            "n": len(datos), "rho": float(resultado.statistic), "p": float(resultado.pvalue)}
