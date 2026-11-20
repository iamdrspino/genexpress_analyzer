# Módulo de visualización


import numpy as np
import pandas as pd
import plotly.express as px


def _asteriscos(p):
    if not np.isfinite(p):
        return None
    if p < 0.0001:
        return '****'
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return 'ns'


def _anotar_comparaciones(fig, datos, grupo, expresion, tabla_pares,
                          mostrar_ns=False, columna_p='p_ajustado_holm',
                          escala_y='lineal'):
    if tabla_pares is None or tabla_pares.empty:
        return fig
    requeridas = {'Grupo_1', 'Grupo_2', columna_p}
    if not requeridas.issubset(tabla_pares.columns):
        raise ValueError('La tabla de comparaciones no contiene las columnas necesarias.')

    niveles = list(pd.unique(datos[grupo]))
    posiciones = {str(v): i for i, v in enumerate(niveles)}
    pares = []
    for _, fila in tabla_pares.iterrows():
        a, b = str(fila['Grupo_1']), str(fila['Grupo_2'])
        if a not in posiciones or b not in posiciones or a == b:
            continue
        p = pd.to_numeric(pd.Series([fila[columna_p]]), errors='coerce').iloc[0]
        if not np.isfinite(p):
            continue
        etiqueta = _asteriscos(float(p))
        if etiqueta == 'ns' and not mostrar_ns:
            continue
        i, j = sorted((posiciones[a], posiciones[b]))
        pares.append((i, j, etiqueta))

    if not pares:
        return fig

    valores = datos[expresion].to_numpy(dtype=float)
    if escala_y == 'log10':
        transformados = np.log10(valores)
    else:
        transformados = valores
    minimo = float(np.min(transformados))
    maximo = float(np.max(transformados))
    rango = max(maximo - minimo, 0.15 if escala_y == 'log10' else max(abs(maximo) * 0.1, 0.01))
    paso = rango * 0.15
    base = maximo + paso * 0.6
    altura = paso * 0.2

    # Ubicar intervalos cortos primero y evitar solapamientos horizontales.
    pares.sort(key=lambda t: (t[1] - t[0], t[0]))
    niveles_ocupados = []
    for i, j, etiqueta in pares:
        nivel = 0
        while any(nivel == nivel_prev and not (j < izq or i > der)
                  for izq, der, nivel_prev in niveles_ocupados):
            nivel += 1
        niveles_ocupados.append((i, j, nivel))
        y = base + nivel * paso
        def coord(v):
            return 10.0 ** v if escala_y == 'log10' else v
        fig.add_shape(type='path',
                      path=f'M {i},{coord(y-altura)} L {i},{coord(y)} L {j},{coord(y)} L {j},{coord(y-altura)}',
                      xref='x', yref='y', line=dict(color='#333333', width=1.5))
        fig.add_annotation(x=(i+j)/2, y=coord(y + altura * 0.3),
                           text=etiqueta, showarrow=False,
                           xref='x', yref='y', font=dict(size=15, color='#222222'),
                           yanchor='bottom')

    superior = base + (max(n for _, _, n in niveles_ocupados) + 1) * paso + altura * 2
    inferior = minimo - rango * 0.08
    if escala_y == 'log10':
        fig.update_yaxes(range=[inferior, superior])
    else:
        fig.update_yaxes(range=[inferior, superior])
    return fig


def generar_grafica(df, columna_grupo, mirna, tipo='boxplot',
                    escala_y='lineal', comparaciones=None,
                    mostrar_significancia=False, mostrar_ns=False):
#Grafica una columna 2^-dCt ya calculada, sin transformarla de nuevo.
# comparaciones: DataFrame filtrado al miRNA seleccionado, con Grupo_1,
#  Grupo_2 y p_ajustado_holm. No se calculan pruebas aquí.
 
    if columna_grupo not in df.columns or mirna not in df.columns:
        raise ValueError('No se encontró la Feature o la columna de expresión.')
    if escala_y not in ('lineal', 'log10'):
        raise ValueError("escala_y debe ser 'lineal' o 'log10'.")
    datos = df[[columna_grupo, mirna]].copy()
    datos[mirna] = pd.to_numeric(datos[mirna], errors='coerce').replace([np.inf, -np.inf], np.nan)
    datos = datos.dropna()
    if escala_y == 'log10':
        datos = datos.loc[datos[mirna] > 0].copy()
    if datos.empty:
        raise ValueError('No hay valores válidos para la escala seleccionada.')
    # Convertir etiquetas a texto para hacer coincidir las categorías con las pruebas.
    datos[columna_grupo] = datos[columna_grupo].astype(str)
    orden = list(pd.unique(datos[columna_grupo]))
    opciones = dict(data_frame=datos, x=columna_grupo, y=mirna,
                    category_orders={columna_grupo: orden},
                    title=f'Expresión normalizada: {mirna.removeprefix("2^-dCt_")}')
    if tipo == 'boxplot':
        fig = px.box(**opciones, points='all')
    elif tipo == 'violin':
        fig = px.violin(**opciones, box=True, points='all')
    elif tipo == 'barras':
        resumen = (datos.groupby(columna_grupo, observed=True, sort=False)[mirna]
                   .agg(['mean', 'std']).reset_index())
        # Las barras de error que cruzan cero no son representables en escala log.
        if escala_y == 'log10':
            resumen['std'] = np.minimum(resumen['std'].fillna(0), resumen['mean'] * 0.999)
        fig = px.bar(resumen, x=columna_grupo, y='mean', error_y='std',
                     category_orders={columna_grupo: orden},
                     title=f'Media ± DE: {mirna.removeprefix("2^-dCt_")}')
    else:
        raise ValueError('Tipo de gráfica no válido.')
    fig.update_layout(template='plotly_white', xaxis_title=columna_grupo,
                      yaxis_title='Expresión normalizada (2^-ΔCt)', showlegend=False)
    fig.update_xaxes(type='category', categoryorder='array', categoryarray=orden)
    fig.update_yaxes(type='log' if escala_y == 'log10' else 'linear')
    if mostrar_significancia:
        if comparaciones is None:
            raise ValueError('Ejecuta primero el análisis estadístico para mostrar significancia.')
        _anotar_comparaciones(fig, datos, columna_grupo, mirna, comparaciones,
                              mostrar_ns=mostrar_ns, escala_y=escala_y)
    return fig


def generar_dispersion(df, columna_feature, columna_expresion, escala_y='lineal'):
    if columna_feature not in df.columns or columna_expresion not in df.columns:
        raise ValueError('No se encontraron las columnas seleccionadas.')
    if escala_y not in ('lineal', 'log10'):
        raise ValueError('Escala no reconocida.')
    datos = df[[columna_feature, columna_expresion]].copy()
    for columna in (columna_feature, columna_expresion):
        datos[columna] = pd.to_numeric(datos[columna], errors='coerce')
    datos = datos.replace([np.inf, -np.inf], np.nan).dropna()
    if escala_y == 'log10':
        datos = datos.loc[datos[columna_expresion] > 0]
    if len(datos) < 2:
        raise ValueError('No hay suficientes pares válidos para graficar.')
    fig = px.scatter(datos, x=columna_feature, y=columna_expresion,
                     title='Asociación entre Feature y expresión')
    fig.update_layout(template='plotly_white')
    fig.update_yaxes(type='log' if escala_y == 'log10' else 'linear')
    return fig
