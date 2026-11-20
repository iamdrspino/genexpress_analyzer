#Genexpress Analyzer 
# Version 1.1 
#Aplicación para el análisis de expresión diferencial por qPCR
#cDr. Israel Martínez-Espinosa
# Programa de Doctorado en Ciencias Biomédicas-UNAM 2027-I
# Programación y modelado computacional usando Python

##### Módulos --------------------------------------------------------------####

import pandas as pd    # Manejo de datos
import streamlit as st # Interfaz gráfica de la aplicación
import hashlib         # Asignar hash al conjunto de datos

#### Módulos exclusivos de Genexpress Analyzer ----------------------------#####

#1. importacion.py, adquiere los datos desde un archivo externo
from modules.importacion import (
  cargar_archivo,
  validar_datos
  )
  
#2. validacion.py, regresa un informe de la estructura de los datos
from modules.validacion import (
  validar_integridad,
  generar_reporte_calidad)
  
#3. identificacion.py, identifica el tipo de dato de cada columna
from modules.identificacion import (
  obtener_columnas,
  identificar_variables
  )
  
#4. normalizacion.py, normaliza los datos de expresión con base en un control
from modules.normalizacion import (
  normalizar_qpcr,
  evaluar_referencia
  )

#5. analisis.py, análisis estadísticos de los valores normalizados
from modules.analisis import (
  calcular_expresion_relativa,
  calcular_fold_change,
  analizar_estadistica,
  correlacion_spearman
  )
  
#6. visualización.py, gráfica los resultados.  
from modules.visualizacion import (
  generar_grafica,
  generar_dispersion
  )
  
#7. exportacion.py, genera un reporte de resultados.
from modules.exportacion import (
  generar_html,
  generar_pdf
  )

#### Funciones ------------------------------------------------------------#####

def formatear_valor_pq(numero):
# Formato de números en valores de pruebas estadísticas
# Si < 0.0005 lo representa en notación científica
# Si > 0.0005 utiliza cuatro cifras decimales

    if pd.isna(numero):
        return "NA"

    if numero < 0.0005:
        return f"{numero:.3e}"

    return f"{numero:.4f}"


def formatear_tabla_estadistica(df):
# Dar formato únicamente a las columnas correspondientes a valores p y q.

    tabla = df.copy()

    for columna in tabla.columns:
        nombre = columna.lower()
        if (nombre.startswith("p_")
            or nombre.startswith("q_")
            or nombre in ["p", "q", "pvalue", "qvalue"]):

            tabla[columna] = tabla[columna].apply(formatear_valor_pq)
    return tabla

#### Configuración de Streamlit -------------------------------------------#####

st.set_page_config(
    page_title="Genexpress Analyzer v 1.1 por Israel MArtínez, 2026",
    page_icon="🧬",
    layout="wide"
    )

#### Encabezado de Genexpress Analyzer ------------------------------------#####

st.title("🧬 Genexpress Analyzer v.1.1", text_alignment="center")

st.write("Plataforma para el análisis de expresión a partir de datos de PCR en tiempo real.")
st.write("Los datos para realizar el análisis deben presentarse en una tabla en formato .csv, .xls o .xlsx.")
st.write('''Por el momento, solo puedo analizar un análisis por muestra. El desarrollo de esta herramienta
 seguirá hasta poder realizar triplicados''')

st.divider()

#### 1. Importación de datos ---------------------------------------------######

st.header("1. Importación de datos")
st.subheader("Carga de archivos")
st.write("Formatos aceptados .csv, .xls y .xlsx.")

archivo = st.file_uploader(
    "Archivo de datos",
    type=["csv", "xlsx", "xls"]
)

#### 2. Validación de los datos ------------------------------------------######

if archivo is not None:
    try:
        df = cargar_archivo(archivo)
        st.success("Archivo importado correctamente.")
        st.subheader("Vista previa de los datos")
        st.write("Se muestran las 5 filas superiores e inferiores.")
        st.dataframe(
            df.head(6),
            use_container_width=True
        )
        st.dataframe(
            df.tail(6),
            use_container_width=True
        )

        # Validación de los datos

        resultados = validar_datos(df)
        st.subheader("Resumen de los datos cargados.")

        c1, c2, c3 = st.columns(3)

        c1.metric(
            "Muestras",
            int(resultados["filas"])
        )
        c2.metric(
            "Varibles",
            int(resultados["numero_columnas"])
        )

        c3.metric(
            "Valores faltantes",
            int(resultados["valores_faltantes"])
        )

        st.subheader("Evaluación de los datos")

        if resultados["filas_vacias"] > 0:
            st.warning(
                f"Se encontraron {resultados['filas_vacias']} filas vacías."
            )

        if resultados["columnas_vacias"] > 0:
            st.warning(
                f"Se encontraron {resultados['columnas_vacias']} columnas vacías."
            )

        if resultados["columnas_duplicadas"] > 0:
            st.warning(
                "Se detectaron nombres de columnas duplicados."
            )
        else:
          st.success("No se encontraron filas o columnas vacías.")
          st.success("No se encontraron filas o columnas duplicadas.")

        #### 3. Identificación de variables-------------------------------######
        st.divider()
        st.header("2. Identificación de variables")
        st.write("Clasifica manualmente cada columna según el tipo de información que contiene.")
        
        columnas = obtener_columnas(df)
        
# Con la librería hashlib creamos un hash que permite identificar nuestro conjunto
# de datos original para evitar reutilizar clasificaciones de otro conjunto de datos.

        firma_archivo = hashlib.sha256(
            archivo.getvalue()
        ).hexdigest()
        
        firma = (
            firma_archivo,
            tuple(columnas)
            )

        if st.session_state.get("firma_datos") != firma:
        
            st.session_state["firma_datos"] = firma
            st.session_state["clasificacion"] = {}
            st.session_state["variables"] = None
            st.session_state["df_normalizado"] = None
            st.session_state["resumen_norm"] = None

# Garantizar que la identificación del tipo de variables no tenga errores.        
            for clave in list(st.session_state):
                if clave.startswith("tipo_var_"):
                    del st.session_state[clave]
        
                if clave.startswith("subtipo_var_"):
                    del st.session_state[clave]
        
        
#### 3.1. Clasificación manual de variables ------------------------------######
        
        clasificacion = {}
        
        for i, columna in enumerate(columnas):
        
            st.markdown(f"{columna}")
        
            col1, col2 = st.columns(2)
        
            with col1:
        
                categoria = st.selectbox(
                    "Categoría",
                    [
                        "Sin clasificar",
                        "id",
                        "caracteristica",
                        "control",
                        "expresion"
                    ],
                    format_func=lambda x: {
                        "Sin clasificar": "Sin clasificar",
                        "id": "ID / Nombre",
                        "caracteristica": "Característica",
                        "control": "Control",
                        "expresion": "Expresión"
                    }[x],
                    key=f"tipo_var_{i}"
                )
        
            subtipo = None
        
            with col2:
        
                if categoria == "caracteristica":
        
                    subtipo = st.selectbox(
                        "Tipo de característica",
                        [
                            "categorica",
                            "numerica"
                        ],
                        format_func=lambda x: {
                            "categorica": "Categórica",
                            "numerica": "Numérica"
                        }[x],
                        key=f"subtipo_var_{i}_feature"
                    )
        
                elif categoria == "control":
        
                    subtipo = st.selectbox(
                        "Tipo de control",
                        [
                            "endogeno",
                            "exogeno"
                        ],
                        format_func=lambda x: {
                            "endogeno": "Endógeno",
                            "exogeno": "Exógeno"
                        }[x],
                        key=f"subtipo_var_{i}_control"
                    )
        
            if categoria != "Sin clasificar":
        
                clasificacion[columna] = {
                    "categoria": categoria,
                    "subtipo": subtipo
                }
        
            st.divider()
        
        confirmar = st.button(
            "Confirmar clasificación",
            type="primary",
            key="confirmar_variables"
        )


        if confirmar:
        
            try:
        
                variables = identificar_variables(
                    df,
                    clasificacion
                )
        
                st.session_state["clasificacion"] = (
                    clasificacion.copy()
                )
        
                st.session_state["variables"] = variables
        
                # Eliminar normalizaciones anteriores
                st.session_state["df_normalizado"] = None
                st.session_state["resumen_norm"] = None
        
                st.success(
                    "Clasificación guardada correctamente."
                )
        
            except ValueError as error:
        
                st.error(str(error))
        
        
# Variables que se identificaron
        
        variables = st.session_state.get("variables")
        
        if variables is not None:
        
            st.subheader("Resumen de variables identificadas.")
        
            c1, c2, c3, c4 = st.columns(4)
        
            c1.metric("ID", len(variables["id"]))
            c2.metric("Características", len(variables["caracteristica"]))
            c3.metric("Control", len(variables["control"]))
            c4.metric("Expresión", len(variables["expresion"]))
        
            st.write("Clasificación detallada")
        
            tabla = []
        
            for columna, datos in (
                st.session_state["clasificacion"].items()
            ):
        
                tabla.append({
                    "Columna": columna,
                    "Categoría": datos["categoria"],
                    "Subtipo": datos["subtipo"] or "—"
                })
        
            st.dataframe(
                pd.DataFrame(tabla),
                use_container_width=True,
                hide_index=True
            )
        
        else:
            st.info(
                "Confirma la clasificación para continuar."
            )
#### 4. Validación de los datos ------------------------------------------######

        st.divider()
        st.header("3. Validación del conjunto de datos.")
        st.subheader("Rango de valores de CT como criterio de aceptación.")
        st.text("Nota: La aplicación no excluirá valores fuera del rango de CT, solo mostrará alertas.")

        # Intervalo de CT de los valores de expresión  para informar de datos fuera de rango.
      
        col1, col2 = st.columns(2)
        
        with col1:
            ct_min = st.number_input(
                "CT mínimo",
                value=0.0,
                step=1.0
            )
        
        with col2:
            ct_max = st.number_input(
                "CT máximo",
                value=40.0,
                step=1.0
            )
        
        if ct_min >= ct_max:
            st.error("El CT mínimo debe ser menor que el máximo.")
            st.stop()
        
        
        reporte = validar_integridad(
            df,
            variables,
            ct_min=ct_min,
            ct_max=ct_max
        )
        
        # Recuperar la clasificación confirmada
        variables = st.session_state.get("variables")
        
        if variables is None:
            st.warning(
                "Primero debes confirmar la clasificación de variables."
            )
            st.stop()
        
       
        
        st.subheader("3.1 Resumen general.")
        
        col1, col2, col3 = st.columns(3)
        
        col1.metric(
            "Total de muestras",
            reporte["total_muestras"]
        )
        
        col2.metric(
            "Total de columnas",
            reporte["total_columnas"]
        )
        
        col3.metric(
            "Valores faltantes",
            reporte["total_faltantes"]
        )
        

        st.subheader("3.2 Identificadores de muestras. ")
        
        datos_id = reporte["identificadores"]
        
        if datos_id:
        
            tabla_id = pd.DataFrame.from_dict(
                datos_id,
                orient="index"
            )
        
            tabla_id.index.name = "Columna ID"
        
            st.dataframe(
                tabla_id.reset_index(),
                use_container_width=True,
                hide_index=True
            )
        
        else:
            st.info("No se encontraron identificadores.")
        
        st.subheader("3.3 Características asociadas con las muestras.")
        
        caracteristicas = reporte["caracteristicas"]
        
        if caracteristicas:
        
            for columna, resultado in caracteristicas.items():
        
                with st.expander(f"Característica: {columna}"):
        
                    st.write(
                        "Tipo de característica:",
                        resultado["tipo"]
                    )
        
                    if resultado["tipo"] == "categorica":
        
                        st.metric(
                            "Número de categorías",
                            resultado["numero_categorias"]
                        )
        
                        st.metric(
                            "Valores faltantes",
                            resultado["faltantes"]
                        )
        
                        frecuencias = pd.DataFrame(
                            list(resultado["frecuencias"].items()),
                            columns=[
                                "Categoría",
                                "Frecuencia"
                            ]
                        )
        
                        st.dataframe(
                            frecuencias,
                            use_container_width=True,
                            hide_index=True
                        )
        
                    elif resultado["tipo"] == "numerica":
        
                        resumen = {
                            clave: valor
                            for clave, valor in resultado.items()
                            if clave != "tipo"
                        }
        
                        st.dataframe(
                            pd.DataFrame([resumen]),
                            use_container_width=True,
                            hide_index=True
                        )
        
        else:
            st.info("No se encontraron Características.")
        
        st.subheader("3.4 Evaluación de los datos de expresión de RT-PCR.")
        
        datos_ct = reporte["expresion"]
        
        if datos_ct:
        
            tabla_ct = pd.DataFrame.from_dict(
                datos_ct,
                orient="index"
            )
        
            tabla_ct.index.name = "Columna"
        
            st.dataframe(
                tabla_ct.reset_index(),
                use_container_width=True,
                hide_index=True
            )
        
        else:
            st.info(
                "No se encontraron columnas de Control o Expresión."
            )
        
        st.subheader("3.5 Advertencias.")
        
        if reporte["advertencias"]:
        
            for advertencia in reporte["advertencias"]:
                st.warning(advertencia)
        
        else:
            st.success(
                "No se detectaron problemas con el conjunto de datos."
            )
        

        st.info(
            "Validación finalizada. El conjunto de datos no ha sido modificado."
        )
        
#### 5. Normalización de los datos de expresión---------------------------######

        st.divider()
        st.header("4. Normalización de expresión")

        variables = st.session_state.get("variables")

        if variables is None:
            st.warning(
                "Confirma primero la clasificación."
            )
            st.stop()

        controles = variables.get("control", [])
        mirnas = variables.get("expresion", [])

        if not controles or not mirnas:
            st.warning(
                "Se necesita al menos una columna "
                "Control y una de Expresión."
            )
            st.stop()

        st.subheader("4.1 Selección de transcrito control de  referencia.")

        referencia = st.selectbox(
            "Control de referencia",
            controles,
            key="referencia_normalizacion"
        )

        st.subheader("4.2 Selección de transcritos por analizar.")

        seleccionados = st.multiselect(
            "Transcritos a normalizar",
            mirnas,
            default=mirnas,
            key="mirnas_normalizacion"
        )

        st.subheader("4.3 Evaluación del transcrito control de referencia.")

        evaluacion = evaluar_referencia(
            df,
            referencia
        )

        st.dataframe(
            pd.DataFrame([evaluacion]),
            hide_index=True,
            use_container_width=True
        )

        st.subheader("4.4 Cálculo de ΔCT. ")

        if st.button(
            "Ejecutar normalización",
            type="primary",
            key="ejecutar_normalizacion"
        ):

            if not seleccionados:
                st.warning(
                    "Selecciona al menos un transcrito."
                )

            else:
                try:
                    df_normalizado, resumen_norm = (
                        normalizar_qpcr(
                            df,
                            variables,
                            columna_referencia=referencia,
                            mirnas_seleccionados=seleccionados
                        )
                    )

                    st.session_state[
                        "df_normalizado"
                    ] = df_normalizado

                    st.session_state[
                        "resumen_norm"
                    ] = resumen_norm

                    st.session_state[
                        "config_normalizacion"
                    ] = (
                        referencia,
                        tuple(seleccionados)
                    )

                    st.success(
                        "Normalización realizada correctamente."
                    )

                except Exception as error:
                    st.error(
                        f"Error en normalización: {error}"
                    )

# Resultados de normalización

        df_norm = st.session_state.get(
            "df_normalizado"
        )

        resumen = st.session_state.get(
            "resumen_norm"
        )

        configuracion_actual = (
            referencia,
            tuple(seleccionados)
        )

        configuracion_guardada = (
            st.session_state.get(
                "config_normalizacion"
            )
        )

        if (
            df_norm is not None
            and resumen is not None
            and configuracion_actual
            == configuracion_guardada
        ):

            st.subheader(
                "4.5 Resultados de normalización"
            )

            st.write(
                "Transcrito de referencia:",
                resumen["referencia"]
            )

            st.dataframe(
                df_norm,
                use_container_width=True,
                hide_index=True
            )

            tabla_resumen = pd.DataFrame.from_dict(
                resumen["resultados"],
                orient="index"
            ).reset_index()

            tabla_resumen = tabla_resumen.rename(
                columns={"index": "microRNA"}
            )

            st.subheader(
                "Resumen de muestras normalizadas"
            )

            st.dataframe(
                tabla_resumen,
                hide_index=True,
                use_container_width=True
            )

            csv = df_norm.to_csv(
                index=False
            ).encode("utf-8-sig")

            st.download_button(
                "Descargar matriz normalizada",
                data=csv,
                file_name="mirna_normalizado.csv",
                mime="text/csv"
            )

        elif df_norm is not None:
            st.info(
                "La selección cambió. Ejecuta nuevamente "
                "la normalización para actualizar los resultados."
            )
          
#### 6. Análisis bioestadístico _______________---------------------------######


        st.divider()
        st.header("5. Análisis bioestadístico")

        df_norm = st.session_state.get("df_normalizado")
        resumen_norm = st.session_state.get("resumen_norm")

        if df_norm is None or resumen_norm is None:
            st.info("Ejecuta primero la normalización.")
            st.stop()

        features = variables.get(
            "caracteristica_categorica", []
        )

        columnas_dct = [
            datos["columna_dct"]
            for datos in resumen_norm["resultados"].values()
        ]

        if not features:
            st.warning(
                "No hay características categóricas para comparar."
            )
            st.stop()

        grupo = st.selectbox(
            "Característica para comparar grupos",
            features,
            key="grupo_analisis"
        )

        calibradores = df_norm[grupo].dropna().unique().tolist()

        if not calibradores:
            st.warning("La característica no contiene grupos válidos.")
            st.stop()

        calibrador = st.selectbox(
            "Grupo calibrador",
            calibradores,
            key="calibrador_analisis"
        )

        metodo = st.radio(
            "Enfoque estadístico",
            ["no_parametrico", "parametrico"],
            format_func=lambda x: {
                "no_parametrico": "No paramétrico",
                "parametrico": "Paramétrico"
            }[x]
        )

        if st.button("Ejecutar análisis estadístico"):

            tabla_fc = calcular_fold_change(
                df_norm,
                grupo,
                columnas_dct,
                calibrador
            )

            tabla_global, tabla_pares = analizar_estadistica(
                df_norm,
                grupo,
                columnas_dct,
                metodo=metodo
            )
            
            st.session_state["tabla_global"] = tabla_global
            st.session_state["tabla_pares"] = tabla_pares
            
            st.session_state["config_estadistica"] = (
                grupo,
                metodo,
                tuple(columnas_dct)
            )

            st.subheader("Fold change")
            st.dataframe(tabla_fc, hide_index=True)
            st.session_state["tabla_fc"] = tabla_fc

            st.subheader("Pruebas globales")

            st.dataframe(
                formatear_tabla_estadistica(tabla_global),
                hide_index=True,
                use_container_width=True
            )
            
            st.subheader("Comparaciones por pares")
            
            st.dataframe(
                formatear_tabla_estadistica(tabla_pares),
                hide_index=True,
                use_container_width=True
            )


#### 7. Visualización de datos en gráficas -------------------------------######
        
        st.divider()
        st.header("6. Visualización")

        columnas_exp = [
            datos["columna_expresion"]
            for datos in resumen_norm["resultados"].values()
        ]

        mirna_grafica = st.selectbox(
            "Transcrito para graficar",
            columnas_exp
        )

        tipo_grafica = st.selectbox(
            "Tipo de gráfica",
            ["boxplot", "violin", "barras"]
        )

        escala_y = st.radio(
            "Escala del eje Y",
            ["lineal", "log10"],
            format_func=lambda x: {
                "lineal": "Decimal (lineal)",
                "log10": "Logarítmica (base 10)"
            }[x],
            horizontal=True
        )

        mostrar_significancia = st.checkbox(
            "Mostrar significancia estadística",
            value=False
        )

        mostrar_ns = False

        if mostrar_significancia:
            mostrar_ns = st.checkbox(
                "Incluir comparaciones no significativas (ns)",
                value=False
            )

        # Recuperar resultados estadísticos
        tabla_pares = st.session_state.get(
            "tabla_pares"
        )

        config_guardada = st.session_state.get(
            "config_estadistica"
        )

        config_actual = (
            grupo,
            metodo,
            tuple(columnas_dct)
        )

        # Identificar miRNA seleccionado
        nombre_mirna = mirna_grafica.removeprefix(
            "2^-ΔCT"
        )

        comparaciones = None

        if mostrar_significancia:

            if (
                tabla_pares is None
                or config_guardada != config_actual
            ):
                st.warning(
                    "Ejecuta el análisis estadístico "
                    "con la configuración actual antes "
                    "de mostrar la significancia."
                )

                mostrar_significancia = False

            else:
                comparaciones = tabla_pares.loc[
                    tabla_pares["miRNA"] == nombre_mirna
                ].copy()

                if comparaciones.empty:
                    st.info(
                        "No hay comparaciones estadísticas "
                        "disponibles para este transcrito."
                    )

        try:
            figura = generar_grafica(
                df_norm,
                grupo,
                mirna_grafica,
                tipo=tipo_grafica,
                escala_y=escala_y,
                comparaciones=comparaciones,
                mostrar_significancia=mostrar_significancia,
                mostrar_ns=mostrar_ns
            )

            st.plotly_chart(
                figura,
                use_container_width=True
            )
            st.session_state["figura_exportacion"] = figura

        except ValueError as error:
            st.warning(str(error))
        
#### 8. Exportación de resultados  ---------------------------------------######

        st.divider()
        st.header("7. Exportación de resultados")

        st.write(
            "Genera un informe con los resultados del análisis de expresión diferencial."
        )

        formato = st.radio(
            "Formato del informe",
            ["HTML", "PDF"],
            horizontal=True,
            key="formato_exportacion"
        )

        incluir_graficos = st.checkbox(
            "Incluir gráfica actual",
            value=True,
            key="incluir_graficos_exportacion"
        )

        titulo_informe = st.text_input(
            "Título del informe",
            value="Informe de expresión de RT-PCR"
        )

        # Recuperar resultados
        datos_exportacion = {
            "df_normalizado": st.session_state.get(
                "df_normalizado"
            ),
            "resumen_norm": st.session_state.get(
                "resumen_norm"
            ),
            "tabla_fc": st.session_state.get(
                "tabla_fc"
            ),
            "tabla_global": st.session_state.get(
                "tabla_global"
            ),
            "tabla_pares": st.session_state.get(
                "tabla_pares"
            ),
            "grupo": grupo,
            "metodo": metodo
        }

        figuras = []

        if incluir_graficos:
            figura_guardada = st.session_state.get(
                "figura_exportacion"
            )

            if figura_guardada is not None:
                figuras.append(figura_guardada)

        # Generar informe
        if st.button(
            "Preparar informe",
            type="primary",
            key="preparar_informe"
        ):

            try:
                if formato == "HTML":

                    archivo = generar_html(
                        datos_exportacion,
                        figuras=figuras,
                        titulo=titulo_informe
                    )

                    st.session_state[
                        "archivo_exportacion"
                    ] = archivo

                    st.session_state[
                        "extension_exportacion"
                    ] = "html"

                else:

                    archivo = generar_pdf(
                        datos_exportacion,
                        figuras=figuras,
                        titulo=titulo_informe
                    )

                    st.session_state[
                        "archivo_exportacion"
                    ] = archivo

                    st.session_state[
                        "extension_exportacion"
                    ] = "pdf"

                st.success(
                    "Informe generado correctamente."
                )

            except Exception as error:
                st.error(
                    f"Error al generar el informe: {error}"
                )

        # Descargar informe
        archivo = st.session_state.get(
            "archivo_exportacion"
        )

        extension = st.session_state.get(
            "extension_exportacion"
        )

        if archivo is not None and extension is not None:

            st.download_button(
                label="Descargar informe",
                data=archivo,
                file_name=(
                    f"analisis_de_expresion.{extension}"
                ),
                mime=(
                    "text/html"
                    if extension == "html"
                    else "application/pdf"
                ),
                key="descargar_informe"
            )


    except Exception as error:
        st.error(f"Error durante la ejecución: {error}")
        
       
