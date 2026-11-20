"""Exportación de informes de miRNA-Lab Analyzer a HTML y PDF."""
from __future__ import annotations

from datetime import datetime
from html import escape
from io import BytesIO
import math
import pandas as pd


def _tabla_html(tabla, max_filas=None):
    if tabla is None or not isinstance(tabla, pd.DataFrame) or tabla.empty:
        return '<p class="vacio">Sin resultados disponibles.</p>'
    vista = tabla if max_filas is None else tabla.head(max_filas)
    return vista.to_html(index=False, classes="tabla", border=0, escape=True, na_rep="NA")


def _secciones(datos):
    return [
        ("Matriz normalizada", datos.get("df_normalizado")),
        ("Fold change (ΔΔCt)", datos.get("tabla_fc")),
        ("Pruebas estadísticas globales", datos.get("tabla_global")),
        ("Comparaciones por pares", datos.get("tabla_pares")),
    ]


def generar_html(datos, figuras=None, titulo="Informe de miRNA-Lab Analyzer"):
    """Devuelve un HTML autocontenido en bytes UTF-8.

    figuras: lista de figuras Plotly (requiere plotly instalado).
    """
    figuras = figuras or []
    partes = [f'''<!doctype html><html lang="es"><head><meta charset="utf-8">
<title>{escape(titulo)}</title><style>
body{{font-family:Arial,sans-serif;color:#243246;margin:32px auto;max-width:1200px;padding:0 20px}}
h1{{color:#193957}}h2{{color:#28577a;margin-top:32px;border-bottom:1px solid #dbe3ea;padding-bottom:8px}}
p{{line-height:1.5}}.meta{{background:#f1f5f8;padding:16px;border-radius:8px}}
.tabla{{border-collapse:collapse;width:100%;font-size:12px}}
.tabla th,.tabla td{{border:1px solid #d9e1e8;padding:7px;text-align:left}}
.tabla th{{background:#eaf0f5;position:sticky;top:0}}.scroll{{overflow-x:auto;max-height:580px;overflow-y:auto}}
.vacio{{color:#777;font-style:italic}}.nota{{font-size:12px;color:#526276}}
@media print{{.scroll{{max-height:none;overflow:visible}}body{{margin:10mm}}}}
</style></head><body><h1>{escape(titulo)}</h1>
<div class="meta"><b>Generado:</b> {escape(datetime.now().strftime('%Y-%m-%d %H:%M'))}<br>
<b>Referencia:</b> {escape(str((datos.get('resumen_norm') or {}).get('referencia', 'No especificada')))}<br>
<b>Grupo de comparación:</b> {escape(str(datos.get('grupo', 'No especificado')))}<br>
<b>Método estadístico:</b> {escape(str(datos.get('metodo', 'No especificado')))}</div>
<p class="nota">Los valores de las tablas se exportan con su precisión numérica original. Las pruebas y los ajustes se calcularon previamente en la aplicación.</p>''']
    for nombre, tabla in _secciones(datos):
        partes.append(f'<h2>{escape(nombre)}</h2><div class="scroll">{_tabla_html(tabla)}</div>')
    if figuras:
        partes.append('<h2>Visualizaciones</h2>')
        for indice, figura in enumerate(figuras):
            partes.append(figura.to_html(full_html=False, include_plotlyjs=True if indice == 0 else False))
    partes.append('</body></html>')
    return '\n'.join(partes).encode('utf-8')


def generar_pdf(datos, figuras=None, titulo="Informe de miRNA-Lab Analyzer"):
    """Devuelve PDF en bytes. Para incluir figuras Plotly requiere kaleido.

    Si una figura no puede convertirse, se registra la omisión en el informe.
    """
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, LongTable, TableStyle, Image, KeepTogether

    def texto_seguro(x, max_chars=75):
        if pd.isna(x):
            return 'NA'
        if isinstance(x, (float, int)) and not isinstance(x, bool):
            if isinstance(x, float) and math.isfinite(x):
                return f'{x:.6g}'
        return str(x)[:max_chars]

    buffer = BytesIO()
    pagina = landscape(A4)
    doc = SimpleDocTemplate(buffer, pagesize=pagina, leftMargin=16*mm,
                            rightMargin=16*mm, topMargin=15*mm, bottomMargin=15*mm)
    estilos = getSampleStyleSheet()
    elementos = [Paragraph(escape(titulo), estilos['Title']), Spacer(1, 6*mm)]
    resumen = datos.get('resumen_norm') or {}
    for clave, valor in [('Fecha', datetime.now().strftime('%Y-%m-%d %H:%M')),
                         ('Referencia', resumen.get('referencia', 'No especificada')),
                         ('Grupo', datos.get('grupo', 'No especificado')),
                         ('Método', datos.get('metodo', 'No especificado'))]:
        elementos.append(Paragraph(f'<b>{escape(clave)}:</b> {escape(str(valor))}', estilos['Normal']))
    elementos.append(Spacer(1, 5*mm))
    elementos.append(Paragraph('Las tablas contienen resultados calculados previamente; no se recalculan pruebas.', estilos['Italic']))

    for nombre, tabla in _secciones(datos):
        elementos.extend([Spacer(1, 6*mm), Paragraph(escape(nombre), estilos['Heading2'])])
        if tabla is None or not isinstance(tabla, pd.DataFrame) or tabla.empty:
            elementos.append(Paragraph('Sin resultados disponibles.', estilos['Normal']))
            continue
        # Dividir tablas anchas en bloques para evitar recortes en PDF.
        columnas = list(tabla.columns)
        ancho_util = pagina[0] - 32*mm
        max_cols = 6
        for inicio in range(0, len(columnas), max_cols):
            bloque = columnas[inicio:inicio+max_cols]
            elementos.append(Paragraph(f'Columnas {inicio+1}–{inicio+len(bloque)} de {len(columnas)}', estilos['Italic']))
            filas = [[Paragraph(escape(str(c)), estilos['BodyText']) for c in bloque]]
            for fila in tabla[bloque].itertuples(index=False, name=None):
                filas.append([Paragraph(escape(texto_seguro(v)), estilos['BodyText']) for v in fila])
            t = LongTable(filas, colWidths=[ancho_util/len(bloque)]*len(bloque), repeatRows=1, hAlign='LEFT')
            t.setStyle(TableStyle([
                ('BACKGROUND', (0,0),(-1,0), colors.HexColor('#e8eff5')),
                ('GRID',(0,0),(-1,-1),0.3,colors.HexColor('#c8d2dc')),
                ('VALIGN',(0,0),(-1,-1),'TOP'),
                ('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),
                ('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4),
            ]))
            elementos.append(t)
    if figuras:
        elementos.extend([Spacer(1, 6*mm), Paragraph('Visualizaciones', estilos['Heading2'])])
        for i, fig in enumerate(figuras, start=1):
            try:
                png = fig.to_image(format='png', width=1100, height=650, scale=1)
                imagen = Image(BytesIO(png), width=245*mm, height=145*mm)
                elementos.append(KeepTogether([Paragraph(f'Figura {i}', estilos['Heading3']), imagen]))
            except Exception:
                elementos.append(Paragraph(f'Figura {i}: no se pudo incorporar (verifica Kaleido y sus dependencias).', estilos['Italic']))
    doc.build(elementos)
    return buffer.getvalue()
