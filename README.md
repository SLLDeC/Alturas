# Alturas
Experimento de calibración temporal de las perturbaciones espaciales.

## Experimento

**Alturas_MANAGER.m**: experimento, genera sujetos_alturas.mat (var. sujeto)

## Procesado

- **todos_los_trialsXsujeto.m**: grafica los trials-baseline y promedios por condición de un sujeto en particular
- **Promedia_dentro_sujeto**: calcula las series temporales promedio de cada condición de cada sujeto (var. sujeto, condicion)
- **Promedia_entre_sujeto**: calcula las series temporales de cada condicion entre todos los sujetos, genera promedios.mat (var. condicion)

## Fiteos

- **Rodrigo.m**: distintos fiteos de los resultados de ALTURAS.
- **Fiteo_mars.m**: fitea por sujeto usando MARS
- **Fiteo_por_sujeto.m**: hace el fiteo por sujeto y genera el .mat adecuado.
- **data_fit.mat**: lista completa de trials del experimento ALTURAS para fitear.
- **datos_para_fiteo_por_sujeto.mat**: lista completa de trials del experimento ALTURAS para fitear por sujeto.

