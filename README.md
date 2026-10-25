# Simulación y análisis de señales con la transformada de Fourier

## Propósito

Este proyecto representa un pulso rectangular, una función escalón y una senoide en tiempo discreto. Calcula y grafica sus espectros bilaterales de magnitud y fase mediante la FFT de NumPy. También compara las propiedades de linealidad, desplazamiento temporal y escalamiento de frecuencia usando una señal gaussiana de prueba.

## Requisitos

- Python 3.9 o posterior
- NumPy
- Matplotlib

## Ejecución

```bash
python -m pip install numpy matplotlib
python simulacion_fourier.py
```

El programa abre las figuras y guarda cuatro imágenes PNG en `figuras/`. Las imágenes pueden volver a generarse al ejecutar el programa.

## Contenido

- `simulacion_fourier.py`: definición de señales, FFT, gráficos y cálculos de propiedades.
- `figuras/01_dominio_tiempo.png`: señales muestreadas.
- `figuras/02_espectros.png`: magnitud y fase para cada señal.
- `figuras/03_propiedades.png`: comprobaciones numéricas de propiedades.
- `figuras/04_tiempo_frecuencia.png`: comparación de duración y ancho espectral.
- `Presentacion_Fourier.pptx`: material de apoyo para una exposición breve.

## Cómo interpretar los resultados

La FFT representa un bloque finito de muestras y supone que ese bloque se repite periódicamente. El escalón, aunque se define aquí sobre una ventana finita, no es una señal periódica infinita; su espectro aparece afectado por la ventana. La magnitud indica cuánto aporta cada frecuencia y la fase indica su relación temporal. Cerca de magnitudes nulas, la fase no es estable y se omite en la gráfica.

Un pulso rectangular produce un espectro tipo sinc: al reducir su duración temporal, aumenta el ancho del espectro. Una senoide de 50 Hz muestra picos alrededor de ±50 Hz en el espectro bilateral. Para señales reales, magnitud y fase tienen simetría conjugada.

Las propiedades ideales son: linealidad, `x(t−t₀) ↔ X(f)e^(−j2πft₀)` y escalamiento `x(at) ↔ (1/|a|)X(f/a)`. En la computadora se comparan muestras finitas; por ello los errores pequeños dependen de la resolución, la ventana y el redondeo del retardo.

## Publicación en GitHub

Crea un repositorio vacío en GitHub, copia estos archivos y ejecuta en esta carpeta:

```bash
git init
git add .
git commit -m "Agregar simulación y análisis de Fourier"
git branch -M main
git remote add origin https://github.com/USUARIO/REPOSITORIO.git
git push -u origin main
```

Reemplaza `USUARIO/REPOSITORIO` con la dirección real de tu repositorio. No se ha creado ni publicado un repositorio desde este paquete.

## Presentación oral

Personaliza el campo `[Nombre del estudiante]` en la portada. Explica primero las muestras en tiempo, después la magnitud y la fase; termina con las propiedades y la limitación de analizar un escalón mediante una ventana finita.
