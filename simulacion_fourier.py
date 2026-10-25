"""Simulación didáctica de señales y transformada discreta de Fourier.

Ejecutar: python simulacion_fourier.py
Genera cuatro PNG en la carpeta figuras/ y muestra los gráficos.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figuras"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"figure.figsize": (10, 5), "font.size": 11,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white", "axes.facecolor": "#fbfcfe"})
BLUE, ORANGE, GREEN, PURPLE = "#2864dc", "#ef7d32", "#179b78", "#7b61a8"

# Malla temporal con frecuencia de muestreo Fs y N muestras.
Fs, N = 1000.0, 1000
t = np.arange(N) / Fs - 0.5
freq = np.fft.fftshift(np.fft.fftfreq(N, d=1/Fs))

# Señales elementales: pulso rectangular centrado, escalón y senoide de 50 Hz.
pulso = (np.abs(t) <= 0.05).astype(float)
escalon = (t >= 0).astype(float)
seno = np.sin(2 * np.pi * 50 * t)
senales = [("Pulso rectangular", pulso, BLUE), ("Función escalón", escalon, ORANGE),
           ("Senoide de 50 Hz", seno, GREEN)]

# Espectro bilateral de la DFT: magnitud y fase centradas en frecuencia cero.
def espectro(x):
    X = np.fft.fftshift(np.fft.fft(x)) / N
    mag = np.abs(X)
    fase = np.angle(X)
    fase[mag < 1e-3 * mag.max()] = np.nan  # La fase no se interpreta donde la magnitud es casi cero.
    return X, mag, fase

fig, axs = plt.subplots(3, 1, figsize=(10, 7), sharex=True, constrained_layout=True)
for ax, (nombre, x, color) in zip(axs, senales):
    ax.plot(t, x, color=color, lw=2)
    ax.set_ylabel("Amplitud")
    ax.set_title(nombre, loc="left", fontsize=12, fontweight="bold")
    ax.grid(alpha=.22)
axs[-1].set_xlabel("Tiempo (s)")
fig.suptitle("Señales en el dominio del tiempo", fontsize=17, fontweight="bold")
fig.savefig(FIG / "01_dominio_tiempo.png", dpi=180, bbox_inches="tight")

fig, axs = plt.subplots(3, 2, figsize=(12, 9), constrained_layout=True)
for row, (nombre, x, color) in enumerate(senales):
    _, mag, fase = espectro(x)
    axs[row, 0].plot(freq, mag, color=color, lw=1.5)
    axs[row, 0].set_xlim(-150, 150)
    axs[row, 0].set_ylabel("|X[k]|")
    axs[row, 0].set_title(f"{nombre}: magnitud", loc="left")
    axs[row, 0].grid(alpha=.22)
    axs[row, 1].plot(freq, fase, color=PURPLE, lw=1.1)
    axs[row, 1].set_xlim(-150, 150)
    axs[row, 1].set_ylim(-np.pi, np.pi)
    axs[row, 1].set_yticks([-np.pi, 0, np.pi], ["−π", "0", "π"])
    axs[row, 1].set_title("Fase (rad)", loc="left")
    axs[row, 1].grid(alpha=.22)
axs[-1, 0].set_xlabel("Frecuencia (Hz)")
axs[-1, 1].set_xlabel("Frecuencia (Hz)")
fig.suptitle("Espectros bilaterales calculados con FFT", fontsize=17, fontweight="bold")
fig.savefig(FIG / "02_espectros.png", dpi=180, bbox_inches="tight")

# Propiedades en una señal gaussiana: los errores se evalúan en una malla finita.
dt = 0.001
tp = np.arange(-0.5, 0.5, dt)
a = 5.0
g = np.exp(-np.pi * tp**2)
h = np.exp(-np.pi * (tp - 0.08)**2)
alpha, beta = 1.7, -0.4
lhs = np.fft.fft(alpha*g + beta*h)
rhs = alpha*np.fft.fft(g) + beta*np.fft.fft(h)
err_lin = np.max(np.abs(lhs-rhs)) / np.max(np.abs(lhs))
# Retardo por fase: desplaza el espectro por exp(-j 2π f t0).
fpos = np.fft.fftfreq(len(tp), d=dt)
t0 = 0.037
G = np.fft.fft(g)
Gdelay = np.fft.fft(np.roll(g, int(round(t0/dt))))
err_delay = np.max(np.abs(Gdelay - G*np.exp(-2j*np.pi*fpos*t0))) / np.max(np.abs(G))
# Escalamiento temporal: g(a t) corresponde a (1/|a|)G(f/a); ambos lados se evalúan analíticamente.
fplot = np.linspace(-12, 12, 1200)
Gf = np.exp(-np.pi*fplot**2)
scaled_t = np.exp(-np.pi*(a*tp)**2)
Gscaled = np.fft.fftshift(np.fft.fft(scaled_t))*dt
fshift = np.fft.fftshift(fpos)
Gscaled_interp = np.interp(fplot, fshift, np.real(Gscaled))
Gscale_ideal = (1/abs(a))*np.exp(-np.pi*(fplot/a)**2)
err_scale = np.max(np.abs(Gscaled_interp-Gscale_ideal)) / np.max(Gscale_ideal)

fig, axs = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
axs[0].bar([0, 1], [err_lin, err_lin], color=[BLUE, GREEN])
axs[0].set_xticks([0, 1], ["αx+βy", "αX+βY"])
axs[0].set_title(f"Linealidad\nerror relativo = {err_lin:.1e}")
axs[1].plot(fpos, np.unwrap(np.angle(Gdelay/G)), color=ORANGE, lw=1)
axs[1].set_xlim(-80, 80)
axs[1].set_title("Retardo: fase lineal\n−2π f t₀")
axs[1].set_xlabel("Frecuencia (Hz)")
axs[2].plot(fplot, Gscale_ideal, "--", color="black", label="1/|a| · G(f/a)")
axs[2].plot(fplot, Gscaled_interp, color=PURPLE, label="FFT de g(at)")
axs[2].set_title(f"Escalamiento temporal\nerror relativo ≈ {err_scale:.2e}")
axs[2].set_xlabel("Frecuencia (Hz)")
axs[2].legend(frameon=False, fontsize=9)
for ax in axs:
    ax.grid(alpha=.22)
    if ax is axs[0]: ax.set_ylabel("Error relativo")
fig.suptitle("Propiedades de la transformada de Fourier", fontsize=17, fontweight="bold")
fig.savefig(FIG / "03_propiedades.png", dpi=180, bbox_inches="tight")

# Comparación de duración del pulso: pulso más estrecho, espectro más ancho.
tlocal = np.arange(-0.5, 0.5, 1/Fs)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
for width, color in [(0.04, ORANGE), (0.12, BLUE)]:
    x = (np.abs(tlocal) <= width/2).astype(float)
    X = np.fft.fftshift(np.fft.fft(x))/len(x)
    f = np.fft.fftshift(np.fft.fftfreq(len(x), d=1/Fs))
    axs[0].plot(tlocal, x, color=color, label=f"T = {width:.2f} s")
    axs[1].plot(f, np.abs(X), color=color, label=f"T = {width:.2f} s")
axs[0].set_xlim(-.2, .2); axs[0].set_xlabel("Tiempo (s)"); axs[0].set_ylabel("Amplitud")
axs[1].set_xlim(-100, 100); axs[1].set_xlabel("Frecuencia (Hz)"); axs[1].set_ylabel("|X[k]|")
axs[0].set_title("Duración del pulso"); axs[1].set_title("Ancho del espectro")
for ax in axs: ax.grid(alpha=.22); ax.legend(frameon=False)
fig.suptitle("Relación tiempo–frecuencia", fontsize=17, fontweight="bold")
fig.savefig(FIG / "04_tiempo_frecuencia.png", dpi=180, bbox_inches="tight")

print(f"Gráficos guardados en: {FIG}")
print(f"Error relativo linealidad: {err_lin:.3e}")
print(f"Error relativo retardo (muestreo entero): {err_delay:.3e}")
print(f"Error relativo escalamiento (malla finita): {err_scale:.3e}")
plt.show()
