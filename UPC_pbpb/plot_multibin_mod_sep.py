import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from pathlib import Path

# ==============================================================
#  CONFIGURAZIONE — modifica qui i nomi dei file e le etichette
# ==============================================================

BINS = [
    {"label": r"$0.1 < z_1 < 0.2$", "file": "sep/kT_max_50_kT_min_3_Vs_13000_thetac_0.03_z1_0.1_0.2_separated.txt"},
    {"label": r"$0.2 < z_1 < 0.3$", "file": "sep/kT_max_50_kT_min_3_Vs_13000_thetac_0.03_z1_0.2_0.3_separated.txt"},
    {"label": r"$0.3 < z_1 < 0.4$", "file": "sep/kT_max_50_kT_min_3_Vs_13000_thetac_0.03_z1_0.3_0.4_separated.txt"},
    {"label": r"$0.4 < z_1 < 0.5$", "file": "sep/kT_max_50_kT_min_3_Vs_13000_thetac_0.03_z1_0.4_0.5_separated.txt"},
    {"label": r"$0.5 < z_1 < 0.7$", "file": "sep/kT_max_50_kT_min_3_Vs_13000_thetac_0.03_z1_0.5_0.7_separated.txt"},
    {"label": r"$0.7 < z_1 < 0.9$", "file": "sep/kT_max_50_kT_min_3_Vs_13000_thetac_0.03_z1_0.7_0.9_separated.txt"},
]

# Indici delle colonne (0-based) nel file .txt
COL_Z2 = 1
COL_PM = 2    # pi+ pi-
COL_MP = 6    # pi- pi+
COL_PP = 10   # pi+ pi+
COL_MM = 14   # pi- pi-

# Informazioni fisiche stampate nel pannello in alto a sinistra (i==0)
VS_TEXT = r"$\sqrt{s} = 92\ \mathrm{GeV}$"
Q0_TEXT = r"$Q_0^2 = 3\ \mathrm{GeV}^2$"

# Etichetta dell'asse y (mostrata come fig.text verticale)
YLABEL = r"$\langle d\sigma_{unp}\vert \cos\phi_{12} \rangle$"

# Colori e stili delle quattro curve
COLOR_PM = "green"
COLOR_MP = "red"
COLOR_PP = "blue"
COLOR_MM = "magenta"
LINE_STYLE = "solid"
LINE_STYLE_L = "dashed"
LINE_WIDTH = 2

# Labels per la legenda
LABELS = [r"$\pi^+ \pi^-$", r"$\pi^- \pi^+$", r"$\pi^+ \pi^+$", r"$\pi^- \pi^-$"]

# Layout della griglia
N_ROWS, N_COLS = 2, 3

# Limiti degli assi (None = automatico)
XLIM = (0.05, 0.85)
YLIM = (-0.08, 0.08)

# Posizioni testi nel pannello i==0 (coordinate axes, 0–1)
Y_SQRTS  = 0.30
Y_Q20    = 0.20

# Posizione testo range z1 in ogni pannello (coordinate axes, 0–1)
Y_ZRANGE = 0.95

# Aggiunge filigrana "Preliminary" su ogni pannello
SHOW_PRELIMINARY = False
PRELIMINARY_ALPHA = 0.4   # metti 0.0 per nasconderla

# ==============================================================
#  FINE CONFIGURAZIONE
# ==============================================================

plt.rcParams["font.family"] = "Arial"
plt.rcParams["text.usetex"] = True
plt.rcParams["hatch.linewidth"] = 1.2

BASE_DIR = Path(__file__).resolve().parent


def load_data(filepath, col_x, col_pm, col_mp, col_pp, col_mm):
    """Carica colonne da file CSV o whitespace-separated."""
    filepath = Path(filepath)
    if not filepath.is_absolute():
        filepath = BASE_DIR / filepath

    with filepath.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            first_data_line = stripped
            break
        else:
            raise ValueError(f"File vuoto o senza dati validi: {filepath}")

    if "," in first_data_line:
        data = np.genfromtxt(filepath, comments="#", delimiter=",",
                             skip_header=1, invalid_raise=False)
    else:
        data = np.loadtxt(filepath, comments="#")

    data = np.atleast_2d(data)
    data = data[~np.all(np.isnan(data), axis=1)]
    return (data[:, col_x], data[:, col_pm], data[:, col_mp],
            data[:, col_pp], data[:, col_mm])


def make_figure(bins, n_rows, n_cols):
    fig, axs = plt.subplots(nrows=n_rows, ncols=n_cols,
                            figsize=(4 * n_cols, 4 * n_rows),
                            sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0.0, hspace=0.0, left=0.1)

    # Etichetta y globale
    fig.text(0.04, 0.5, YLABEL, va="center", ha="center",
             rotation="vertical", fontsize=20)

    axes = axs.flatten()
    legend_handles = []

    for i, bin_cfg in enumerate(bins):
        ax = axes[i]

        # Limiti assi
        if XLIM:
            ax.set_xlim(XLIM)
        if YLIM:
            ax.set_ylim(YLIM)

        # Tick style
        ax.tick_params("both", direction="in", top=True, right=True,
                       labelsize=15, length=5, which="major")
        ax.tick_params("both", direction="in", top=True, right=True,
                       labelsize=15, length=3, which="minor")
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.set_axisbelow(False)

        # Asse x solo nella riga in basso
        if i >= (n_rows - 1) * n_cols:
            ax.set_xlabel(r"$z_2$", size=20)

        # Linea a zero
        ax.axhline(0, color="black", lw=1, zorder=0)

        # Testi fisici nel pannello i==0 (coordinate axes)
        if i == 0:
            ax.text(0.1, Y_SQRTS, VS_TEXT, size=15, transform=ax.transAxes,
                    ha="left", va="top")
            ax.text(0.1, Y_Q20,   Q0_TEXT, size=15, transform=ax.transAxes,
                    ha="left", va="top")

        # Testo range z1 (coordinate axes)
        ax.text(0.1, Y_ZRANGE, bin_cfg["label"], size=15,
                transform=ax.transAxes, ha="left", va="top")

        # --- Dati ---
        try:
            z2, y_PM, y_MP, y_PP, y_MM = load_data(
                bin_cfg["file"], COL_Z2, COL_PM, COL_MP, COL_PP, COL_MM)
            lPM, = ax.plot(z2, y_PM, lw=LINE_WIDTH, ls=LINE_STYLE,
                           color=COLOR_PM, label=LABELS[0])
            lMP, = ax.plot(z2, y_MP, lw=LINE_WIDTH, ls=LINE_STYLE_L,
                           color=COLOR_MP, label=LABELS[1])
            lPP, = ax.plot(z2, y_PP, lw=LINE_WIDTH, ls=LINE_STYLE,
                           color=COLOR_PP, label=LABELS[2])
            lMM, = ax.plot(z2, y_MM, lw=LINE_WIDTH, ls=LINE_STYLE_L,
                           color=COLOR_MM, label=LABELS[3])
            if i == 0:
                legend_handles = [lPM, lMP, lPP, lMM]

        except FileNotFoundError:
            print(f"ATTENZIONE: File non trovato: {bin_cfg['file']}")
            z2 = np.linspace(0.1, 0.9, 80)
            rng = np.random.default_rng(i)
            handles_placeholder = []
            for color, label in zip(
                    [COLOR_PM, COLOR_MP, COLOR_PP, COLOR_MM], LABELS):
                h, = ax.plot(z2, 0.01 * rng.standard_normal(len(z2)),
                             lw=LINE_WIDTH, ls=LINE_STYLE,
                             color=color, label=label)
                handles_placeholder.append(h)
            if i == 0:
                legend_handles = handles_placeholder
            ax.text(0.5, 0.5, "FILE NON TROVATO",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="red", alpha=0.5)

        # Filigrana "Preliminary"
        if SHOW_PRELIMINARY:
            ax.text(0.5, 0.5, "Preliminary",
                    transform=ax.transAxes, fontsize=36,
                    color="gray", alpha=PRELIMINARY_ALPHA,
                    ha="center", va="center", rotation=30)

    # Legenda globale in alto al centro con le label corrette
    if legend_handles:
        fig.legend(legend_handles, LABELS,
                   loc="outside upper center",
                   bbox_to_anchor=(0.55, 0.99),
                   ncol=4, frameon=False,
                   prop={"size": 18})

    return fig


if __name__ == "__main__":
    fig = make_figure(BINS, N_ROWS, N_COLS)
    plt.savefig("separated.pdf", dpi=150, bbox_inches="tight")
    plt.savefig("separated.png", dpi=150, bbox_inches="tight")
    plt.show()
    print("Salvato: separated.pdf / separated.png")