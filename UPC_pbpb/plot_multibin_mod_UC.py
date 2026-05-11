import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from pathlib import Path

# ==============================================================
#  CONFIGURAZIONE — modifica qui i nomi dei file e le etichette
# ==============================================================

BINS = [
    {"label": r"$0.1 < z_1 < 0.2$", "file": "kT_max_50_kT_min_3_Vs_13000_thetac_0.03_z1_0.1_0.2_UC.txt"},
    {"label": r"$0.2 < z_1 < 0.3$", "file": "kT_max_50_kT_min_3_Vs_13000_thetac_0.03_z1_0.2_0.3_UC.txt"},
    {"label": r"$0.3 < z_1 < 0.4$", "file": "kT_max_50_kT_min_3_Vs_13000_thetac_0.03_z1_0.3_0.4_UC.txt"},
    {"label": r"$0.4 < z_1 < 0.5$", "file": "kT_max_50_kT_min_3_Vs_13000_thetac_0.03_z1_0.4_0.5_UC.txt"},
    {"label": r"$0.5 < z_1 < 0.7$", "file": "kT_max_50_kT_min_3_Vs_13000_thetac_0.03_z1_0.5_0.7_UC.txt"},
    {"label": r"$0.7 < z_1 < 0.9$", "file": "kT_max_50_kT_min_3_Vs_13000_thetac_0.03_z1_0.7_0.9_UC.txt"},
]

# Indici delle colonne (0-based) nel file .txt
COL_Z2 = 1   # colonna con z2 (asse x)
COL_Y  = 2   # colonna con il valore da plottare (curva U)
COL_Y2 = 6   # colonna con il valore da plottare (curva L)

# Informazioni fisiche stampate nel pannello in alto a sinistra (i==0)
VS_TEXT  = r"$\sqrt{s} = 92\ \mathrm{GeV}$"
Q0_TEXT  = r"$Q_0^2 = 3\ \mathrm{GeV}^2$"

# Etichetta dell'asse y (mostrata come fig.text verticale)
YLABEL = r"$\langle d\sigma_{unp}\vert \cos\phi_{12} \rangle$"

# Layout della griglia
N_ROWS, N_COLS = 2, 3

# Colori e stili delle curve
COLOR_U = "green"
COLOR_L = "orange"
LINE_STYLE = "solid"
LINE_STYLE_L = "dashed"
LINE_WIDTH = 2

# Limiti degli assi (None = automatico)
XLIM = (0.05, 0.85)
YLIM = (0.01, 0.04)   # None per automatico

# Posizioni verticali per i testi nel pannello i==0
Y_SQRTS  = 0.8
Y_Q20    = 0.7
# Posizione verticale per il testo del range z1 in ogni pannello
Y_ZRANGE =  0.95
# Posizioni delle linee U/L nel pannello i==2
Y_UC       = -0.06
Y_LC       = -0.045
Y_U_LABEL = -0.0625
Y_L_LABEL = -0.0475

# Aggiunge filigrana "Preliminary" su ogni pannello
SHOW_PRELIMINARY = True

# ==============================================================
#  FINE CONFIGURAZIONE
# ==============================================================

plt.rcParams["font.family"] = "Arial"
plt.rcParams["text.usetex"] = True
plt.rcParams["hatch.linewidth"] = 1.2

BASE_DIR = Path(__file__).resolve().parent


def load_data(filepath, col_x, col_y, col_y2):
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
    return data[:, col_x], data[:, col_y], data[:, col_y2]


def make_figure(bins, n_rows, n_cols):
    fig, axs = plt.subplots(nrows=n_rows, ncols=n_cols,
                            figsize=(4 * n_cols, 4 * n_rows),
                            sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0.0, hspace=0.0)

    # Etichetta y globale
    fig.text(0.05, 0.5, YLABEL, va="center", ha="center",
             rotation="vertical", fontsize=30)

    axes = axs.flatten()
    legend_handles = []  # per la legenda in alto (U / L)

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

        # Testi nel pannello i==0
        if i == 0:
            ax.text(0.1, Y_SQRTS, VS_TEXT, size=15, transform=ax.transAxes,
                    ha="left", va="top")
            ax.text(0.1, Y_Q20,   Q0_TEXT, size=15, transform=ax.transAxes,
                    ha="left", va="top")

        # Testo range z1
        ax.text(0.1, Y_ZRANGE, bin_cfg["label"], size=15,
                transform=ax.transAxes, ha="left", va="top")

        # --- Dati ---
        try:
            z2, y_uc, y_lc = load_data(bin_cfg["file"], COL_Z2, COL_Y, COL_Y2)
            lu, = ax.plot(z2, y_uc, lw=LINE_WIDTH, ls=LINE_STYLE,   color=COLOR_U)
            #ll, = ax.plot(z2, y_lc, lw=LINE_WIDTH, ls=LINE_STYLE_L, color=COLOR_L)
            if i == 0:
                legend_handles = [lu]  # , ll]
        except FileNotFoundError:
            print(f"ATTENZIONE: File non trovato: {bin_cfg['file']}")
            z2 = np.linspace(0.1, 0.9, 80)
            rng = np.random.default_rng(i)
            lu, = ax.plot(z2, 0.015 + 0.005 * rng.standard_normal(len(z2)),
                          lw=LINE_WIDTH, ls=LINE_STYLE, color=COLOR_U,
                          label="U [placeholder]")
            #ll, = ax.plot(z2, 0.010 + 0.005 * rng.standard_normal(len(z2)),
            #              lw=LINE_WIDTH, ls=LINE_STYLE_L, color=COLOR_L,
            #              label="L [placeholder]")
            if i == 0:
                legend_handles = [lu]  # , ll]
            ax.text(0.5, 0.5, "FILE NON TROVATO",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="red", alpha=0.5)

        # Legenda manuale U / L nel pannello i==2 (in alto a destra)
        if i == 2:
            x_seg = [0.1, 0.25]
            ax.plot(x_seg, [Y_UC] * 2, lw=2, ls=LINE_STYLE,   color="black")
            #ax.plot(x_seg, [Y_LC] * 2, lw=2, ls=LINE_STYLE_L, color="black")
            ax.text(0.3, Y_U_LABEL, r"$UC$", size=15)
            #ax.text(0.3, Y_L_LABEL, r"$L$", size=15)

        # Filigrana "Preliminary"
        if SHOW_PRELIMINARY:
            ax.text(0.5, 0.5, "Preliminary",
                    transform=ax.transAxes, fontsize=42,
                    color="gray", alpha=0.0,
                    ha="center", va="center", rotation=30)

    # Legenda globale in alto al centro
    if legend_handles:
        fig.legend(legend_handles, [r"$UC$"],  # , r"$L$"],
                   loc="outside upper center",
                   bbox_to_anchor=(0.4, 0.99),
                   ncol=2, frameon=False,
                   prop={"size": 26})

    return fig


if __name__ == "__main__":
    fig = make_figure(BINS, N_ROWS, N_COLS)
    plt.savefig("ratio_UC.pdf", dpi=150, bbox_inches="tight")
    plt.savefig("ratio_UC.png", dpi=150, bbox_inches="tight")
    #plt.show()
    print("Salvato: ratio_UC.pdf / ratio_UC.png")
