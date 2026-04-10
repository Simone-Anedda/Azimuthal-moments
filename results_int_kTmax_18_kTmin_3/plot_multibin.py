import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

# ==============================================================
#  CONFIGURAZIONE — modifica qui i nomi dei file e le etichette
# ==============================================================

# Ogni elemento corrisponde a un subplot (ordine: riga per riga)
# "label"  : stringa mostrata nell'angolo in alto a sinistra
# "file"   : percorso al file .txt con i dati
#            Il file deve avere (almeno) due colonne: z2  valore
#            Puoi aggiungere altre colonne e selezionarle con col_z2 / col_y

BINS = [
    {"label": r"$0.1 < z_1 < 0.2$", "file": "kT_max_18_kT_min_3_Vs_92_thetac_0.03_z1_0.1_0.2_U_L.txt"},
    {"label": r"$0.2 < z_1 < 0.3$", "file": "kT_max_18_kT_min_3_Vs_92_thetac_0.03_z1_0.2_0.3_U_L.txt"},
    {"label": r"$0.3 < z_1 < 0.4$", "file": "kT_max_18_kT_min_3_Vs_92_thetac_0.03_z1_0.3_0.4_U_L.txt"},
    {"label": r"$0.4 < z_1 < 0.5$", "file": "kT_max_18_kT_min_3_Vs_92_thetac_0.03_z1_0.4_0.5_U_L.txt"},
    {"label": r"$0.5 < z_1 < 0.6$", "file": "kT_max_18_kT_min_3_Vs_92_thetac_0.03_z1_0.5_0.7_U_L.txt"},
    {"label": r"$0.6 < z_1 < 0.7$", "file": "kT_max_18_kT_min_3_Vs_92_thetac_0.03_z1_0.7_0.9_U_L.txt"},
]

# Indici delle colonne (0-based) nel file .txt
COL_Z2 = 1   # colonna con z2 (asse x)
COL_Y  = 2   # colonna con il valore da plottare (asse y)
COL_Y2  = 6
# Informazioni fisiche stampate nel pannello centrale in alto
INFO_TEXT = r"$\sqrt{s} = 92\ \mathrm{GeV}$" + "\n" + r"$Q_0^2 = 3\ \mathrm{GeV}^2$"

# Etichetta dell'asse y
YLABEL = r"$\langle d\sigma | 1, 1 \rangle$"

# Layout della griglia (righe x colonne)
N_ROWS, N_COLS = 2, 3

# Stile della curva
LINE_COLOR  = "steelblue"
LINE_STYLE  = "-"
LINE_WIDTH  = 1.8

# Limiti degli assi (None = automatico)
XLIM = (0.1, 1.0)
YLIM = None          # es. (-0.03, 0.05)

# ==============================================================
#  FINE CONFIGURAZIONE
# ==============================================================


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
        data = np.genfromtxt(
            filepath,
            comments="#",
            delimiter=",",
            skip_header=1,
            invalid_raise=False,
        )
    else:
        data = np.loadtxt(filepath, comments="#")

    data = np.atleast_2d(data)
    data = data[~np.all(np.isnan(data), axis=1)]
    return data[:, col_x], data[:, col_y], data[:, col_y2]


def make_figure(bins, n_rows, n_cols):
    fig = plt.figure(figsize=(4.5 * n_cols, 3.5 * n_rows))
    gs  = gridspec.GridSpec(n_rows, n_cols, figure=fig,
                            hspace=0.08, wspace=0.08)

    axes = []
    for idx, bin_cfg in enumerate(bins):
        row, col = divmod(idx, n_cols)
        ax = fig.add_subplot(gs[row, col])
        axes.append(ax)

        # --- carica dati ---
        try:
            z2, y, y2 = load_data(bin_cfg["file"], COL_Z2, COL_Y, COL_Y2)

            ax.plot(z2, y, color=LINE_COLOR,
                    ls=LINE_STYLE, lw=LINE_WIDTH, label = "U")
            ax.plot(z2, y2, color="orange",
                    ls=LINE_STYLE, lw=LINE_WIDTH,
                    label="L")
        except FileNotFoundError:
            # Se il file non esiste genera dati fittizi come placeholder
            print(f"ATTENZIONE: File non trovato: {bin_cfg['file']}")
            z2 = np.linspace(0.1, 0.9, 80)
            y  = 0.015 * np.ones_like(z2) + 0.005 * np.random.randn(len(z2))
            y2 = 0.010 * np.ones_like(z2) + 0.005 * np.random.randn(len(z2))
            ax.plot(z2, y, color=LINE_COLOR,
                    label="[placeholder — file non trovato]")   
            ax.plot(z2, y2, color="orange",
                    ls=LINE_STYLE, lw=LINE_WIDTH,
                    label="[placeholder — file non trovato]")
            ax.text(0.5, 0.5, "FILE NON TROVATO",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="red", alpha=0.5)

        # --- linea a zero ---
        ax.axhline(0, color="gray", lw=0.7, ls="--")

        # --- etichetta bin ---
        ax.text(0.05, 0.95, bin_cfg["label"],
                transform=ax.transAxes, ha="left", va="top",
                fontsize=10)

        # --- limiti assi ---
        if XLIM:
            ax.set_xlim(XLIM)
        if YLIM:
            ax.set_ylim(YLIM)

        # --- etichette assi (solo bordo esterno) ---
        if row == n_rows - 1:
            ax.set_xlabel(r"$z_2$", fontsize=12)
        else:
            ax.set_xticklabels([])

        if col == 0:
            ax.set_ylabel(YLABEL, fontsize=12)
        else:
            ax.set_yticklabels([])

    # --- testo informativo nel pannello vuoto in alto al centro ---
    # (pannello indice 1 nella griglia 2x3, posizione centrale)
    info_ax_idx = 1   # cambia se usi layout diverso
    if info_ax_idx < len(axes):
        axes[info_ax_idx].text(
            0.05, 0.08, INFO_TEXT,
            transform=axes[info_ax_idx].transAxes,
            ha="left", va="bottom", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.7)
        )

    fig.tight_layout()
    return fig


if __name__ == "__main__":
    fig = make_figure(BINS, N_ROWS, N_COLS)
    plt.savefig("plot_multibin.pdf", dpi=150, bbox_inches="tight")
    plt.savefig("plot_multibin.png", dpi=150, bbox_inches="tight")
    plt.show()
    print("Salvato: plot_multibin.pdf / plot_multibin.png")
