import pandas as pd
import matplotlib.pyplot as plt

def crea_grafici_integrali(file_path):
    """
    Crea grafici separati per ogni integrale vs Q2
    
    Args:
        file_path (str): Percorso del file di dati
    """
    # Carica i dati
    data = pd.read_csv(file_path, sep=r'\s+')
    
    # Rimuove la prima riga e le ultime 3 righe
    # Rimuove le righe che hanno valori NaN
    data = data.dropna()

    print(f"Dati caricati: {len(data)} righe")

    # Prendi il nome della prima colonna (es: 'Q2', 'xB', ecc.)
    x_col = data.columns[0]
    print(f"Colonna x usata per l'asse X: {x_col}")

    # Trova tutte le colonne degli integrali (escludi la prima colonna e colonne di errore)
    int_columns = [col for col in data.columns if col != x_col and not col.startswith('err[')]

    print(f"Trovate {len(int_columns)} colonne di integrali")

    # Crea un grafico per ogni integrale
    for i, int_col in enumerate(int_columns, start=1):
        plt.figure(figsize=(8, 6))
        err_col = f'err[{i}]'
        err_col = f'err[{i}]'

        # Rimuove i valori NaN
        if err_col in data.columns:
            mask = ~(data[x_col].isna() | data[int_col].isna() | data[err_col].isna())
            x_data = data[x_col][mask]
            y_data = data[int_col][mask]
            err_data = data[err_col][mask]
            plt.errorbar(x_data, y_data, yerr=err_data, fmt='o', capsize=3, 
                         capthick=1, elinewidth=1, markersize=4, alpha=0.7)
        else:
            mask = ~(data[x_col].isna() | data[int_col].isna())
            x_data = data[x_col][mask]
            y_data = data[int_col][mask]
            plt.scatter(x_data, y_data, alpha=0.7, s=50)

        plt.xlabel(x_col)
        plt.ylabel(int_col)
        plt.title(f'{int_col} vs {x_col}')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{int_col}_vs_{x_col}.png")
        plt.show()

# Uso del programma
if __name__ == "__main__":
    # Sostituisci con il percorso del tuo file
    file_path = "Prova2.txt"
    
    crea_grafici_integrali(file_path)
    # Per Prova.txt: variabile xB
    crea_grafici_integrali("Prova.txt")

    # Per Prova1.txt: variabile y
    crea_grafici_integrali("Prova1.txt")