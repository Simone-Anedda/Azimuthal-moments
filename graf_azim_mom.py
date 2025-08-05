import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def read_data_file(filename):
    """
    Legge il file di dati e organizza i dati in formato utilizzabile
    Gestisce automaticamente i valori NaN e non numerici
    """
    # Leggi il file saltando le righe vuote e gestendo gli spazi
    data = []
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    # Trova la riga con gli header (quella che contiene "err[1]", "err[2]", ecc.)
    header_line = None
    data_start = None
    
    for i, line in enumerate(lines):
        # Cerca una riga che contenga almeno due pattern "err[" 
        # (questo indica che è la riga degli header)
        if line.count('err[') >= 2:
            header_line = i
            data_start = i + 1
            break
    
    if header_line is None:
        raise ValueError("Non riesco a trovare la riga degli header nel file. Assicurati che ci siano almeno due colonne con pattern 'err[n]'")
    
    # Estrai gli header
    headers = lines[header_line].split()
    
    # Leggi i dati numerici
    data_lines = []
    for line_num, line in enumerate(lines[data_start:], start=data_start+1):
        line = line.strip()
        if line and not line.startswith('#'):  # Salta righe vuote e commenti
            # Dividi la riga in numeri
            values = line.split()
            if len(values) >= len(headers):  # Assicurati che ci siano abbastanza valori
                try:
                    # Converti in float, gestendo NaN e valori non numerici
                    numeric_values = []
                    for val in values[:len(headers)]:
                        try:
                            # Gestisci diversi formati di NaN
                            if val.lower() in ['nan', 'na', 'null', '']:
                                numeric_values.append(np.nan)
                            else:
                                numeric_values.append(float(val))
                        except ValueError:
                            # Se non riesce a convertire, metti NaN
                            numeric_values.append(np.nan)
                    
                    data_lines.append(numeric_values)
                except Exception as e:
                    print(f"Attenzione: Errore nella riga {line_num}, saltata: {e}")
                    continue
    
    # Converti in array numpy
    data_array = np.array(data_lines)
    
    # Statistiche sui NaN
    total_values = data_array.size
    nan_count = np.isnan(data_array).sum()
    if nan_count > 0:
        print(f"Trovati {nan_count} valori NaN su {total_values} totali ({nan_count/total_values*100:.1f}%)")
    
    return headers, data_array

def create_scatter_plots(filename, save_plots=True, show_plots=True):
    """
    Crea scatter plot separati per ogni coppia y-errore vs x
    Esclude automaticamente i punti con valori NaN
    """
    # Leggi i dati
    headers, data = read_data_file(filename)
    
    # La prima colonna è x
    x_label = headers[0]
    x_data = data[:, 0]
    
    # Identifica le coppie y-errore
    y_pairs = []
    i = 1
    while i < len(headers) - 1:
        if 'err[' in headers[i+1]:  # Il prossimo è un errore
            y_label = headers[i]
            err_label = headers[i+1]
            y_pairs.append((i, i+1, y_label, err_label))
            i += 2
        else:
            i += 1
    
    print(f"Trovate {len(y_pairs)} coppie di dati y-errore")
    print(f"Variabile x: {x_label}")
    
    # Crea i grafici
    for idx, (y_col, err_col, y_label, err_label) in enumerate(y_pairs):
        plt.figure(figsize=(10, 6))
        
        y_data = data[:, y_col]
        err_data = data[:, err_col]
        
        # Filtra i valori NaN
        # Crea una maschera per identificare i valori validi (non NaN)
        valid_mask = ~(np.isnan(x_data) | np.isnan(y_data) | np.isnan(err_data))
        
        # Applica la maschera per ottenere solo i valori validi
        x_valid = x_data[valid_mask]
        y_valid = y_data[valid_mask]
        err_valid = err_data[valid_mask]
        
        # Conta i punti esclusi
        excluded_points = len(x_data) - len(x_valid)
        if excluded_points > 0:
            print(f"  {y_label}: Esclusi {excluded_points} punti con valori NaN")
        
        # Verifica se ci sono ancora dati da plottare
        if len(x_valid) == 0:
            print(f"  ATTENZIONE: {y_label} - Tutti i punti contengono NaN, salto questo grafico")
            plt.close()
            continue
        
        # Crea scatter plot con error bars
        plt.errorbar(x_valid, y_valid, yerr=err_valid, 
                    fmt='o', capsize=3, capthick=1, 
                    markersize=4, alpha=0.7,
                    label=f'{y_label} ({len(x_valid)} punti)')
        
        plt.xlabel(x_label, fontsize=12)
        plt.ylabel(y_label, fontsize=12)
        plt.title(f'{y_label} vs {x_label}', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Migliora il layout
        plt.tight_layout()
        
        if save_plots:
            # Pulisci il nome del file per il salvataggio
            clean_y_label = y_label.replace('<', '').replace('>', '').replace('|', '_').replace(',', '_')
            import os
            plot_dir = os.path.join('plots', os.path.splitext(filename)[0])
            os.makedirs(plot_dir, exist_ok=True)
            plt.savefig(os.path.join(plot_dir, f'plot_{idx+1}_{clean_y_label}_vs_{x_label}.png'), 
                       dpi=300, bbox_inches='tight')
            print(f"Salvato: plot_{idx+1}_{clean_y_label}_vs_{x_label}.png")
        
        if show_plots:
            plt.show()
        else:
            plt.close()

def analyze_data_structure(filename):
    """
    Analizza la struttura del file per debug
    Include statistiche sui valori NaN
    """
    headers, data = read_data_file(filename)
    
    print("=== ANALISI STRUTTURA DATI ===")
    print(f"Numero di colonne: {len(headers)}")
    print(f"Numero di righe di dati: {data.shape[0]}")
    print("\nHeaders trovati:")
    for i, header in enumerate(headers):
        nan_count = np.isnan(data[:, i]).sum()
        nan_percent = (nan_count / len(data)) * 100
        print(f"  Colonna {i}: {header} ({nan_count} NaN, {nan_percent:.1f}%)")
    
    print(f"\nPrime 5 righe di dati:")
    print(data[:5])
    
    # Statistiche generali sui NaN
    total_nan = np.isnan(data).sum()
    total_values = data.size
    print(f"\nStatistiche NaN generali:")
    print(f"  Totale valori NaN: {total_nan}")
    print(f"  Percentuale NaN: {total_nan/total_values*100:.2f}%")
    
    # Identifica righe completamente valide
    rows_with_nan = np.any(np.isnan(data), axis=1).sum()
    valid_rows = len(data) - rows_with_nan
    print(f"  Righe con almeno un NaN: {rows_with_nan}")
    print(f"  Righe completamente valide: {valid_rows}")
    
    return headers, data

# Esempio di utilizzo
if __name__ == "__main__":
    filenames = ["Fixed_xB.txt", "Fixed_y.txt", "Fixed_Q2.txt"]  # Sostituisci con i nomi dei tuoi file

    for filename in filenames:
        print(f"\n=== Elaborazione file: {filename} ===")
        try:
            # Analizza prima la struttura (opzionale, per debug)
            print("Analizzando la struttura del file...")
            analyze_data_structure(filename)
            
            print("\n" + "="*50)
            print("Creando i grafici...")
            
            # Crea i grafici
            create_scatter_plots(filename, save_plots=True, show_plots=False)
            
            print("Processo completato!")
            
        except FileNotFoundError:
            print(f"Errore: File '{filename}' non trovato!")
            print("Assicurati che il file sia nella stessa directory dello script.")
        except Exception as e:
            print(f"Errore durante l'elaborazione: {e}")
        
    # Esempio di come usare solo alcune funzioni
    # analyze_data_structure(filename)  # Solo per analizzare
    # create_scatter_plots(filename, save_plots=False, show_plots=True)  # Solo mostra