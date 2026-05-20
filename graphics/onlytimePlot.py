import json
import matplotlib.pyplot as plt

nome_file = 'dati.json'

try:
    with open(nome_file, 'r') as f:
        data = json.load(f)
except FileNotFoundError:
    print(f"Errore: Il file '{nome_file}' non è stato trovato.")
    data = {}

if data:
    # 1. Raggruppamento rigoroso
    groups = {}
    for v in data.values():
        p_count = v['particles_per_type']
        if p_count not in groups:
            groups[p_count] = []
        groups[p_count].append((v['threads'], v['time_seconds']))

    # 2. Creazione Grafici
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    for p_count in sorted(groups.keys()):
        # ORDINAMENTO FONDAMENTALE: i thread devono essere in ordine (1, 2, 4...)
        values = sorted(groups[p_count], key=lambda x: x[0])

        threads = [x[0] for x in values]
        times = [x[1] for x in values]

        # CERCHIAMO IL TEMPO A 1 THREAD (T1)
        # Se non esiste il test con 1 thread per questo gruppo, lo speedup non è calcolabile correttamente
        t1_list = [t for n, t in values if n == 1]

        if not t1_list:
            print(f"Attenzione: Manca il test a 1 thread per {p_count} particelle. Speedup saltato.")
            continue

        t1 = t1_list[0]

        # CALCOLO SPEEDUP: T1 / Tn
        # Se Tn è quasi uguale a T1, lo speedup sarà ~1 (scarso)
        speedup = [t1 / tn for tn in times]

        label = f"{p_count} particles"
        ax1.plot(threads, times, marker='o', label=label, linewidth=2)
        ax2.plot(threads, speedup, marker='s', label=label, linewidth=2)

    # 3. Estetica e Verifica Ideale
    ax1.set_title('Execution Time ')
    ax1.set_ylabel('Seconds')
    ax1.set_xlabel('Threads')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Linea Ideale: Se raddoppio i thread, il tempo dovrebbe dimezzarsi (Speedup = n_threads)
    all_threads = [n for n, t in [item for sublist in groups.values() for item in sublist]]
    max_th = max(all_threads)
    ax2.plot([1, max_th], [1, max_th], color='gray', linestyle='--', label='Ideal (Linear)')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_title('Speedup')
    ax2.set_ylabel('Speedup Ratio (T1/Tn)')
    ax2.set_xlabel('Threads')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    plt.savefig('SAopt_onlytimes.png', dpi=300)
    plt.show()