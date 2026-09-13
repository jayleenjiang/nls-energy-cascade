import pandas as pd
import matplotlib.pyplot as plt
import re
import os

files = ['energy_n200_no_summation.csv']
plot_filename = 'new.png'

plt.figure(figsize=(10, 7))
ax = plt.gca()

found_files = 0
for filename in files:
    if os.path.exists(filename):
        df = pd.read_csv(filename)

        n_value = re.search(r'n(\d+)', filename).group(1)
        label = f'n = {n_value}'
        
        ax.plot(df['mode'][:-1], df['energy'][:-1], 'o-', label=label, linewidth=2, markersize=6)
        found_files += 1
    else:
        print(f"Warning: Could not find file '{filename}'. Skipping.")

if found_files > 0:
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_title('Energy Distribution', fontsize=16)
    ax.set_xlabel('Mode Number (j)', fontsize=12)
    ax.set_ylabel('Energy (Σt Ij)', fontsize=12)
    ax.legend(loc='lower left', fontsize=12)
    ax.grid(True, which="both", linestyle='--', linewidth=0.5)

    plt.savefig(plot_filename)
    print(f"\nSuccess! Plot saved as '{plot_filename}'")
else:
    print("\nError: No data files were found. Could not generate plot.")