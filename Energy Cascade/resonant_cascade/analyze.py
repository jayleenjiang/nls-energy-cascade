import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.optimize import curve_fit
import seaborn as sns

def analyze_cascade_results(evolution_file='cascade_evolution.csv', 
                          distribution_file='energy_distribution.csv'):
    """
    Analyze the cascade simulation results and create visualizations
    """
    
    # Load data
    evolution_df = pd.read_csv(evolution_file)
    distribution_df = pd.read_csv(distribution_file)
    
    # Set up plot style
    plt.style.use('seaborn-v0_8-darkgrid')
    fig = plt.figure(figsize=(16, 10))
    
    # 1. Energy evolution by generation
    ax1 = plt.subplot(2, 3, 1)
    generations = evolution_df['generation'].unique()
    times = evolution_df['time'].unique()
    
    for gen in generations:
        gen_data = evolution_df[evolution_df['generation'] == gen]
        ax1.plot(gen_data['time'], gen_data['total_energy'], 
                label=f'Gen {gen}', linewidth=2)
    
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Total Energy')
    ax1.set_title('Energy Evolution by Generation')
    ax1.legend(loc='best', ncol=2)
    ax1.set_yscale('log')
    
    # 2. Energy distribution at final time
    ax2 = plt.subplot(2, 3, 2)
    final_time = evolution_df['time'].max()
    final_data = evolution_df[evolution_df['time'] == final_time]
    
    gen_energies = []
    gen_labels = []
    for gen in generations:
        gen_energy = final_data[final_data['generation'] == gen]['total_energy'].values
        if len(gen_energy) > 0:
            gen_energies.append(gen_energy[0])
            gen_labels.append(f'Gen {gen}')
    
    ax2.bar(gen_labels, gen_energies, color='steelblue', edgecolor='black')
    ax2.set_ylabel('Total Energy')
    ax2.set_title(f'Energy Distribution at t={final_time:.1f}')
    ax2.set_yscale('log')
    ax2.tick_params(axis='x', rotation=45)
    
    # 3. Frequency-Energy relationship
    ax3 = plt.subplot(2, 3, 3)
    
    # Group by generation and plot
    colors = plt.cm.viridis(np.linspace(0, 1, len(generations)))
    for i, gen in enumerate(generations):
        gen_dist = distribution_df[distribution_df['generation'] == gen]
        if len(gen_dist) > 0:
            ax3.scatter(gen_dist['frequency_magnitude'], 
                       gen_dist['energy'],
                       label=f'Gen {gen}', 
                       alpha=0.7, 
                       s=50,
                       color=colors[i])
    
    ax3.set_xlabel('Frequency Magnitude |k|')
    ax3.set_ylabel('Energy')
    ax3.set_title('Energy vs Frequency')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.legend(loc='best', ncol=2)
    
    # 4. Energy cascade rate
    ax4 = plt.subplot(2, 3, 4)
    
    # Calculate energy flux between generations
    time_points = np.array(times[::10])  # Sample every 10th point
    flux_data = []
    
    for t in time_points:
        t_data = evolution_df[evolution_df['time'] == t]
        fluxes = []
        for gen in range(1, len(generations)):
            e_curr = t_data[t_data['generation'] == gen]['total_energy'].values
            e_next = t_data[t_data['generation'] == gen+1]['total_energy'].values
            if len(e_curr) > 0 and len(e_next) > 0:
                fluxes.append(e_next[0] - e_curr[0])
        flux_data.append(np.mean(fluxes) if fluxes else 0)
    
    ax4.plot(time_points, flux_data, 'b-', linewidth=2)
    ax4.set_xlabel('Time')
    ax4.set_ylabel('Mean Energy Flux')
    ax4.set_title('Energy Transfer Rate')
    ax4.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    ax4.grid(True, alpha=0.3)
    
    # 5. Power spectrum analysis
    ax5 = plt.subplot(2, 3, 5)
    
    # Fit power law to energy distribution
    freq_mags = distribution_df['frequency_magnitude'].values
    energies = distribution_df['energy'].values
    
    # Remove zeros and sort
    mask = (freq_mags > 0) & (energies > 1e-10)
    freq_mags = freq_mags[mask]
    energies = energies[mask]
    
    if len(freq_mags) > 10:
        # Bin the data for cleaner power law fit
        bins = np.logspace(np.log10(freq_mags.min()), 
                          np.log10(freq_mags.max()), 20)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        bin_energies = np.zeros(len(bin_centers))
        
        for i in range(len(bin_centers)):
            mask = (freq_mags >= bins[i]) & (freq_mags < bins[i+1])
            if np.any(mask):
                bin_energies[i] = np.mean(energies[mask])
        
        # Fit power law to non-zero bins
        valid = bin_energies > 0
        if np.sum(valid) > 3:
            def power_law(x, a, b):
                return a * x**b
            
            try:
                popt, _ = curve_fit(power_law, 
                                   bin_centers[valid], 
                                   bin_energies[valid],
                                   p0=[1, -2])
                
                ax5.scatter(bin_centers[valid], bin_energies[valid], 
                           s=50, alpha=0.7, label='Data')
                
                fit_x = np.logspace(np.log10(bin_centers[valid].min()),
                                   np.log10(bin_centers[valid].max()), 100)
                ax5.plot(fit_x, power_law(fit_x, *popt), 
                        'r-', linewidth=2,
                        label=f'Fit: E ~ k^{popt[1]:.2f}')
                
                # Add KZ prediction lines
                ax5.plot(fit_x, popt[0] * fit_x**(-7/3), 
                        'g--', linewidth=1.5, alpha=0.7,
                        label='KZ direct: k^{-7/3}')
                ax5.plot(fit_x, popt[0] * fit_x**(-5/3), 
                        'b--', linewidth=1.5, alpha=0.7,
                        label='KZ inverse: k^{-5/3}')
                
                ax5.set_xscale('log')
                ax5.set_yscale('log')
                ax5.set_xlabel('Frequency |k|')
                ax5.set_ylabel('Energy E(k)')
                ax5.set_title('Power Spectrum')
                ax5.legend()
                ax5.grid(True, alpha=0.3)
            except:
                print("Could not fit power law")
    
    # 6. Cascade efficiency over time
    ax6 = plt.subplot(2, 3, 6)
    
    # Calculate efficiency metric: ratio of high-freq to low-freq energy
    efficiency = []
    eff_times = []
    
    for t in times[::10]:
        t_data = evolution_df[evolution_df['time'] == t]
        
        # Low frequency: first 2 generations
        low_freq_e = 0
        for gen in range(1, min(3, len(generations)+1)):
            gen_e = t_data[t_data['generation'] == gen]['total_energy'].values
            if len(gen_e) > 0:
                low_freq_e += gen_e[0]
        
        # High frequency: last 2 generations  
        high_freq_e = 0
        for gen in range(max(1, len(generations)-1), len(generations)+1):
            gen_e = t_data[t_data['generation'] == gen]['total_energy'].values
            if len(gen_e) > 0:
                high_freq_e += gen_e[0]
        
        if low_freq_e > 1e-10:
            efficiency.append(high_freq_e / low_freq_e)
            eff_times.append(t)
    
    ax6.plot(eff_times, efficiency, 'k-', linewidth=2)
    ax6.set_xlabel('Time')
    ax6.set_ylabel('High-Freq / Low-Freq Energy')
    ax6.set_title('Cascade Efficiency')
    ax6.set_yscale('log')
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('cascade_analysis.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    # Print summary statistics
    print("\n=== Cascade Analysis Summary ===")
    print(f"Number of generations: {len(generations)}")
    print(f"Total simulation time: {final_time:.1f}")
    
    if len(efficiency) > 0:
        print(f"Final cascade efficiency: {efficiency[-1]:.3e}")
        print(f"Maximum efficiency reached: {max(efficiency):.3e}")
    
    # Calculate theoretical vs actual norm growth
    if len(distribution_df) > 0:
        early_gen = distribution_df[distribution_df['generation'] <= 3]
        late_gen = distribution_df[distribution_df['generation'] >= len(generations)-2]
        
        if len(early_gen) > 0 and len(late_gen) > 0:
            early_norm = np.sum(early_gen['frequency_magnitude']**4 * early_gen['energy'])
            late_norm = np.sum(late_gen['frequency_magnitude']**4 * late_gen['energy'])
            
            if early_norm > 0:
                actual_factor = late_norm / early_norm
                theoretical_factor = 2**(len(generations) - 5)
                
                print(f"\nNorm explosion analysis (s=2):")
                print(f"Theoretical factor: {theoretical_factor:.2f}")
                print(f"Actual factor: {actual_factor:.3e}")
                print(f"Ratio (actual/theoretical): {actual_factor/theoretical_factor:.3f}")

def create_cascade_animation(evolution_file='cascade_evolution.csv'):
    """
    Create an animation showing energy cascade through generations
    """
    evolution_df = pd.read_csv(evolution_file)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    generations = sorted(evolution_df['generation'].unique())
    times = sorted(evolution_df['time'].unique())
    
    # Initialize plots
    bars = ax1.bar(range(len(generations)), np.zeros(len(generations)))
    ax1.set_xlabel('Generation')
    ax1.set_ylabel('Total Energy (log scale)')
    ax1.set_title('Energy Cascade Animation')
    ax1.set_yscale('log')
    ax1.set_ylim(1e-6, 10)
    ax1.set_xticks(range(len(generations)))
    ax1.set_xticklabels([f'G{g}' for g in generations])
    
    lines = []
    for gen in generations:
        line, = ax2.plot([], [], label=f'Gen {gen}', linewidth=2)
        lines.append(line)
    
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Energy')
    ax2.set_title('Energy Evolution')
    ax2.set_xlim(0, times[-1])
    ax2.set_ylim(1e-6, 10)
    ax2.set_yscale('log')
    ax2.legend(loc='best', ncol=2)
    
    time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)
    
    def animate(frame):
        t = times[frame]
        t_data = evolution_df[evolution_df['time'] == t]
        
        # Update bar chart
        energies = []
        for gen in generations:
            gen_e = t_data[t_data['generation'] == gen]['total_energy'].values
            energies.append(gen_e[0] if len(gen_e) > 0 else 1e-10)
        
        for bar, e in zip(bars, energies):
            bar.set_height(max(e, 1e-10))
        
        # Update line plots
        for i, gen in enumerate(generations):
            gen_history = evolution_df[(evolution_df['generation'] == gen) & 
                                      (evolution_df['time'] <= t)]
            lines[i].set_data(gen_history['time'], gen_history['total_energy'])
        
        time_text.set_text(f'Time: {t:.2f}')
        
        return list(bars) + lines + [time_text]
    
    anim = FuncAnimation(fig, animate, frames=range(0, len(times), 5),
                        interval=50, blit=True)
    
    plt.tight_layout()
    anim.save('cascade_animation.gif', writer='pillow', fps=20)
    plt.show()
    
    print("Animation saved as 'cascade_animation.gif'")

if __name__ == "__main__":
    # Analyze results
    analyze_cascade_results()
    
    # Create animation
    create_cascade_animation()