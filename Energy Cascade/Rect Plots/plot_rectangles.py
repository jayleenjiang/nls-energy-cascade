import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import re
import numpy as np

def format_rational(x):
    if isinstance(x, str) and '/' in x:
        try:
            num, den = map(float, x.split('/'))
            return num / den if den != 0 else 0 
        except ValueError:
            return 0 
    try:
        return float(x)
    except (ValueError, TypeError):
        return 0


def create_labeled_separate_plots(filenames):
    for filename in filenames:
        try:
            df = pd.read_csv(filename)
        except pd.errors.EmptyDataError:
            print(f"Warning: CSV file '{filename}' is empty. Skipping.")
            continue
        except FileNotFoundError:
            print(f"Warning: CSV file '{filename}' not found. Skipping.")
            continue
            
        df['real'] = df['real'].apply(format_rational)
        df['imag'] = df['imag'].apply(format_rational)

        plt.figure(figsize=(12, 12))
        ax = plt.gca()

        # Group by family to draw each rectangle
        for family_id, group in df.groupby('family_id'):
            # Ensure there are enough points for parents and children
            if len(group) < 4: continue
            parents = group[group['point_type'] == 'parent']
            children = group[group['point_type'] == 'child']
            if len(parents) < 2 or len(children) < 2: continue

            # Plot the points
            ax.scatter(parents['real'], parents['imag'], c='royalblue', s=50, zorder=5, label='Parents' if family_id == 0 else "")
            ax.scatter(children['real'], children['imag'], c='firebrick', s=50, zorder=5, label='Children' if family_id == 0 else "")

            p1, p2 = parents.iloc[0], parents.iloc[1]
            c1, c2 = children.iloc[0], children.iloc[1]
            
            rect_x = [p1['real'], c2['real'], p2['real'], c1['real'], p1['real']]
            rect_y = [p1['imag'], c2['imag'], p2['imag'], c1['imag'], p1['imag']]
            
            ax.plot(rect_x, rect_y, color='black', alpha=0.6, linewidth=1)

            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            x_offset = (xlim[1] - xlim[0]) * 0.01 if (xlim[1] - xlim[0]) != 0 else 0.1
            y_offset = (ylim[1] - ylim[0]) * 0.015 if (ylim[1] - ylim[0]) != 0 else 0.1
            
            for _, point in group.iterrows():
                label = f"({point['real']:.1f}, {point['imag']:.1f}i)"
                ax.text(point['real'] + x_offset, point['imag'] + y_offset, label, fontsize=8, ha='center', color='#333')

        base_name = os.path.splitext(os.path.basename(filename))[0]
        match = re.search(r'N(\d+)_R(\d+)', filename)
        title_params = f" (N={match.group(1)}, R={match.group(2)})" if match else ""
        base_title = ' '.join(part for part in base_name.split('_') if part.startswith('f')).replace('f','F').replace('To',' To ') # Improved title extraction
        plt.title(f'Nuclear Family Rectangles {base_title}{title_params}', fontsize=16)
        
        plt.xlabel('Re', fontsize=12)
        plt.ylabel('Im', fontsize=12)
        plt.axhline(0, color='grey', linewidth=0.5); plt.axvline(0, color='grey', linewidth=0.5)
        plt.grid(True, linestyle='--', alpha=0.5); 
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            plt.legend()
            
        plt.gca().set_aspect('equal', adjustable='box')
        plt.margins(0.1) 
        
        save_filename = f"{base_name}.png"
        plt.savefig(save_filename, dpi=150, bbox_inches='tight') 
        plt.close()
        print(f"Saved separate plot to {save_filename}")

def create_labeled_combined_plot(filenames):

    plt.figure(figsize=(14, 14))
    ax = plt.gca()
    color_cycle = plt.cm.viridis(np.linspace(0, 0.9, len(filenames))) 
    plotted_labels = set()

    all_points_df = pd.DataFrame()

    for i, filename in enumerate(filenames):
        try:
            df = pd.read_csv(filename)
        except pd.errors.EmptyDataError:
            print(f"Warning: CSV file '{filename}' is empty. Skipping.")
            continue
        except FileNotFoundError:
            print(f"Warning: CSV file '{filename}' not found. Skipping.")
            continue
            
        df['real'] = df['real'].apply(format_rational)
        df['imag'] = df['imag'].apply(format_rational)
        df['color'] = [color_cycle[i]] * len(df)
        gen_label_match = re.search(r'(f\d+_to_f\d+)', filename)
        gen_label = gen_label_match.group(1).replace('f','F').replace('_to_',' To ') if gen_label_match else f"Gen {i+1}"
        df['gen_label'] = gen_label
        all_points_df = pd.concat([all_points_df, df])

        for family_id, group in df.groupby('family_id'):
             # Ensure there are enough points for parents and children
            if len(group) < 4: continue
            parents = group[group['point_type'] == 'parent']
            children = group[group['point_type'] == 'child']
            if len(parents) < 2 or len(children) < 2: continue
            
            p1, p2 = parents.iloc[0], parents.iloc[1]
            c1, c2 = children.iloc[0], children.iloc[1]

            # Draw the RECTANGLE
            rect_x = [p1['real'], c2['real'], p2['real'], c1['real'], p1['real']]
            rect_y = [p1['imag'], c2['imag'], p2['imag'], c1['imag'], p1['imag']]
            
            # Use generation label for legend
            current_label = df['gen_label'].iloc[0]
            if current_label not in plotted_labels:
                ax.plot(rect_x, rect_y, color=color_cycle[i], alpha=0.6, linewidth=1.2, label=current_label, zorder=i) # Use zorder for layering
                plotted_labels.add(current_label)
            else:
                ax.plot(rect_x, rect_y, color=color_cycle[i], alpha=0.6, linewidth=1.2, zorder=i)

    # Plot all points on top of the lines, maybe make them slightly transparent
    ax.scatter(all_points_df['real'], all_points_df['imag'], c=all_points_df['color'], s=40, zorder=len(filenames)+1, edgecolors='black', linewidth=0.5, alpha=0.9)
        
    # Calculate dynamic offsets based on final plot range
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_offset = (xlim[1] - xlim[0]) * 0.01 if (xlim[1] - xlim[0]) != 0 else 0.1
    y_offset = (ylim[1] - ylim[0]) * 0.01 if (ylim[1] - ylim[0]) != 0 else 0.1

    # Keep track of label positions to avoid overlap (basic approach)
    label_positions = []
    min_dist_sq = ((xlim[1]-xlim[0])*0.03)**2 + ((ylim[1]-ylim[0])*0.03)**2 # Min distance to draw label

    for _, point in all_points_df.iterrows():
        label = f"({point['real']:.0f}, {point['imag']:.0f}i)"
        pos = (point['real'], point['imag'] - y_offset)
        
        # Check if too close to existing labels
        too_close = False
        for lx, ly in label_positions:
            dist_sq = (pos[0]-lx)**2 + (pos[1]-ly)**2
            if dist_sq < min_dist_sq:
                too_close = True
                break
        
        if not too_close:
            ax.text(pos[0], pos[1], label, fontsize=6, ha='center', color='#111', zorder=len(filenames)+2) # Make labels smaller
            label_positions.append(pos)


    # Final plot styling
    # Extract N and R from the first filename for the title
    match = re.search(r'N(\d+)_R(\d+)', filenames[0])
    title_params = f" (N={match.group(1)}, R={match.group(2)})" if match else ""
    plt.title(f'Combined {title_params}', fontsize=18)
    
    plt.xlabel('Re', fontsize=14); plt.ylabel('Im', fontsize=14)
    plt.axhline(0, color='grey', linewidth=0.5, zorder=0); plt.axvline(0, color='grey', linewidth=0.5, zorder=0)
    plt.grid(True, linestyle='--', alpha=0.5); 
    plt.legend(title="Transition", fontsize=10, loc='upper left') # Adjust legend position
    plt.gca().set_aspect('equal', adjustable='box')
    plt.margins(0.1) # Add margins

    save_filename = f"combined_plot_N{match.group(1)}_R{match.group(2)}.png" if match else "labeled_rectangles_combined_plot.png"
    plt.savefig(save_filename, dpi=150, bbox_inches='tight') # Use tight layout
    plt.close()
    print(f"Saved combined plot to {save_filename}")


if __name__ == '__main__':
    # Find all CSV files that match the naming convention
    csv_files = sorted(glob.glob('rectangles_N*_R*_f*_to_f*.csv'))
    if not csv_files:
        print("No CSV data files found. Please run the C++ simulation first.")
        # Fallback for old file names
        csv_files_old = sorted(glob.glob('rectangles_f*.csv'))
        if not csv_files_old:
             print("Also could not find old-format 'rectangles_f*.csv' files.")
             exit()
        else:
             print("Found old-format CSV files. Using those...")
             csv_files = csv_files_old
            
    if not csv_files:
         exit() # Exit if still no files found

    try:
        import numpy as np
    except ImportError:
        print("Error: numpy is required for this script. Please install it using 'pip install numpy'")
        exit()
        
    print(f"Found {len(csv_files)} data files. Generating plots...")
    create_labeled_separate_plots(csv_files) # Generate separate plots first
    create_labeled_combined_plot(csv_files) # Then generate the combined plot
    print("\nAll plots have been saved as PNG files in your project folder.")

