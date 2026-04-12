import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from pathlib import Path
import seaborn as sns
from typing import Dict, List, Tuple
from collections import defaultdict

def read_single_file(file_path: str,
                     type_label: str) -> pd.DataFrame:
    """
        read_single_file:
            Read single-site files, parse data, and tag types
        Parameters:
            file_path: [str]
                Input file path.
            type_label: [str]
                Wild or mutant type name.
        Returns:
            pd.DataFrame:["locus", "adjP_value", "log2FC", "avg_Abundance", "sig"]
    """
    data = []
    file_path = Path(file_path)
    if not file_path.exists():
        print(f"[Warning] File does not exist, skipping : {file_path}")
        return pd.DataFrame(columns=["locus", "adjP_value", "log2FC", "avg_Abundance", "sig"])

    try:
        with open(file_path, "r") as f:
            next(f)
            lines = f.readlines()

        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                # Parse core fields and handle non-numeric exceptions
                locus = parts[0]
                p_value = float(parts[1])
                fold_change = float(parts[2])
                avg_abundance = float(parts[3])
                data.append({
                    'locus': locus,
                    'adjP_value': p_value,
                    'log2FC': fold_change,
                    'avg_Abundance': avg_abundance,
                    'sig': type_label
                })
    except Exception as e:
        print(f"[Warning] File does not exist, skipping: {file_path}, error: {str(e)}")
        return pd.DataFrame(columns=["locus",
                                     "adjP_value",
                                     "log2FC",
                                     "avg_Abundance",
                                     "sig"])

    df = pd.DataFrame(data)
    return df

# Define batch read function: Iterate through file list and merge data of the same type.
def batch_read_files(file_list: List[str],
                     type_label: str) -> pd.DataFrame:
    """
        batch_read_files:
            Batch-read a list of files of the same type and merge them into a single DataFrame.
        Parameters:
            file_list: List[str]
                List of file paths requiring batch processing.
            type_label: [str]
                Wild or mutant type name.
        Returns:
            pd.DataFrame
    """
    df_list = []
    for file in file_list:
        df = read_single_file(file,
                              type_label)
        if not df.empty:
            df_list.append(df)
    return pd.concat(df_list, ignore_index=True) \
        if df_list else pd.DataFrame()

# ============= Volcano Map =============
def Volcano_Plot(up_files: List[str],
                 up2fold_files: List[str],
                 down_files: List[str],
                 down2fold_files: List[str],
                 mut_type: str,
                 outdir: Path):
    """
        Volcano_Plot:
                Visualize results using volcano plots to quickly filter and reveal
            the quantitative distribution of differentially expressed clusters.
        Parameters:
            up_files: List[str]
                Pass in the list of files resulting from the 1-fold up.
            up2fold_files: List[str]
                Pass in the list of files resulting from the 2-fold up.
            down_files: List[str]
                Pass in the list of files resulting from the 1-fold down.
            down2fold_files: List[str]
                Pass in the list of files resulting from the 2-fold down.
            mut_type: str
                Pass the mutant type.
            outdir: Path
                Output path for the result file.
        Returns:
            volcanoPlot_wt_vs_{mut_type}.png
    """
    # 1. Read the four files separately and mark their corresponding types.
    df_up = batch_read_files(up_files, "Up")  # 1x Up
    df_up2fold = batch_read_files(up2fold_files, "Up2Fold")  # 2x Up
    df_down = batch_read_files(down_files, "Down")  # 1x Down
    df_down2fold = batch_read_files(down2fold_files, "Down2Fold")  # 2x Down

    # 2. Merge all valid data
    all_dfs = [df for df in [df_up,
                             df_up2fold,
                             df_down,
                             df_down2fold]
               if not df.empty]
    if not all_dfs:
        print("[Warning] No valid data available. Volcano chart generation terminated.")
        return
    df_combined = pd.concat(all_dfs, ignore_index=True)

    epsilon = 1e-10
    # Replace values where adjP_value equals 0 with epsilon.
    df_combined['adjP_value'] = df_combined['adjP_value'].replace(0, epsilon)
    # 3. Calculate core metrics for volcano plots (log2FC and -log10p)
    df_combined['-log10p'] = -np.log10(df_combined['adjP_value'])

    # 4. Define the color palette
    palette = {
        'Up': '#FC8D59',
        'Up2Fold': '#D73027',
        'Down': '#4575B4',
        'Down2Fold': '#313695'
    }

    # 5. Plot a Volcano Diagram
    fig, ax = plt.subplots(figsize=(6, 5))

    # Plot Scatter Plots by Type
    for sig_type in palette.keys():
        subset = df_combined[df_combined['sig'] == sig_type]
        if not subset.empty:
            ax.scatter(subset['log2FC'],
                       subset['-log10p'],
                       c=palette[sig_type],
                       s=15,
                       alpha=0.8,
                       edgecolor='none',
                       label=sig_type)

    # Add threshold line
    ax.axhline(-np.log10(0.05),
               color='gray',
               ls='--',
               lw=0.8)  # Adjust P-values
    ax.axvline(1,
               color='#969696',
               ls='--',
               lw=0.8)  # log2(2)=1
    ax.axvline(-1,
               color='#969696',
               ls='--',
               lw=0.8)  # log2(0.5)=-1
    ax.axvline(0,
               color='gray',
               ls='-',
               lw=0.5)  # log2(1)=0

    # Set axes, titles, and styles
    ax.set_xlabel('log2(Fold Change)',
                  fontsize=11,
                  fontweight='bold')
    ax.set_ylabel('-log10(adjust P-value)',
                  fontsize=11,
                  fontweight='bold')
    ax.set_title(f'Volcano Plot of Differentially Expressed Small RNA Clusters - {mut_type}',
                 fontsize=12,
                 pad=12)
    ax.grid(alpha=0.3, ls=':')
    sns.despine()

    # Add legend
    ax.legend(title='Small RNA Status',
              loc='upper right',
              frameon=False,
              fontsize=6)

    save_path = outdir / f"volcanoPlot_wt_vs_{mut_type}.png"

    fig.savefig(save_path,
                dpi=600,
                bbox_inches='tight',
                transparent=False)
    plt.close(fig)
    print(f"The volcano diagram has been saved as {save_path.name}.")

    return str(save_path)

# ================ Bar-shaped Differential Distribution Diagram =====================
def BarshapedDistribution(up_files: List[str],
                          up2fold_files: List[str],
                          down_files: List[str],
                          down2fold_files: List[str],
                          mut_type: str,
                          outdir: Path):
    """
        BarshapedDistribution:
                Visualize using a bar-shaped differential analysis plot to reflect relative
            position information and the degree of significant difference at each site.
            All chromosomes are displayed in one figure.
        Parameters:
            up_files: List[str]
                Pass in the list of files resulting from the 1-fold up.
            up2fold_files: List[str]
                Pass in the list of files resulting from the 2-fold up.
            down_files: List[str]
                Pass in the list of files resulting from the 1-fold down.
            down2fold_files: List[str]
                Pass in the list of files resulting from the 2-fold down.
            mut_type: str
                Pass the mutant type.
            outdir: Path
                Output path for the result file.
        Returns:
            Barshaped_wt_vs_{mut_type}.png
    """

    # Read data
    df_up = batch_read_files(up_files, "Up")
    df_up2fold = batch_read_files(up2fold_files, "Up2Fold")
    df_down = batch_read_files(down_files, "Down")
    df_down2fold = batch_read_files(down2fold_files, "Down2Fold")

    # Merge all data
    all_dfs = [df for df in [df_up,
                             df_up2fold,
                             df_down,
                             df_down2fold]
               if not df.empty]
    if not all_dfs:
        print("[Warning] No valid data available. Bar chart generation terminated.")
        return
    df_combined = pd.concat(all_dfs, ignore_index=True)

    # Analyzing Chromosomal Location Information
    locus_split = df_combined['locus'].str.extract(r'^([^:]+):(\d+)-(\d+)$')
    df_combined[['chrom', 'start', 'stop']] = locus_split
    df_combined = df_combined.dropna(subset=['chrom',
                                             'start',
                                             'stop'])
    if df_combined.empty:
        print("[Warning] No valid chromosome position data available. Terminating plotting.")
        return

    # Calculate the midpoint
    df_combined['start'] = df_combined['start'].astype(int)
    df_combined['stop'] = df_combined['stop'].astype(int)
    df_combined['mid'] = (df_combined['start'] + df_combined['stop']) / 2

    # Define color mapping
    colors = {
        'Up2Fold': '#D73027',
        'Up': '#FC8D59',
        'Down': '#4575B4',
        'Down2Fold': '#313695'
    }
    # Grouping and sorting by chromosome
    chr_groups = defaultdict(list)
    for _, row in df_combined.iterrows():
        chr_groups[row['chrom']].append({
            'chr': row['chrom'],
            'mid': row['mid'],
            'log2FC': row['log2FC'],
            'direction': row['sig'],
            'norm_log2FC': abs(row['log2FC'])  # Normalization using absolute values
        })

    chrs = sorted(df_combined['chrom'].unique(),
                  key=lambda x: (not x.lstrip('Chr').isdigit(), x))
    n_chr = len(chrs)

    if n_chr == 0:
        print("[Warning] No valid chromosome data")
        return

    # Normalization
    all_log2fc = [abs(item['log2FC']) for chrom_list in chr_groups.values()
                  for item in chrom_list]
    max_val = max(all_log2fc) \
        if all_log2fc else 1.0

    for chrom_list in chr_groups.values():
        for item in chrom_list:
            item['norm_log2FC'] = item['norm_log2FC'] / max_val

    fig_height = 4.0 * n_chr
    fig, axes = plt.subplots(n_chr,
                             1,
                             figsize=(14, fig_height),
                             sharey=True)

    if n_chr == 1:
        axes = [axes]

    for ax, chr_name in zip(axes, chrs):
        group = chr_groups[chr_name]
        group = sorted(group,
                       key=lambda x: x["mid"])

        mids = [g["mid"] for g in group]
        heights = []
        plot_colors = []

        # Prepare drawing data
        for g in group:
            norm_val = g["norm_log2FC"]
            if g["direction"] in ["Up", "Up2Fold"]:
                heights.append(norm_val)
            else:
                heights.append(-norm_val)

            if g["direction"] == "Up2Fold":
                plot_colors.append(colors["Up2Fold"])
            elif g["direction"] == "Up":
                plot_colors.append(colors["Up"])
            elif g["direction"] == "Down2Fold":
                plot_colors.append(colors["Down2Fold"])
            else:  # "Down"
                plot_colors.append(colors["Down"])

        # Set column width
        if len(mids) > 1:
            span = max(mids) - min(mids)
            width = span * 0.004 \
                if span > 0 else 1.0
        else:
            width = 1.0

        # Plot a bar chart
        ax.bar(mids,
               heights,
               width=width,
               color=plot_colors,
               edgecolor="black",
               linewidth=0.3,
               alpha=0.8)

        # Set axes and labels
        ax.set_xlabel(f"{chr_name}",
                      labelpad=10,
                      fontweight='bold')
        ax.set_ylabel("Normalized |log2(FC)|",
                      fontweight='bold')

        # Add Grid and Guides
        ax.grid(axis='y',
                linestyle='--',
                alpha=0.3)
        ax.axhline(y=0,
                   color='black',
                   linewidth=1.0)

        # Set the x-axis range to the data range
        if len(mids) > 0:
            ax.set_xlim(min(mids) * 0.99,
                        max(mids) * 1.01)

    fig.suptitle(f'Bar-shaped Distribution of sRNA Clusters - {mut_type}',
                 fontsize=16,
                 fontweight='bold',
                 y=0.98)

    legend_handles = [
        Patch(color="#D73027",
              label="Up2-regulated",
              alpha=0.8),
        Patch(color="#FC8D59",
              label="Up-regulated",
              alpha=0.8),
        Patch(color="#4575B4",
              label="Down-regulated",
              alpha=0.8),
        Patch(color="#313695",
              label="Down2-regulated",
              alpha=0.8)
    ]

    fig.legend(handles=legend_handles,
               loc="upper right",
               bbox_to_anchor=(0.98, 0.98),
               frameon=True,
               fancybox=True,
               shadow=True,
               fontsize=10)

    # Adjust the layout
    plt.tight_layout(rect=(0, 0, 0.95, 0.96))

    # Save png
    save_path = outdir / f"Barshaped_wt_vs_{mut_type}.png"
    fig.savefig(save_path,
                dpi=300,
                bbox_inches='tight',
                facecolor='white')
    plt.close(fig)

    print(f"The bar distribution diagram has been saved as {save_path.name}.")

    return str(save_path)

# ================= Sliding Windows ==================
def Slide_Windows(up_files: List[str],
                  up2fold_files: List[str],
                  down_files: List[str],
                  down2fold_files: List[str],
                  mut_type: str,
                  outdir: Path):
    """
        Slide_Windows:
            Implement sliding window visualization.
        Parameters:
            up_files: List[str]
                Pass in the list of files resulting from the 1-fold up.
            up2fold_files: List[str]
                Pass in the list of files resulting from the 2-fold up.
            down_files: List[str]
                Pass in the list of files resulting from the 1-fold down.
            down2fold_files: List[str]
                Pass in the list of files resulting from the 2-fold down.
            mut_type: str
                Pass the mutant type.
            outdir: Path
                Output path for the result file.
        Returns:
            SlidingWindow_{mut_type}.png
    """
    # 1. Batch read four file types and mark their corresponding types.
    df_up = batch_read_files(up_files, "Up")
    df_up2fold = batch_read_files(up2fold_files, "Up2Fold")
    df_down = batch_read_files(down_files, "Down")
    df_down2fold = batch_read_files(down2fold_files, "Down2Fold")

    # 2. Merge all data
    all_dfs = [df for df in [df_up, df_up2fold, df_down, df_down2fold]
               if not df.empty]
    if not all_dfs:
        print("[Warning] No valid differential data detected; sliding window chart plotting terminated.")
        return None
    df_combined = pd.concat(all_dfs, ignore_index=True)

    # 3. Analyzing Chromosomal Location Information
    locus_split = df_combined['locus'].str.extract(r'^([^:]+):(\d+)-(\d+)$')
    df_combined[['chrom', 'start', 'stop']] = locus_split
    df_combined = df_combined.dropna(subset=['chrom', 'start', 'stop'])
    if df_combined.empty:
        print("[Warning] No valid chromosome position data available. Terminating plotting.")
        return None

    df_combined['start'] = df_combined['start'].astype(int)
    df_combined['stop'] = df_combined['stop'].astype(int)

    colors = {
        'Up2Fold': '#D73027',
        'Up': '#FC8D59',
        'Down': '#4575B4',
        'Down2Fold': '#313695',
    }

    # 4. Sliding Window Parameters
    WINDOW_SIZE = 1.0  # Window size = 1Mbp
    STEP_SIZE = 0.1  # step = 0.1Mbp

    # 5. Chromosome Sorting Function
    def chr_sort_key(chr_name):
        name = chr_name.replace("chr", "").replace("CHR", "").replace("Chr", "")
        try:
            return 0, int(name)
        except ValueError:
            return 1, name

    # 6. Grouping and sorting by chromosome
    chromosomes = sorted(df_combined['chrom'].unique(),
                         key=chr_sort_key)
    n_chr = len(chromosomes)

    if n_chr == 0:
        print("[Warning] No valid chromosome data")
        return

    # 7. Compute the global maximum value across all chromosomes.
    global_max = 0

    # Iterate through all chromosomes, compute the window count for each chromosome, and find the global maximum value.
    for chrom in chromosomes:
        chrom_data = df_combined[df_combined['chrom'] == chrom].copy()
        if chrom_data.empty:
            continue

        # Calculate chromosome length and the midpoint position within sRNA clusters
        chrom_length = chrom_data['stop'].max()
        chrom_length_mb = chrom_length / 1e6
        chrom_data['midpoint'] = (chrom_data['start'] + chrom_data['stop']) / 2 / 1e6

        # Generate sliding window
        window_starts = np.arange(0, chrom_length_mb, STEP_SIZE)
        window_ends = window_starts + WINDOW_SIZE

        # Count the number of windows on the current chromosome
        for start, end in zip(window_starts, window_ends):
            window_mask = (chrom_data['midpoint'] >= start) & (chrom_data['midpoint'] < end)
            window_data = chrom_data[window_mask]

            # Calculate the total count for the current window
            window_max = max(len(window_data[window_data['sig'] == 'Up']),
                            len(window_data[window_data['sig'] == 'Up2Fold']),
                            len(window_data[window_data['sig'] == 'Down']),
                            len(window_data[window_data['sig'] == 'Down2Fold']))

            # Update the global maximum value
            if window_max > global_max:
                global_max = window_max

    # If all windows are 0, set global_max to 1 to avoid division by zero.
    if global_max == 0:
        global_max = 1

    fig_height = 4.0 * n_chr
    fig, axes = plt.subplots(n_chr,
                             1,
                             figsize=(14, fig_height),
                             sharey=True)

    if n_chr == 1:
        axes = [axes]

    # 8. Calculate sliding window density for each chromosome (using global normalization)
    for ax, chrom in zip(axes, chromosomes):
        chrom_data = df_combined[df_combined['chrom'] == chrom].copy()
        if chrom_data.empty:
            continue

        # Calculate chromosome length and the midpoint position within sRNA clusters
        chrom_length = chrom_data['stop'].max()
        chrom_length_mb = chrom_length / 1e6
        chrom_data['midpoint'] = (chrom_data['start'] + chrom_data['stop']) / 2 / 1e6

        # Generate sliding window
        window_starts = np.arange(0, chrom_length_mb, STEP_SIZE)
        window_ends = window_starts + WINDOW_SIZE
        # The midpoint of the window serves as the X-axis.
        x_positions = (window_starts + window_ends) / 2

        # Initialize count lists for each type
        up_counts, up2_counts, down_counts, down2_counts = [], [], [], []

        # Count the number of clusters of each type within each window.
        for start, end in zip(window_starts, window_ends):
            window_mask = (chrom_data['midpoint'] >= start) & (chrom_data['midpoint'] < end)
            window_data = chrom_data[window_mask]

            up_counts.append(len(window_data[window_data['sig'] == 'Up']))
            up2_counts.append(len(window_data[window_data['sig'] == 'Up2Fold']))
            down_counts.append(len(window_data[window_data['sig'] == 'Down']))
            down2_counts.append(len(window_data[window_data['sig'] == 'Down2Fold']))

        # Normalize using the global maximum value (scale down to negative values)
        up_density = [count / global_max
                      for count in up_counts]
        up2_density = [count / global_max
                       for count in up2_counts]
        down_density = [-count / global_max
                        for count in down_counts]
        down2_density = [-count / global_max
                         for count in down2_counts]

        # Plot the sliding window curve
        ax.plot(x_positions,
                up_density,
                color=colors['Up'],
                linewidth=1.2,
                alpha=0.8,
                label='Up(1x)')
        ax.plot(x_positions,
                up2_density,
                color=colors['Up2Fold'],
                linewidth=1.5,
                alpha=0.9,
                label='Up(2x)')
        ax.plot(x_positions,
                down_density,
                color=colors['Down'],
                linewidth=1.2,
                alpha=0.8,
                label='Down(1x)')
        ax.plot(x_positions,
                down2_density,
                color=colors['Down2Fold'],
                linewidth=1.5,
                alpha=0.9,
                label='Down(2x)')

        # Set axes and labels
        ax.set_xlabel(f"{chrom}",
                      labelpad=10,
                      fontweight='bold')
        ax.set_ylabel("Normalized Density",
                      fontweight='bold')

        # Add guides and grids
        ax.axhline(y=0,
                   color='black',
                   linewidth=1.0,
                   alpha=0.8)
        ax.grid(axis='y',
                linestyle='--',
                alpha=0.3)

        # Set the X-axis range
        if len(x_positions) > 0:
            ax.set_xlim(0,
                        max(x_positions) * 1.02)

        ax.set_ylim(-1.1, 1.1)

    fig.suptitle(f'Sliding Window Analysis of sRNA Clusters - {mut_type}',
                 fontsize=16,
                 fontweight='bold',
                 y=0.98)

    legend_handles = [
        Patch(color="#D73027", label="Up2-regulated", alpha=0.8),
        Patch(color="#FC8D59", label="Up-regulated", alpha=0.8),
        Patch(color="#4575B4", label="Down-regulated", alpha=0.8),
        Patch(color="#313695", label="Down2-regulated", alpha=0.8)
    ]

    fig.legend(handles=legend_handles,
               loc="upper right",
               bbox_to_anchor=(0.98, 0.98),
               frameon=True,
               fancybox=True,
               shadow=True,
               fontsize=10,
               title="Significance Category",
               title_fontsize=11)

    plt.tight_layout(rect=(0, 0, 0.95, 0.96))

    save_path = outdir / f"SlidingWindow_{mut_type}.png"
    fig.savefig(save_path,
                dpi=300,
                bbox_inches='tight',
                facecolor='white')
    plt.close(fig)

    print(f"The sliding window diagram has been saved as {save_path.name}.")

    return str(save_path)

def Visualization(inputFiles, outdir: Path, mut_type: str):
    """
    Visualization:
        Provide a visual interface
    Parameters:
        inputFiles:
            Input file
        outdir:
            Output the result file to the specified path.
        mut_type:
            Mutant type name.
    Output:
        Visualization/{mut_type}
    """
    out_dir = outdir / f"{mut_type}"
    out_dir.mkdir(parents=True,
                  exist_ok=True)

    visualizationRst = []

    vRst = Volcano_Plot(inputFiles['up'],
                        inputFiles['up2fold'],
                        inputFiles['down'],
                        inputFiles['down2fold'],
                        mut_type,
                        out_dir)
    bRst = BarshapedDistribution(inputFiles['up'],
                                 inputFiles['up2fold'],
                                 inputFiles['down'],
                                 inputFiles['down2fold'],
                                 mut_type,
                                 out_dir)
    sRst = Slide_Windows(inputFiles['up'],
                         inputFiles['up2fold'],
                         inputFiles['down'],
                         inputFiles['down2fold'],
                         mut_type,
                         out_dir)

    visualizationRst.append(vRst)
    visualizationRst.append(bRst)
    visualizationRst.append(sRst)

    return visualizationRst