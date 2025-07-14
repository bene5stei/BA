import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string
from matplotlib.patches import Rectangle
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from ipywidgets import interact, widgets
import os
from scipy.interpolate import interp1d, griddata
from scipy.signal import savgol_filter
from scipy.interpolate import make_interp_spline
depth_bins = [(0, 100), (100, 700)]
# Hilfsfunktion zum Erzeugen eines Polarplots
def plot_polar_subplot(ax, df, title):
    ratio_cols = [col for col in df.columns if col.startswith('ratio_')]
    freqs = [float(col.split('_')[1][:-2]) for col in ratio_cols]
    angles = np.deg2rad(df['backazimuth'].values)

    for freq, col in zip(freqs, ratio_cols):
        r = np.full_like(angles, freq)
        values = df[col].values
        sc = ax.scatter(angles, r, c=values, cmap='viridis', alpha=0.75, s=10)

    ax.set_rlabel_position(170)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_title(title)
   # Frequenzen sortieren und nur jedes 10. Label auswählen
    all_freqs = sorted(set(freqs))
    label_freqs = all_freqs[5::8]  # jedes 10. Element
    
    # Nur diese anzeigen und beschriften
    ax.set_yticks(label_freqs)
    ax.set_yticklabels([f"{f:.1f} Hz" for f in label_freqs], fontsize=9)
    return sc

    
def plot_component(ax, df, title):
    ratio_cols = [col for col in df.columns if col.startswith("ratio_")]
    frequencies = [float(col.split("_")[1][:-2]) for col in ratio_cols]

    # Plot nach Tiefenbereichen
    for (dmin, dmax), label in zip(depth_bins, depth_labels):
        df_subset = df[(df["depth"] >= dmin) & (df["depth"] < dmax)]
        num_events = len(df_subset)
        if df_subset.empty:
            continue

        mean_ratios = df_subset[ratio_cols].mean()
        sem_ratios = df_subset[ratio_cols].sem()
        legend_label = f"depth: {label} (n = {num_events})"

        if dmin == 0 and dmax == 700:
            ax.plot(frequencies, mean_ratios.values, linestyle='--', linewidth=2, label=legend_label)
        else:
            ax.plot(frequencies, mean_ratios.values, linestyle='-', label=legend_label)

        ax.fill_between(frequencies,
                        mean_ratios.values - sem_ratios.values,
                        mean_ratios.values + sem_ratios.values,
                        alpha=0.2)

    # Gesamtmittelwert über alle Tiefen (ohne SEM)
    overall_mean = df[ratio_cols].mean()
    ax.plot(frequencies,
            overall_mean.values,
            color='black',
            linestyle=':',
            linewidth=2,
            label='Overall mean (all depths)')

    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Mean ratio (FUR/WET)")
    ax.set_title(title)
    ax.grid(True)
    ax.legend()

def add_break_marks(ax_left, ax_right):
    d = .015
    kwargs = dict(transform=ax_left.transAxes, color='k', clip_on=False)
    ax_left.spines['right'].set_visible(False)
    ax_right.spines['left'].set_visible(False)
    ax_left.tick_params(labelright=False)
    ax_right.yaxis.tick_right()
    ax_left.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax_left.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
    kwargs.update(transform=ax_right.transAxes)
    ax_right.plot((-d, +d), (-d, +d), **kwargs)
    ax_right.plot((-d, +d), (1 - d, 1 + d), **kwargs)    


# Funktion zum Erkennen von Ausreißern via IQR
def detect_outliers_iqr(data):
    Q1 = np.percentile(data, 25)
    Q3 = np.percentile(data, 75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return (data < lower_bound) | (data > upper_bound)

# Funktion zur Vorbereitung der Daten
def prepare_plot_data(df):
    ratio_cols = [col for col in df.columns if col.startswith("ratio_")]
    frequencies = [float(col.split("_")[1][:-2]) for col in ratio_cols]
    return ratio_cols, frequencies

# Hilfsfunktion zum Zeichnen eines Rechtecks
def add_region(ax, lon_min, lon_max, lat_min, lat_max, label, color='red', alpha=0.2):
    rect = Rectangle((lon_min, lat_min), lon_max - lon_min, lat_max - lat_min,
                     linewidth=1, edgecolor=color, facecolor=color, alpha=alpha,
                     transform=ccrs.PlateCarree(), label=label)
    ax.add_patch(rect)

def smooth_line(x, y, points=500):
    x = np.array(x)
    y = np.array(y)
    x_new = np.linspace(x.min(), x.max(), points)
    spline = make_interp_spline(x, y, k=4)  # kubische Interpolation
    y_smooth = spline(x_new)
    return x_new, y_smooth

# Funktion zur Erstellung der Interpolationsdaten
def prepare_interpolation_data(df):
    ratio_cols = [col for col in df.columns if col.startswith('ratio_')]
    freqs = [float(col.split('_')[1][:-2]) for col in ratio_cols]

    theta_list = []
    r_list = []
    z_list = []

    for freq, col in zip(freqs, ratio_cols):
        theta = np.deg2rad(df['backazimuth'].values)
        r = np.full_like(theta, freq)
        z = df[col].values

        theta_list.extend(theta)
        r_list.extend(r)
        z_list.extend(z)

    theta_array = np.array(theta_list)
    r_array = np.array(r_list)
    z_array = np.array(z_list)

    # Periodizität sicherstellen
    theta_array_extended = np.concatenate([theta_array, theta_array + 2 * np.pi])
    r_array_extended = np.concatenate([r_array, r_array])
    z_array_extended = np.concatenate([z_array, z_array])

    return theta_array_extended, r_array_extended, z_array_extended, min(freqs), max(freqs)

def plot_component_with_sem(ax, df, color, label):
    # Frequenzen extrahieren
    ratio_cols = [col for col in df.columns if col.startswith("ratio_")]
    frequencies = [float(col.split("_")[1][:-2]) for col in ratio_cols]

    # Berechnungen
    mean_ratios = df[ratio_cols].mean()
    sem_ratios = df[ratio_cols].sem()
    std_ratios = df[ratio_cols].std()
    num_events = len(df)

    # Originale Werte als Arrays
    freq_array = np.array(frequencies)
    mean_array = mean_ratios.values
    sem_array = sem_ratios.values
    std_array = std_ratios.values

    # Glätten
    x_smooth, mean_smooth = smooth_line(freq_array, mean_array)
    _, sem_smooth = smooth_line(freq_array, sem_array)
    _, std_smooth = smooth_line(freq_array, std_array)

    # Plot: Glatte Linie
    ax.plot(x_smooth, mean_smooth, color=color, linewidth=2, label=f"{label} (n = {num_events})")
    ax.fill_between(x_smooth,
                    mean_smooth - sem_smooth,
                    mean_smooth + sem_smooth,
                    color=color,
                    alpha=0.3,
                    label=f"{label} ± SEM")