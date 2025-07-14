def makeplotDifferenceROMYADRsubplotslag(st, iadr, oadr, aadr, config):
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    from datetime import datetime

    adr_locations = ['21','22', '23']
    adr_station = 'ROMY'
    adr_network = 'BW'
    adr_channel = 'BJZ'

    font = 10
    num_plots = len(adr_locations) * 2  # zwei Subplots pro Location
    fig, axes = plt.subplots(nrows=num_plots, figsize=(18, 4 * num_plots), sharex=True)

    if num_plots == 1:
        axes = [axes]

    # Hole Originalspur
    original_trace = None
    for tr in st:
        if tr.stats.station == adr_station and tr.stats.channel == adr_channel:
            original_trace = tr
            break

    if original_trace is None:
        raise ValueError("Originalspur ROMY.BJZ nicht gefunden")

    for i, loc in enumerate(adr_locations):
        if loc == '21':
           adr_trace = iadr.select(station=adr_station, location=loc, network=adr_network, channel=adr_channel)
        elif loc == '22':
            adr_trace = oadr.select(station=adr_station, location=loc, network=adr_network, channel=adr_channel)
        elif loc == '23':
            adr_trace = aadr.select(station=adr_station, location=loc, network=adr_network, channel=adr_channel)
        else:
            continue

        if not adr_trace:
            print(f"Keine Daten gefunden für Location {loc}")
            continue

        tr_adr = adr_trace[0]
        #um Zeiten nachträglich anpassen zukönnen
        
                # Synchronisieren (vorläufig trimmen für Cross-Korrelation)
        tr1_full, tr2_full = original_trace.copy(), tr_adr.copy()
        tr1_full.trim(starttime=max(tr1_full.stats.starttime, tr2_full.stats.starttime),
                      endtime=min(tr1_full.stats.endtime, tr2_full.stats.endtime))

        tr2_full.trim(starttime=tr1_full.stats.starttime, endtime=tr1_full.stats.endtime)

        if tr1_full.stats.sampling_rate != tr2_full.stats.sampling_rate:
            tr2_full.resample(tr1_full.stats.sampling_rate)

        # Länge angleichen
        minlen = min(len(tr1_full.data), len(tr2_full.data))
        tr1_full.data = tr1_full.data[:minlen]
        tr2_full.data = tr2_full.data[:minlen]

        # Kreuzkorrelation und optimalen Lag bestimmen
        tr1_norm = (tr1_full.data - np.mean(tr1_full.data)) / np.std(tr1_full.data)
        tr2_norm = (tr2_full.data - np.mean(tr2_full.data)) / np.std(tr2_full.data)
        corr = np.correlate(tr1_norm, tr2_norm, mode='full')
        lag = np.argmax(corr) - (len(tr1_norm) - 1)
        max_corr = np.max(corr) / len(tr1_norm)

        time_shift = lag / tr1_full.stats.sampling_rate
        print(f"Lag für {loc}: {lag} Samples → {time_shift:.3f} s")

        # ADR-Daten um den optimalen Lag verschieben
        tr2_shifted = tr_adr.copy()
        tr2_shifted.stats.starttime += time_shift

        # Neu trimmen nach Zeitverschiebung
        tr1, tr2 = original_trace.copy(), tr2_shifted
        tr1.trim(starttime=max(tr1.stats.starttime+, tr2.stats.starttime),
                 endtime=min(tr1.stats.endtime, tr2.stats.endtime))
        tr2.trim(starttime=tr1.stats.starttime, endtime=tr1.stats.endtime)

        if tr1.stats.sampling_rate != tr2.stats.sampling_rate:
            tr2.resample(tr1.stats.sampling_rate)

        if len(tr1.data) != len(tr2.data):
            minlen = min(len(tr1.data), len(tr2.data))
            tr1.data = tr1.data[:minlen]
            tr2.data = tr2.data[:minlen]


        # Original vs ADR Plot
        ax_data = axes[i * 2]
        ax_data.plot(tr1.times(), tr1.data, label='dircet recording ROMY.BJZ', lw=1)
        ax_data.plot(tr2.times(), tr2.data, label=f'ADR ROMY.{loc}.BJZ', lw=1)
        ax_data.set_ylabel("rotation rate (rad/s)", fontsize=font)
        ax_data.legend(fontsize=8)
        ax_data.set_title(f"direct recordings of ROMY Z vs ADR: ROMY.{loc}.BJZ", fontsize=12)

        # Differenzplot
        diff = np.abs(tr1.data - tr2.data)
        max_diff = np.max(diff)
        ax_diff = axes[i * 2 + 1]
        # Kreuzkorrelation berechnen (z-normalisiert)
        tr1_norm = (tr1.data - np.mean(tr1.data)) / np.std(tr1.data)
        tr2_norm = (tr2.data - np.mean(tr2.data)) / np.std(tr2.data)
        corr = np.correlate(tr1_norm, tr2_norm, mode='valid')
        max_corr = np.max(corr) / len(tr1_norm)  # Normierung auf [−1, 1]
        ax_diff.plot(tr1.times(), diff, label=f'magnitude of the amplitude difference between ADR and ROMY.{loc}.BJZ', lw=1, color='red')
        ax_diff.set_ylabel("", fontsize=font)
        #ax_diff.legend(fontsize=8)
        ax_diff.set_title(f"magnitude of the amplitude difference between ADR and ROMY.{loc}.BJZ", fontsize=12)
        ax_diff.text(0.98, 0.90,f"Max diff: {max_diff:.2e}\nCC: {max_corr:.3f}\nLag: {time_shift:.3f} s",
                     transform=ax_diff.transAxes,
                     ha='right', va='top', fontsize=font, color='black',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    axes[-1].set_xlabel(f"Time (s) from {tr1.stats.starttime.strftime('%Y-%m-%d %H:%M:%S')}", fontsize=font)
    fig.suptitle("comparism of direct ROMY recordings and array derived rotations", fontsize=14)

    # Speicherpfad vorbereiten
    timestamp = config['tbeg'].strftime("%Y%m%d_%H%M%S")
    os.makedirs(config['outpath_figs'], exist_ok=True)
    save_path = os.path.join(config['outpath_figs'], f"difference_and_comparison_ROMY_BJZ_{timestamp}_lag_wi.png")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Plot gespeichert unter: {save_path}")

    return fig
def makeplotDifferenceROMYADRsubplots(st,iadr, oadr, aadr, config):
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    from datetime import datetime

    adr_locations = ['21', '22', '23']
    adr_station = 'ROMY'
    adr_network = 'BW'
    adr_channel = 'BJZ'

    font = 10
    num_plots = len(adr_locations) * 2  # zwei Subplots pro Location
    fig, axes = plt.subplots(nrows=num_plots, figsize=(18, 4 * num_plots), sharex=True)

    if num_plots == 1:
        axes = [axes]

    # Hole Originalspur
    original_trace = None
    for tr in st:
        if tr.stats.station == adr_station and tr.stats.channel == adr_channel:
            original_trace = tr
            break

    if original_trace is None:
        raise ValueError("Originalspur ROMY.BJZ nicht gefunden")

    for i, loc in enumerate(adr_locations):
        if loc == '21':
            adr_trace = iadr.select(station=adr_station, location=loc, network=adr_network, channel=adr_channel)
        elif loc == '22':
            adr_trace = oadr.select(station=adr_station, location=loc, network=adr_network, channel=adr_channel)
        elif loc == '23':
            adr_trace = aadr.select(station=adr_station, location=loc, network=adr_network, channel=adr_channel)
        else:
            continue

        if not adr_trace:
            print(f"Keine Daten gefunden für Location {loc}")
            continue

        tr_adr = adr_trace[0]
        #um Zeiten nachträglich anpassen zukönnen
        
        # Synchronisieren
        tr1, tr2 = original_trace.copy(), tr_adr.copy()
        tr1.trim(starttime=max(tr1.stats.starttime, tr2.stats.starttime),
                 endtime=min(tr1.stats.endtime, tr2.stats.endtime))
        tr2.trim(starttime=tr1.stats.starttime, endtime=tr1.stats.endtime)

        if tr1.stats.sampling_rate != tr2.stats.sampling_rate:
            tr2.resample(tr1.stats.sampling_rate)

        if len(tr1.data) != len(tr2.data):
            minlen = min(len(tr1.data), len(tr2.data))
            tr1.data = tr1.data[:minlen]
            tr2.data = tr2.data[:minlen]

        # Original vs ADR Plot
        ax_data = axes[i * 2]
        ax_data.plot(tr1.times(), tr1.data, label='dircet recording ROMY.BJZ', lw=1)
        ax_data.plot(tr2.times(), tr2.data, label=f'ADR ROMY.{loc}.BJZ', lw=1)
        ax_data.set_ylabel("rotation rate (rad/s)", fontsize=font)
        ax_data.legend(fontsize=8)
        ax_data.set_title(f"direct recordings of ROMY Z vs ADR: ROMY.{loc}.BJZ", fontsize=12)

        # Differenzplot
        diff = np.abs(tr1.data - tr2.data)
        max_diff = np.max(diff)
        ax_diff = axes[i * 2 + 1]
        # Kreuzkorrelation berechnen (z-normalisiert)
        tr1_norm = (tr1.data - np.mean(tr1.data)) / np.std(tr1.data)
        tr2_norm = (tr2.data - np.mean(tr2.data)) / np.std(tr2.data)
        corr = np.correlate(tr1_norm, tr2_norm, mode='valid')
        max_corr = np.max(corr) / len(tr1_norm)  # Normierung auf [−1, 1]
        ax_diff.plot(tr1.times(), diff, label=f'magnitude of the amplitude difference between ADR and ROMY.{loc}.BJZ', lw=1, color='red')
        ax_diff.set_ylabel("", fontsize=font)
        #ax_diff.legend(fontsize=8)
        ax_diff.set_title(f"magnitude of the amplitude difference between ADR and ROMY.{loc}.BJZ", fontsize=12)
        ax_diff.text(0.98, 0.90,f"Max diff: {max_diff:.2e}\nCC: {max_corr:.3f}",
                     transform=ax_diff.transAxes,
                     ha='right', va='top', fontsize=font, color='black',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    axes[-1].set_xlabel(f"Time (s) from {tr1.stats.starttime.strftime('%Y-%m-%d %H:%M:%S')}", fontsize=font)
    fig.suptitle("comparism of direct ROMY recordings and array derived rotations", fontsize=14)

    # Speicherpfad vorbereiten
    timestamp = config['tbeg'].strftime("%Y%m%d_%H%M%S")
    os.makedirs(config['outpath_figs'], exist_ok=True)
    save_path = os.path.join(config['outpath_figs'], f"difference_and_comparison_ROMY_BJZ_{timestamp}_wi.png")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Plot gespeichert unter: {save_path}")

    return fig
