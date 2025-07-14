def makeplotDifferenceROMYADR(st, oadr, aadr):
    import matplotlib.pyplot as plt
    import numpy as np

    adr_locations = ['23']#, '23']
    adr_station = 'ROMY'
    adr_network = 'BW'
    adr_channel = 'BJZ'

    # Gemeinsamer Plot für Differenz
    fig, ax = plt.subplots(figsize=(15, 6))
    font = 10

    # Hole Originalspur aus st (z. B. ROMY.BJZ)
    original_trace = None
    for tr in st:
        if tr.stats.station == adr_station and tr.stats.channel == adr_channel:
            original_trace = tr
            break

    if original_trace is None:
        raise ValueError("Originalspur ROMY.BJZ nicht gefunden")

    for loc in adr_locations:
        #if loc == '21':
         #   adr_trace = iadr.select(station=adr_station, location=loc, network=adr_network, channel=adr_channel)
        if loc == '22':
            adr_trace = oadr.select(station=adr_station, location=loc, network=adr_network, channel=adr_channel)
        elif loc == '23':
            adr_trace = aadr.select(station=adr_station, location=loc, network=adr_network, channel=adr_channel)
        else:
            continue

        if not adr_trace:
            continue

        tr_adr = adr_trace[0]

        # Synchronisation: bringe beide Datenreihen auf gemeinsame Zeitbasis
        tr1, tr2 = original_trace.copy(), tr_adr.copy()
        tr1.trim(starttime=max(tr1.stats.starttime, tr2.stats.starttime),
                 endtime=min(tr1.stats.endtime, tr2.stats.endtime))
        tr2.trim(starttime=tr1.stats.starttime, endtime=tr1.stats.endtime)

        # Resample falls nötig
        if tr1.stats.sampling_rate != tr2.stats.sampling_rate:
            tr2.resample(tr1.stats.sampling_rate)

        if len(tr1.data) != len(tr2.data):
            minlen = min(len(tr1.data), len(tr2.data))
            tr1.data = tr1.data[:minlen]
            tr2.data = tr2.data[:minlen]

        # Differenz berechnen
        #diff = tr1.data - tr2.data
        # Betrag der Differenz berechnen
        diff = np.abs(tr1.data - tr2.data)

        ax.plot(tr1.times(), diff, label=f'Differenz ROMY.{loc}.BJZ', lw=1)

    ax.set_xlabel("Zeit (s)", fontsize=font)
    ax.set_ylabel("Amplitude Differenz", fontsize=font)
    ax.legend(fontsize=8)
    ax.set_title("Betrag der Differenz zwischen Original und ADR-Daten", fontsize=12)
    from datetime import datetime
    timestamp = config['tbeg'].strftime("%Y%m%d_%H%M%S")
    save_path = config['outpath_figs']+f"difference_plot_ROMY_BJZ_{timestamp}.png"
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Plot gespeichert unter: {save_path}")


    return fig
