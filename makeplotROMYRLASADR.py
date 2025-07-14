def makeplotROMYRLASADR(st, config, oadr, aadr):
    import matplotlib.pyplot as plt
    # Liste der auszuschließenden Stationen oder Channels
    exclude_stations = []# ['FUR', 'WET' ]
    #exclude_channels = ['BJZ']
    st_in = st.copy()
    st_in.sort(keys=['channel'], reverse=True)

    # Vorbereitung: Gemeinsamer Plot
    fig, ax = plt.subplots(figsize=(15, 8))
    font = 10

    # Linke Daten (Originale) – alle gemeinsam
    for tr in st_in:
        if tr.stats.station in exclude_stations:
            continue
        ax.plot(
            tr.times(),
            tr.data,
            label=f"Original: {tr.stats.network}.{tr.stats.station}.{tr.stats.channel}",
            lw=1.0
        )

    # Rechte Daten (ADR) – alle gemeinsam
    adr_locations = ['22', '23']
    adr_station = 'ROMY'
    adr_network = 'BW'
    adr_channel = 'BJZ'

    for loc in adr_locations:
        if loc == '23':
            adr_trace = aadr.select(station=adr_station,
                                location=loc,
                                network=adr_network,
                                channel=adr_channel)
            #print(f"Trace-Auswahl für loc={loc}: {adr_trace}")

            if adr_trace:
                    tr_adr = adr_trace[0]
                    ax.plot(
                        tr_adr.times(),
                        tr_adr.data,
                        label=f"ADR: {adr_network}.{adr_station}.{loc}.{adr_channel}",
                        lw=1.0,
                        linestyle='--'  # zur Unterscheidung
                    )
        if loc == '22':
            adr_trace = oadr.select(station=adr_station,
                                location=loc,
                                network=adr_network,
                                channel=adr_channel)
            if adr_trace:
                    tr_adr = adr_trace[0]
                    ax.plot(
                            tr_adr.times(),
                            tr_adr.data,
                            label=f"ADR: {adr_network}.{adr_station}.{loc}.{adr_channel}",
                            lw=1.0,
                            linestyle='--'  # zur Unterscheidung
            )
        
    ax.set_xlabel(
        f"Time from {tr.stats.starttime.date} {str(tr.stats.starttime.time)[:8]} (s)",
        fontsize=font
    )
    ax.set_ylabel("Amplitude", fontsize=font)
    ax.legend(loc='upper right', fontsize=8)
    ax.set_title(f"Direkte Messung & ADR-Daten", fontsize=12)
    from datetime import datetime
    timestamp = config['tbeg'].strftime("%Y%m%d_%H%M%S")
    save_path = config['outpath_figs'] + f"gemeinsom_{timestamp}.png"
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Plot gespeichert unter: {save_path}")
    return fig
