def __makeplotStreamSpectra2(st, config, fscale=None):

    from scipy import fftpack
    from andbro__fft import __fft
    import matplotlib.pyplot as plt

    st_in = st.copy()

    NN = len(st_in)
    rot_scaling, rot_unit = 1e9, r"nrad/s"
    trans_scaling, trans_unit = 1e6, r"$\mu$m/s"

    fig, axes = plt.subplots(NN,2,figsize=(15,int(NN*2)), sharex='col')

    font = 14

    plt.subplots_adjust(hspace=0.3)

    ## _______________________________________________

    st.sort(keys=['channel'], reverse=True)

    for i, tr in enumerate(st_in):

#         comp_fft = abs(fftpack.fft(tr.data))
#         ff       = fftpack.fftfreq(comp_fft.size, d=1/tr.stats.sampling_rate)
#         comp_fft = fftpack.fftshift(comp_fft)
#         ff, spec = ff[1:len(ff)//2], abs(fftpack.fft(tr.data)[1:len(ff)//2])

        if tr.stats.channel[-2] == "J":
            scaling = rot_scaling
        elif tr.stats.channel[-2] == "H":
            scaling = trans_scaling

        spec, ff, ph = __fft(tr.data*scaling, tr.stats.delta, window=None, normalize=None)


        ## _________________________________________________________________
        if tr.stats.channel[-2] == "J":
            axes[i,0].plot(
                        tr.times(),
                        tr.data*rot_scaling,
                        color='black',
                        label='{} {}'.format(tr.stats.station, tr.stats.channel),
                        lw=1.0,
                        )

        elif tr.stats.channel[-2] == "H":
            axes[i,0].plot(
                        tr.times(),
                        tr.data*trans_scaling,
                        color='black',
                        label='{} {}'.format(tr.stats.station, tr.stats.channel),
                        lw=1.0,
                        )
        ## _________________________________________________________________
        if fscale == "loglog":
            axes[i,1].loglog(ff, spec, color='black', lw=1.0)
        elif fscale == "loglin":
            axes[i,1].semilogx(ff, spec, color='black', lw=1.0)
        elif fscale == "linlog":
            axes[i,1].semilogy(ff, spec, color='black', lw=1.0)
        else:
            axes[i,1].plot(ff, spec, color='black', lw=1.0)         
        

        if tr.stats.channel[1] == "J":
            sym, unit = r"$\Omega$", rot_unit
        elif tr.stats.channel[1] == "H":
            sym, unit = "v", trans_unit
        else:
            unit = "Amplitude", "a.u."

        axes[i,0].set_ylabel(f'{sym} ({unit})',fontsize=font)    
        axes[i,1].set_ylabel(f'ASD \n({unit}/Hz)',fontsize=font)        
        axes[i,0].legend(loc='upper left',bbox_to_anchor=(0.8, 1.10), framealpha=1.0)

#         axes[i,0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
#         axes[i,1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    if "fmin" in config.keys() and "fmax" in config.keys():
        axes[i,1].set_xlim(config['fmin'],config['fmax'])

    axes[NN-1,0].set_xlabel(f"Time from {tr.stats.starttime.date} {str(tr.stats.starttime.time)[:8]} (s)",fontsize=font)     
    axes[NN-1,1].set_xlabel(f"Frequency (Hz)",fontsize=font)     

    del st_in
    return fig

#nur für translation --> nicht für rotation
def __makeplotspectraT(st, config, fscale=None):

    from scipy import fftpack
    from andbro__fft import __fft
    import matplotlib.pyplot as plt

    st_in = st.copy()
    print([tr.stats.channel for tr in st])
    st_in = st_in.select(component="T") 
    if len(st_in) == 0:
        print(" -> keine Transversal-Komponente gefunden.")
        return
    NN = len(st_in)+1

    trans_scaling = 1e6
    trans_unit = r"$\mathrm{\mu m/s}$"
   
    fig, axes = plt.subplots(NN,figsize=(18,NN*4), sharex=False)

    font = 14

    plt.subplots_adjust(hspace=0.3)

    ## _______________________________________________

    st.sort(keys=['channel'], reverse=True)

    for i, tr in enumerate(st_in):
        scaling = trans_scaling
        spec, ff, ph = __fft(tr.data*scaling, tr.stats.delta, window=None, normalize=None)
        
        if tr.stats.channel[-2] == "H":
            if i == 0:
                tr1 = tr
                spec1, ff1, ph1 = __fft(tr.data*scaling, tr.stats.delta, window=None, normalize=None)
            if i == 1:
                tr2 = tr
                spec2, ff2, ph2 = __fft(tr.data*scaling, tr.stats.delta, window=None, normalize=None)     
        
        if tr.stats.channel[1] == "H":
            sym, unit = "v", trans_unit
        else:
            unit = "Amplitude", "a.u."
    envelope1 = np.abs(hilbert(tr1.data * trans_scaling))
    envelope2 = np.abs(hilbert(tr2.data * trans_scaling))
    max_amp1 = np.max(envelope1)
    max_amp2 = np.max(envelope2)
    amp_ratioFURWET = max_amp1/max_amp2
    amp_ratioFURWET_g = np.round(amp_ratioFURWET, 5)
    # Farben definieren
    colors = {
        'tr1': '#1f77b4',        # Blau
        'envelope1': '#2ca02c',  # Grün
        'tr2': '#ff7f0e',        # Orange
        'envelope2': '#17becf',  # Türkis
        'spec1': '#1f77b4',      # Blau
        'spec2': '#ff7f0e',      # Orange
    }
    
    axes[0].plot(tr1.times(), tr1.data*trans_scaling, label=f"{tr1.stats.station} {tr1.stats.channel}", color=colors['tr1'], lw=1.0)
    axes[0].plot(tr1.times(), envelope1, color=colors['envelope1'], label='Envelope ' + tr1.stats.channel) 
    axes[0].set_ylabel(f'envelope and {sym} ({unit})', fontsize=font)
    axes[0].text(0.995, 0.15,f"Max amp (FUR): {max_amp1:.2e} $\mu$m/s", transform=axes[0].transAxes,
                   ha='right', va='top', fontsize=font, color='black')
    axes[1].plot(tr2.times(), tr2.data*trans_scaling, label=f"{tr2.stats.station} {tr2.stats.channel}", color=colors['tr2'], lw=1.0)
    axes[1].plot(tr2.times(), envelope2, color=colors['envelope2'], label='Envelope ' + tr2.stats.channel)
    axes[1].set_ylabel(f'Envelope and a ({unit})', fontsize=font)
    axes[1].set_xlabel(f"Time from {tr.stats.starttime.date} {str(tr.stats.starttime.time)[:8]} (s)", fontsize=font)
    axes[1].text(0.995, 0.2,f"Max amp (WET): {max_amp2:.2e} $\mu$m/s \n Amp ratio (FUR/WET): {amp_ratioFURWET_g}", transform=axes[1].transAxes,
                   ha='right', va='top', fontsize=font, color='black')
    axes[2].plot(ff1, spec1, color=colors['spec1'], lw=1.0, label=f"Spektrum {tr1.stats.station} {tr1.stats.channel}")
    axes[2].plot(ff2, spec2, color=colors['spec2'], lw=1.0, label=f"Spektrum {tr2.stats.station} {tr2.stats.channel}")
    axes[2].set_ylabel(f'ASD \n({unit}/Hz)', fontsize=font)
    
    axes[0].legend()
    axes[1].legend()
    axes[2].legend()


    if "fmin" in config.keys() and "fmax" in config.keys():
        axes[2].set_xlim(config['fmin'],config['fmax'])

    axes[1].set_xlabel(f"Time from {tr.stats.starttime.date} {str(tr.stats.starttime.time)[:8]} (s)",fontsize=font)     
    axes[2].set_xlabel(f"Frequency (Hz)",fontsize=font)     
    fig.suptitle(f"M{magnitude:.1f} - {distance_km:.0f} km @ {ev_depth:.0f} km | {origin_time} UTC", fontsize=font+2)

    del st_in
    return fig


#nur für translation --> nicht für rotation
def __makeplotspectraZ(st, config, fscale=None):

    from scipy import fftpack
    from andbro__fft import __fft
    import matplotlib.pyplot as plt

    st_in = st.copy()
    print([tr.stats.channel for tr in st])
    st_in = st_in.select(component="Z") 
    if len(st_in) == 0:
        print(" -> keine Z-Komponente gefunden.")
        return
    NN = len(st_in)+1

    trans_scaling = 1e6
    trans_unit = r"$\mathrm{\mu m/s}$"
   
    fig, axes = plt.subplots(NN,figsize=(18,NN*4), sharex=False)

    font = 14

    plt.subplots_adjust(hspace=0.3)

    ## _______________________________________________

    st.sort(keys=['channel'], reverse=True)

    for i, tr in enumerate(st_in):
        scaling = trans_scaling
        spec, ff, ph = __fft(tr.data*scaling, tr.stats.delta, window=None, normalize=None)
        
        if tr.stats.channel[-2] == "H":
            if i == 0:
                tr1 = tr
                spec1, ff1, ph1 = __fft(tr.data*scaling, tr.stats.delta, window=None, normalize=None)
            if i == 1:
                tr2 = tr
                spec2, ff2, ph2 = __fft(tr.data*scaling, tr.stats.delta, window=None, normalize=None)     
        
        if tr.stats.channel[1] == "H":
            sym, unit = "v", trans_unit
        else:
            unit = "Amplitude", "a.u."
    envelope1 = np.abs(hilbert(tr1.data * trans_scaling))
    envelope2 = np.abs(hilbert(tr2.data * trans_scaling))
    max_amp1 = np.max(envelope1)
    max_amp2 = np.max(envelope2)
    amp_ratioFURWET = max_amp1/max_amp2
    amp_ratioFURWET_g = np.round(amp_ratioFURWET, 5)
    # Farben definieren
    colors = {
        'tr1': '#1f77b4',        # Blau
        'envelope1': '#2ca02c',  # Grün
        'tr2': '#ff7f0e',        # Orange
        'envelope2': '#17becf',  # Türkis
        'spec1': '#1f77b4',      # Blau
        'spec2': '#ff7f0e',      # Orange
    }
    
    axes[0].plot(tr1.times(), tr1.data*trans_scaling, label=f"{tr1.stats.station} {tr1.stats.channel}", color=colors['tr1'], lw=1.0)
    axes[0].plot(tr1.times(), envelope1, color=colors['envelope1'], label='Envelope ' + tr1.stats.channel) 
    axes[0].set_ylabel(f'envelope and {sym} ({unit})', fontsize=font)
    axes[0].text(0.995, 0.15,f"Max amp (FUR): {max_amp1:.2e} $\mu$m/s", transform=axes[0].transAxes,
                   ha='right', va='top', fontsize=font, color='black')
    axes[1].plot(tr2.times(), tr2.data*trans_scaling, label=f"{tr2.stats.station} {tr2.stats.channel}", color=colors['tr2'], lw=1.0)
    axes[1].plot(tr2.times(), envelope2, color=colors['envelope2'], label='Envelope ' + tr2.stats.channel)
    axes[1].set_ylabel(f'Envelope and a ({unit})', fontsize=font)
    axes[1].set_xlabel(f"Time from {tr.stats.starttime.date} {str(tr.stats.starttime.time)[:8]} (s)", fontsize=font)
    axes[1].text(0.995, 0.2,f"Max amp (WET): {max_amp2:.2e} $\mu$m/s \n Amp ratio (FUR/WET): {amp_ratioFURWET_g}", transform=axes[1].transAxes,
                   ha='right', va='top', fontsize=font, color='black')
    axes[2].plot(ff1, spec1, color=colors['spec1'], lw=1.0, label=f"Spektrum {tr1.stats.station} {tr1.stats.channel}")
    axes[2].plot(ff2, spec2, color=colors['spec2'], lw=1.0, label=f"Spektrum {tr2.stats.station} {tr2.stats.channel}")
    axes[2].set_ylabel(f'ASD \n({unit}/Hz)', fontsize=font)
    
    axes[0].legend()
    axes[1].legend()
    axes[2].legend()


    if "fmin" in config.keys() and "fmax" in config.keys():
        axes[2].set_xlim(config['fmin'],config['fmax'])

    axes[1].set_xlabel(f"Time from {tr.stats.starttime.date} {str(tr.stats.starttime.time)[:8]} (s)",fontsize=font)     
    axes[2].set_xlabel(f"Frequency (Hz)",fontsize=font)     
    fig.suptitle(f"M{magnitude:.1f} - {distance_km:.0f} km @ {ev_depth:.0f} km | {origin_time} UTC", fontsize=font+2)

    del st_in
    return fig