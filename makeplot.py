import matplotlib.pyplot as plt
def __makeplot(config, st):

    import matplotlib.pyplot as plt
    st_in = st.copy()

    #fig, ax = plt.subplots(11,1, figsize=(15,20), sharex=True)    #eig fig, ax = plt.subplots(6,1, figsize=(15,10), sharex=True)
    fig, ax = plt.subplots(len(st), 1, figsize=(15, 2*len(st)), sharex=True)  #test

    font = 14

    time_scaling, time_unit = 1, "sec"
    #Möglichkeit für bestimmte Einheit
    #rot_scaling = 1    #urspr. 1e9
    #trans_scaling = 1    #urspr. 1e6

    for i, tr in enumerate(st_in):

        if i in [0,1,2]:
            #ax[i].set_ylabel(r"$\omega$ (nrad/s)", fontsize=font)
            ax[i].plot(tr.times()/time_scaling, tr.data, 'k', label=tr.stats.station+"."+tr.stats.channel)

        elif i in [3,4,5,6,7,8,9,10,11]:
            #ax[i].set_ylabel(r"u ($\mu$m/s)", fontsize=font)
            ax[i].plot(tr.times()/time_scaling, tr.data, 'k', label=tr.stats.station+"."+tr.stats.channel)

        ax[i].legend(loc=1)

    ax[10].set_xlabel(f"Time ({time_unit}) from {st[0].stats.starttime.date} {str(st[0].stats.starttime.time).split('.')[0]} UTC", fontsize=font)
    ax[0].set_title(config['title']+f" | {config['fmin']} - {config['fmax']} Hz", fontsize=font, pad=10)

    plt.show();
    del st_in
    return fig