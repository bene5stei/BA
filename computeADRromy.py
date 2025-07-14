import timeit
import os
import sys
import obspy as obs
import matplotlib.pyplot as plt
import numpy as np
from obspy.clients.fdsn import Client
import timeit
import matplotlib.colors
from obspy import UTCDateTime, Stream, read_inventory
from obspy.clients import fdsn
from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics import locations2degrees
from obspy.clients.fdsn import Client, RoutingClient
from obspy.signal import array_analysis as AA
from obspy.signal.util import util_geo_km
from obspy.signal.rotate import rotate2zne
from datetime import datetime
import warnings
def __compute_romy_adr(tbeg, tend, submask='all', ref_station='GR.FUR', excluded_stations=[], map_plot=False, verbose=True):

    """
    rotation_X = -u_nz
    rotation_Y =  u_ez
    rotation_Z = 0.5*(u_ne-u_en)
    """

    # start timer for runtime
    start_timer = timeit.default_timer()

    # _____________________________________________________

    # generate configuration object
    config = {}

    # change to UTCDateTime object
    config['tbeg'] = tbeg
    config['tend'] = tend
    # select the fdsn client for the stations
    # config['fdsn_client'] = {"BW":Client('http://jane'), "GR":Client('BGR')}
    config['fdsn_client'] = {"BW":Client('https://jane.geophysik.uni-muenchen.de'), "GR":Client('BGR')}

    # define path to local inventory data
    config['data_to_inventory'] = None #"./data/stationxml_ringlaser/"

    # define output seed
    config['out_seed'] = "BW.ROMY"

    # add submask
    config['submask'] = submask

    # specify location codes for subarrays
    if config['submask'] == "inner":
        config['location'] = "21"
    elif config['submask'] == "outer":
        config['location'] = "22"
    elif config['submask'] == "all":
        config['location'] = "23"

    # select stations to consider:
    if submask is not None:

        if submask == "inner":
            config['subarray_mask'] = [0, 1, 2, 3]
            config['freq1'] = 0.4    #urspr. 0.2
            config['freq2'] = 3.7     #urspr. 1.0

        elif submask == "outer":
            config['subarray_mask'] = [0, 4, 5, 6, 7, 8]
            config['freq1'] = 0.04   #urspr. 0.01
            config['freq2'] = 0.2    #urspr. 0.3

        elif submask == "all":
            config['subarray_mask'] = [0, 1, 2, 3, 4, 5, 6, 7, 8]
            config['freq1'] = 0.04   #urspr. 0.01
            config['freq2'] = 0.2    #urspr. 0.1

    # decide if information is printed while running the code
    config['verbose'] = verbose

    # specify reference station (where to compute adr)
    config['reference_station'] = ref_station

    # ROMY array station information
    config['array_stations'] = ['GR.FUR', 'BW.FFB1', 'BW.FFB2', 'BW.FFB3',
                                'BW.BIB', 'BW.TON', 'BW.GELB', 'BW.ALFT', 'BW.GRMB',
                               ]

    # select stations of subarray
    config['subarray_stations'] = []
    for i in config['subarray_mask']:
        if config['array_stations'][i] not in excluded_stations:
            config['subarray_stations'].append(config['array_stations'][i])

    config['subarray_sta'] = config['subarray_stations']

    # specify if bandpass is applied to data
    # config['prefilt'] = (0.001, 0.01, 5, 10)
    config['apply_bandpass'] = True

    # adr parameters
    config['vp'] = 3400#1 # 5000 #6264. #1700
    config['vs'] = 3400#1 # 3500 #3751. #1000
    config['sigmau'] = 1e-7 # 0.0001


# ab hier test: der noch fehlende code

    # launch a times
    start_timer1 = timeit.default_timer()

    # status of stations loaded
    config['stations_loaded'] = np.ones(len(config['subarray_stations']))

    # request data for pfo array
    st, ref_station, config = __get_data(config)

    # check for reference station
    if config['reference_station'] not in config['subarray_stations']:
        print(f" -> reference station not in station set!!!")

    # processing
    st = st.detrend("linear")
    st = st.detrend("demean")

    print(st)

    # bandpass filter
    if config['apply_bandpass']:
        st = st.taper(0.02, type="cosine")
        st = st.filter('bandpass', freqmin=config['freq1'], freqmax=config['freq2'], corners=4, zerophase=True)
        print(f" -> bandpass: {config['freq1']} - {config['freq2']} Hz")

    ## plot station coordinates for check up
    if map_plot:
        stas = []
        for tr in st:
            if tr.stats.station not in stas:
                stas.append(tr.stats.station)

        # make checkup plot
        import matplotlib.pyplot as plt
        plt.figure()
        for c, s in zip(config['coo'], stas):
            if config['verbose']:
                print(s, c)
            plt.scatter(c[0], c[1], label=s)
            plt.legend()
        plt.show();

    # check if enough stations for ADR are available otherwise continue
    if len(st) < 9:
        print(" -> not enough stations (< 3) for ADR computation!")
        return
    else:
        if config['verbose']:
            print(f" -> continue computing ADR for {int(len(st)/3)} of {len(config['subarray_mask'])} stations ...")

    # homogenize the time line
    st = __adjust_time_line(st, reference=config['reference_station'])

    # check for same amount of samples
    __check_samples_in_stream(st, config)

    print(st)

    # compute array derived rotation (ADR)
    try:
        rot = __compute_ADR(st, config, ref_station)
    except Exception as e:
        print(e)
        return None, None

    # trim to requested interval
    rot = rot.trim(config['tbeg'], config['tend'], nearest_sample=False)

    # stop times
    stop_timer1 = timeit.default_timer()
    print(f"\n -> Runtime: {round((stop_timer1 - start_timer1)/60, 2)} minutes\n")

    return rot
def __get_data(config):

        config['subarray'] = []

        st = Stream()

        coo = []
        print(config['array_stations'])  # sollte nur Strings enthalten
        for k, station in enumerate(config['subarray_stations']):
            print(" -> station entry:", station, type(station))  # Debug-Ausgabe
            net, sta = station.split(".")
            print(net, sta)
            loc, cha = "", "BH*"

            # modify time interval
            t1, t2 = config['tbeg']-3600, config['tend']+3600

            # querry inventory data
            try:
                for stage in ['local', 'jane']:
                    #if stage == 'local' and config['data_to_inventory'] != None:
                       # try:
                            # load local version
                          #  inventory = read_inventory(config['data_to_inventory']+f"stationxml_ringlaser/station_{net}.{sta}.xml")
                          #  if config['verbose']:
                           #     print(f" -> loaded local inventory for {net}.{sta}")
                          #  break
                      #  except:
                            #pass
                    if stage == 'jane':
                        try:
                            inventory = config['fdsn_client'][net].get_stations(
                                                                            network=net,
                                                                            station=sta,
                                                                            location=loc,
                                                                            channel=cha,
                                                                            starttime=t1,
                                                                            endtime=t2,
                                                                            level="response"
                                                                                )
                            if config['verbose']:
                                print(f" -> loaded jane inventory for {net}.{sta}")
                            break
                        except:
                            pass
            except Exception as e:
                print(f" -> failed to load inventory for {net}.{sta} ...!", e)
                inventory = None

            if inventory is not None:

                # add coordinates
                l_lon = float(inventory.get_coordinates('%s.%s.%s.%sZ'%(net,sta,loc,cha[:2]))['longitude'])
                l_lat = float(inventory.get_coordinates('%s.%s.%s.%sZ'%(net,sta,loc,cha[:2]))['latitude'])
                height = float(inventory.get_coordinates('%s.%s.%s.%sZ'%(net,sta,loc,cha[:2]))['elevation'])

                # store reference coordinates
                print(config['reference_station'])
                if sta == config['reference_station'].split(".")[1]:
                    o_lon, o_lat, o_height = l_lon, l_lat, height

                # compute relative distances to reference station in km
                lon, lat = util_geo_km(o_lon, o_lat, l_lon, l_lat)

            # try to get waveform data
            try:

                # load waveforms for station
                stats = config['fdsn_client'][net].get_waveforms(
                                                                network=net,
                                                                station=sta,
                                                                location=loc,
                                                                channel=cha,
                                                                starttime=config['tbeg']-3600,
                                                                endtime=config['tend']+3600,
                                                                attach_response=True,
                                                                )
            except Exception as e:
                print(e)
                print(f" -> failed to load waveforms for {net}.{sta} ...")
                config['stations_loaded'][k] = 0
                continue

            # merge if masked
            if len(stats) > 3:
                print(f" -> merging stream. Length: {len(stats)} -> 3") if config['verbose'] else None
                stats.merge(method=1, fill_value="interpolate")

            # successfully obtained
            if len(stats) == 3:
                if config['verbose']:
                    print(f" -> obtained waveforms for {net}.{sta}\n")

            # remove response [VEL -> rad/s | DISP -> rad]
            stats = stats.remove_sensitivity(inventory)
            # stats = stats.remove_response(inventory, output="VEL", water_level=None, pre_filt=[0.005, 0.008, 8, 10])
            # stats = stats.remove_response(inventory, output="VEL", water_level=1)

            # rotate to ZNE (if FFB* stations are included)
            try:
                if config['submask'] in ["inner", "all"]:
                    stats = stats.rotate(method='->ZNE', inventory=inventory, components=['Z12'])
                else:
                    stats = stats.rotate(method='->ZNE', inventory=inventory, components=['ZNE'])
            except:
                print(f" -> failed to rotate to ZNE for {net}.{sta}")
                continue

            # resampling using decitmate
            stats = stats.detrend("linear");
            stats = stats.detrend("simple");

            # stats = stats.taper(0.01);
            stats = stats.filter("highpass", freq=0.001, corners=4, zerophase=True);
            stats = stats.filter("lowpass", freq=5.0, corners=4, zerophase=True); # added 2025-02-18

            # store stream of reference station as template
            if station == config['reference_station']:
                ref_station = stats.copy()

                if config['submask'] == "inner":
                    # upsample (FUR) from 20 to 40 Hz
                    stats = stats.resample(40.0, no_filter=True)
                    stats = stats.trim(t1, t2)

            # resample to output sampling rate
            if config['submask'] == "inner":
                stats = stats.decimate(2, no_filter=False)

            elif config['submask'] == "all":
                for tr in stats:
                    if "FFB" in tr.stats.station:
                        tr = tr.decimate(2, no_filter=False)

            # store distances as meters
            coo.append([lon*1000, lat*1000, height-o_height])

            # add station data to stream
            st += stats

            # update subarray
            config['subarray'].append(f"{stats[0].stats.network}.{stats[0].stats.station}")

        # trim to interval
        # stats.trim(config['tbeg'], config['tend'], nearest_sample=False)

        # convert to array object
        config['coo'] = np.array(coo)

        config['subarray_stations'] = config['subarray']

        if config['verbose']:
            print(f" -> obtained: {int(len(st)/3)} of {len(config['subarray_stations'])} stations!")

        # check for reference station
        if config['reference_station'] not in config['subarray_stations']:
            print(f" -> reference station not in station set!!!")
            return st, Stream(), config

        if len(st) == 0:
            return st, Stream(), config
        else:
            return st, ref_station, config

def __compute_ADR(st, config, ref_station):

        # prepare data arrays
        tsz, tsn, tse = [], [], []
        for tr in st:
            try:
                if "Z" in tr.stats.channel:
                    tsz.append(tr.data)
                elif "N" in tr.stats.channel:
                    tsn.append(tr.data)
                elif "E" in tr.stats.channel:
                    tse.append(tr.data)
            except:
                print(" -> stream data could not be appended!")

        # make sure input is array type
        tse, tsn, tsz = np.array(tse), np.array(tsn), np.array(tsz)

        # define array for subarray stations with linear numbering
        substations = np.arange(len(config['subarray_stations']))

        try:
            result = AA.array_rotation_strain(substations,
                                              np.transpose(tse),
                                              np.transpose(tsn),
                                              np.transpose(tsz),
                                              config['vp'],
                                              config['vs'],
                                              config['coo'],
                                              config['sigmau'],
                                             )

        except Exception as e:
            print(e)
            print("\n -> failed to compute ADR ...")
            return None

        # create rotation stream and add data
        rotsa = ref_station.copy()

        rotsa[0].data = result['ts_w3']
        rotsa[1].data = result['ts_w2']
        rotsa[2].data = result['ts_w1']

        rotsa[0].stats.channel='BJZ'
        rotsa[1].stats.channel='BJN'
        rotsa[2].stats.channel='BJE'

        rotsa[0].stats.station=config['out_seed'].split(".")[1]
        rotsa[1].stats.station=config['out_seed'].split(".")[1]
        rotsa[2].stats.station=config['out_seed'].split(".")[1]

        rotsa[0].stats.network=config['out_seed'].split(".")[0]
        rotsa[1].stats.network=config['out_seed'].split(".")[0]
        rotsa[2].stats.network=config['out_seed'].split(".")[0]

        rotsa[0].stats.location=config['location']
        rotsa[1].stats.location=config['location']
        rotsa[2].stats.location=config['location']

        rotsa = rotsa.detrend('linear')

        return rotsa


def __adjust_time_line(st0, reference="GR.FUR"):

        Rnet, Rsta = reference.split(".")

        ref_start = st0.select(network=Rnet, station=Rsta)[0].stats.starttime
        ref_times = st0.select(network=Rnet, station=Rsta)[0].times()

        dt = st0.select(network=Rnet, station=Rsta)[0].stats.delta

        for tr in st0:
            times = tr.times(reftime=ref_start)

            tr.data = np.interp(ref_times, times, tr.data)
            tr.stats.starttime = ref_start

        return st0
def __check_samples_in_stream(st, config):

        Rnet, Rsta = config['reference_station'].split(".")

        Rsamples = st.select(network=Rnet, station=Rsta)[0].stats.npts

        for tr in st:
            if tr.stats.npts != Rsamples:
                print(f" -> removing {tr.stats.station} due to improper number of samples ({tr.stats.npts} not {Rsamples})")
                st.remove(tr)

        return st

    