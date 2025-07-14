#!/usr/bin/python
#
# querry seismic traces and station data
#
# by AndBro @2021
# __________________________



def __querrySeismoData(seed_id=None, starttime=None, endtime=None, repository=None, path=None, restitute=True, detail=None, fill_value=None):

    '''

    Querry stream and station data of OBS

    VARIABLES:
        seed_id:    code of seismic stations (e.g. "BW.ROMY..BJU")
        tbeg:       begin of time period
        tend:       temporal length of period
        repository:      location to retrieve data from: 'local', 'online', 'george', 'archive'
        path:       if repository is 'local', path to the data has to be provided. 
                    file names are assumed to be: BW.ROMY..BJV.D.2021.059
        resitute:   if response is removed or not
        detail:     if information is printed at the end


    DEPENDENCIES:
        import sys

        from obspy.clients.fdsn import Client, RoutingClient
        from obspy.core.util import AttribDict
        from obspy import UTCDateTime, Stream, Inventory, read, read_inventory
        from numpy import ma
        from os.path import isfile

    OUTPUT:
        out1: stream
        out2: inventory

    EXAMPLE:

        >>> st, inv = __querrySeismoData(seed_id="BW.DROMY..FJZ",
                             starttime="2021-02-18 12:00",
                             endtime="2021-02-18 12:10",
                             repository='local',
                             path='/home/andbro/Documents/ROMY/data/',
                             restitute=True,
                             detail=True,
                            )

    '''

    ## importing libraries and modules
    import sys

    from obspy.clients.fdsn import Client, RoutingClient
    from obspy.core.util import AttribDict
    from obspy import UTCDateTime, Stream, Inventory, read, read_inventory
    from numpy import ma
    from os.path import isfile, isdir
    from obspy.clients.filesystem.sds import Client as sdsclient


    ## split seed_id string
    net, sta, loc, cha = seed_id.split(".")

    ## convert to datetime format in case provided as string
    ## add +-1 second to allow correct trimming at the end
    doy = UTCDateTime(starttime).julday
    year = UTCDateTime(starttime).year

    tbeg = UTCDateTime(starttime)-1
    tend = UTCDateTime(endtime)+1


    ## check path if provided
    if path:
        if path[-1] == "/":
            path = path[:-1]

    details = []

    st = Stream()
    inv = Inventory()

    ## check if input variables are as expected
    for arg in [net, sta, loc, cha, tbeg, tend]:
        if arg is None and not 'loc':
            raise NameError(print(f"\nwell, {arg} has not been defined after all!"))
            sys.exit()



    ## __________________________________________________________________
    ##

    if repository == 'online':
        ## attempting to get data from either EIDA or IRIS.
        try:
            route = RoutingClient("eida-routing")
            details.append(f"RoutingClient: {route}")

            if route:
                inv = route.get_stations(network=net, station=sta, location=loc, channel=cha,
                                         starttime=tbeg, endtime=tend, level="response")


                st = route.get_waveforms(network=net, station=sta, location=loc, channel=cha,
                                         starttime=tbeg, endtime=tend)

                try:
                    st[0].stats.coordinates = AttribDict({  'latitude':  inv[0][0].latitude,
                                                            'elevation': inv[0][0].elevation,
                                                            'longitude': inv[0][0].longitude,
                                                          })
                except:
                    details.append(f"no coordinates added to {sta[0]}")

        except:
            route = RoutingClient("iris-federator")
            details.append(f"RoutingClient: iris-federator")

            if route:
                inv = route.get_stations(network=net, station=sta, location=loc, channel=cha,
                                         starttime=tbeg, endtime=tend, level="response")

                st = route.get_waveforms(network=net, station=sta, location=loc, channel=cha,
                                         starttime=tbeg, endtime=tend)



    ## __________________________________________________________________
    ## load local data

    if repository == 'local':

        year = tbeg.year
        doy  = tbeg.julday

        ## adjust doy string to 3 chars
        doy = str(doy).rjust(3, "0")

        if not path:
            print("no path provided!")
            sys.exit()
        try:
            st = read(path,
                      starttime=tbeg,
                      endtime=tend,
                      )
            st.merge()
        except:
            print("failed to load mseed")

        try:
            try:
                route = RoutingClient("eida-routing")
                details.append(f"RoutingClient: {route}")

                if route:
                    inv = route.get_stations(network=net,
                                             station=sta,
                                             location=loc,
                                             channel=cha,
                                             starttime=tbeg,
                                             endtime=tend,
                                             level="response",
                                            )
            except:
                inv = read_inventory(path+f"/{sta}.xml")

        except:
            details.append("failed to obtain inventory")

    ## __________________________________________________________________
    ## load data from george

    if repository.lower() in ['george', 'jane']:

        if repository.lower() == 'george':
            waveform_client = Client(base_url='http://george', timeout=200)
        elif repository.lower() == 'jane':
            waveform_client = Client(base_url='https://jane.geophysik.uni-muenchen.de', timeout=200)

        ## WAVEFORMS
        try:
            st = waveform_client.get_waveforms(location=loc,
                                               channel=cha,
                                               network=net,
                                               station=sta,
                                               starttime=tbeg,
                                               endtime=tend,
                                               level='response',
                                              );
        except:
            st = waveform_client.get_waveforms(location=loc,
                                               channel=cha,
                                               network=net,
                                               station=sta,
                                               starttime=tbeg,
                                               endtime=tend,
                                          	);

        ## INVENTORY
        try:
            try:
                route = RoutingClient("eida-routing")

                inv = route.get_stations(network=net,
                                         station=sta,
                                         location=loc,
                                         channel=cha,
                                         starttime=tbeg,
                                         endtime=tend,
                                         level="response",
                                        );
            except:

                route = Client("LMU")

                inv = route.get_stations(network=net,
                                         station=sta,
                                         location=loc,
                                         channel=cha,
                                         starttime=tbeg,
                                         endtime=tend,
                                         level="response",
                                        );

        except:

            try:
                if sta == "ROMY":
                    inv = read_inventory("/home/andbro/Documents/ROMY/data/ROMY.xml")
            except:
                details.append("no inventory found")


    ## __________________________________________________________________
    ## load data from archive

    if repository == 'archive':

        if path is not None:
            path2sds = path
        else:
            path2sds = f"/import/freenas-ffb-01-data/romy_archive/"

        if not isdir(path2sds):
            print(f"no such path: \n {path2sds}")
            sys.exit()

        ## define SDS client
        sds_client = sdsclient(path2sds, sds_type='D', format='MSEED')

        ## stream data from archive
        try:
            st = sds_client.get_waveforms(net, sta, loc, cha, tbeg, tend, merge=-1)
        except:
            print(f"failed to get data from archive: \n {path2sds} \n for seed: {seed_id}")
            # sys.exit()

        ## check if stream is masked / fragmented and merge as consequence
        if len(st) > 1:
            frags = len(st)
            st.merge()
            details.append(f"merged fragmented stream (fragments = {frags})")

        ## inventory
        try:
            try:
                route = RoutingClient("eida-routing")

                inv = route.get_stations(network=net,
                                         station=sta,
                                         location=loc,
                                         channel=cha,
                                         starttime=tbeg,
                                         endtime=tend,
                                         level="response",
                                        );
            except:

                route = Client("LMU")

                inv = route.get_stations(network=net,
                                         station=sta,
                                         location=loc,
                                         channel=cha,
                                         starttime=tbeg,
                                         endtime=tend,
                                         level="response",
                                        );

        except:
            try:
                if sta == "ROMY":
                    inv = read_inventory("/home/andbro/Documents/ROMY/data/ROMY.xml");
            except:
                details.append("failed to find an inventory")


    ## __________________________________________________________________
    ## add coordinate information to each channel, if available in inventory

    try:
        coords = inv.get_coordinates(f'{net}.{sta}.{loc}.{cha}', tbeg)
        for tr in st:
            tr.stats.coordinates = AttribDict(coords)
    except:
        details.append(" -> Coordinates could not be attached to Inventory!")


    ## __________________________________________________________________
    ## remove response of instrument specified in inventory


    if restitute:

        if cha[1] == "H":
            pre_filter = [0.001, 0.005, 95, 100]
            water_level = 50

            out="VEL"  # alternatives: "DISP" "ACC"
            try:      
                st.remove_response(
                                    inventory=inv, 
                                    pre_filt=pre_filter,
                                    output=out,
                                  #  water_level=water_level,   
                                    zero_mean=True,                             
                                    )

                details.append(f"OUT: {out}")
                details.append(f'pre-filter: {pre_filter}')
                    
            except:
                details.append("no response removed")

        elif cha[1] == "J":

            pre_filter = [0.001, 0.005, 95, 100]
            water_level = 50

            out="VEL"  # alternatives: "DISP" "ACC"
            try:      
                st.remove_sensitivity(inventory=inv)
                    
            except:
                details.append("no sensitivity removed")

    ## __________________________________________________________________
    ## 

#    from numpy import nan
#    if fill_value is None:
#        fill_value=nan

#    st.merge();
#    for tr in st:
#        if ma.is_masked(tr.data):
##             tr.data = ma.filled(tr.data, fill_value=-999999)
#            tr.data = ma.filled(tr.data, fill_value=fill_value)
#            details.append(f"trace: {tr.stats.network}.{tr.stats.station}.{tr.stats.location}.{tr.stats.channel} is masked")
                 

    ## __________________________________________________________________
    ##     

    ## print processing details if set true
    if detail is True:
        for det in details:
            print(det)
    
    
    ## final trim to exact time window of request  -> fails with masked streams! 
    st.trim(UTCDateTime(starttime), UTCDateTime(endtime));
    
    return st, inv

## End of File
