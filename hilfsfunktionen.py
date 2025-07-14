from obspy.signal.rotate import rotate_ne_rt
import obspy as obs
import numpy as np
def transform_to_rtz(stream, azimuth):
    """
    Transform stream from NEZ to RTZ coordinate system using azimuth.
    """
    # Extrahiere horizontale Komponenten
    n_comp = stream.select(component="N")[0]
    e_comp = stream.select(component="E")[0]
    z_comp = stream.select(component="Z")[0]  # bleibt gleich für RTZ

    # Drehe N/E → R/T
    r_data, t_data = rotate_ne_rt(n_comp.data, e_comp.data, azimuth)

    # Erzeuge neue Trace-Objekte für R/T/Z
    r_trace = n_comp.copy()
    r_trace.data = r_data
    r_trace.stats.channel = r_trace.stats.channel[:-1] + "R"

    t_trace = e_comp.copy()
    t_trace.data = t_data
    t_trace.stats.channel = t_trace.stats.channel[:-1] + "T"

    # Z bleibt unverändert
    z_trace = z_comp.copy()
    z_trace.stats.channel = z_trace.stats.channel[:-1] + "Z"

    return obs.Stream(traces=[r_trace, t_trace, z_trace])

def safe_divide(spec1_s, spec2_s):
    try:
        # Prüfe ob beide Arrays kompatibel sind
        if not isinstance(spec1_s, np.ndarray) or not isinstance(spec2_s, np.ndarray):
            spec1_s = np.array(spec1_s)
            spec2_s = np.array(spec2_s)

        # Prüfe auf identische Shapes
        if spec1_s.shape != spec2_s.shape:
            print(f"Shape mismatch: {spec1_s.shape} vs {spec2_s.shape}")
            return None

        # Division mit Fehlerunterdrückung
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.divide(spec1_s, spec2_s)
            ratio[~np.isfinite(ratio)] = 0  # optional: NaN und inf zu 0 machen

        return ratio

    except Exception as e:
        print(f"Fehler bei safe_divide: {e}")
        return None