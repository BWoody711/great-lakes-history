"""Export for the rebound / topobathy viewer.

  - the ICE-6G_C, ICE-7G_NA and gauge-calibrated uplift lattices at their
    native 1 degree
  - the modern topo-bathy, decimated to 0.05 deg
  - a high-resolution transect from the Chicago outlet to North Bay, with the
    modern ground, the palaeo ground and the modelled water plane at each step

The transect and the lakes follow model.py, so they carry whatever CALIB is
set to; the three lattices are always all exported, so the viewer can switch
between them without a rebuild.
"""
import numpy as np, geopandas as gpd, json, gc
import model as M
import gia_calib

STEPS = list(range(14500, 500, -500)) + [0]
TRANSECT = [(-87.80, 41.72), (-79.45, 46.32)]      # Chicago outlet -> North Bay
NPT = 320
MARKS = [(-87.80, 41.72, 'Chicago outlet'), (-86.6, 42.6, 'Lake Michigan'),
         (-85.0, 44.4, 'Michigan basin'), (-84.75, 45.75, 'Straits of Mackinac'),
         (-82.6, 45.0, 'Georgian Bay'), (-80.6, 45.9, 'French River'),
         (-79.45, 46.32, 'North Bay sill')]


def lattice(gdbpath, layerfmt, t):
    g = gpd.read_file(gdbpath, layer=layerfmt % t)
    lon = np.where(g['lon'].values > 180, g['lon'].values-360, g['lon'].values)
    lat = g['lat'].values
    m = (lon >= -99) & (lon <= -66) & (lat >= 35) & (lat <= 57)
    lo, la, gi = lon[m], lat[m], g['GIA'].values[m]
    ulo, ula = np.unique(lo), np.unique(la)
    grid = np.full((len(ula), len(ulo)), np.nan)
    grid[np.searchsorted(ula, la), np.searchsorted(ulo, lo)] = gi
    return ulo.tolist(), ula.tolist(), grid


# transect sample points
t0, t1 = TRANSECT
tlon = np.linspace(t0[0], t1[0], NPT)
tlat = np.linspace(t0[1], t1[1], NPT)
trow = ((M.N - tlat)/M.RES).astype(int)
tcol = ((tlon - M.W)/M.RES).astype(int)
dist = np.cumsum(np.r_[0, np.hypot(np.diff(tlon)*111.32*np.cos(np.radians(tlat[:-1])),
                                   np.diff(tlat)*111.32)])

G6, G7 = gia_calib.GDB['ICE6G'], gia_calib.GDB['ICE7G']

# The calibration is a rate field in space times a scalar in time, so the
# calibrated lattice is built from the same fit the model uses rather than
# re-derived here. Evaluated on the 1 deg lattice it is only the correction
# the viewer needs to show; the taper keeps it inside the gauge network.
CF = gia_calib.fit('ICE7G')

lat6 = {}; lat7 = {}; latc = {}
for t in STEPS:
    if t == 0:
        continue
    lo, la, g6 = lattice(*G6, t); _, _, g7 = lattice(*G7, t)
    LO, LA = np.meshgrid(np.array(lo), np.array(la))
    dv = gia_calib.surface(CF, LO, LA)/1e4                  # m/yr
    lat6[t] = np.round(g6, 1).tolist(); lat7[t] = np.round(g7, 1).tolist()
    latc[t] = np.round(g7 + dv*gia_calib.time_factor(t), 1).tolist()
    print('lattice', t, flush=True)
LON1D, LAT1D = lo, la

# decimated modern topo-bathy for the map background
dec = 6                                   # 30" * 6 = 3 arc-minutes
dm = M.dem[:M.dem.shape[0]//dec*dec, :M.dem.shape[1]//dec*dec]
dm = dm.reshape(dm.shape[0]//dec, dec, dm.shape[1]//dec, dec).mean(axis=(1, 3))
print('background grid', dm.shape)

frames = []
for t in STEPS:
    label = None if t == 0 else f'{t/1000:g}ka'
    r = M.run(t, label)
    surf, filled, ice = r['surf'], r['filled'], r['ice']
    water = (filled - np.where(ice, 9000.0, surf)) > 0.05
    wl = np.where(water[trow, tcol], filled[trow, tcol], np.nan)
    frames.append(dict(
        t=t, esl=round(r['esl'], 1),
        ground=np.round(surf[trow, tcol], 1).tolist(),
        water=[None if np.isnan(v) else round(float(v), 1) for v in wl],
        ice=[bool(v) for v in ice[trow, tcol]]))
    del r, surf, filled, ice, water
    gc.collect()
    print('transect', t, flush=True)

json.dump(dict(
    lon1d=LON1D, lat1d=LAT1D, ice6=lat6, ice7=lat7, cal=latc,
    calib=dict(on=bool(M.CALIB), tau=gia_calib.TAU, deg=gia_calib.DEGREE,
               lo=round(CF['lo'], 1), hi=round(CF['hi'], 1),
               wrms_before=round(CF['wrms_before'], 2),
               wrms_after=round(CF['wrms_after'], 2),
               n=len(CF['G']['name'])),
    bg=dict(w=dm.shape[1], h=dm.shape[0], bbox=[M.W, M.S, M.E, M.N],
            z=np.round(dm).astype('int16').ravel().tolist()),
    transect=dict(lon=np.round(tlon, 3).tolist(), lat=np.round(tlat, 3).tolist(),
                  dist=np.round(dist, 1).tolist(),
                  modern=np.round(M.dem[trow, tcol], 1).tolist(),
                  marks=[dict(lon=a, lat=b, name=c) for a, b, c in MARKS]),
    frames=frames), open('out/rebound.json', 'w'), separators=(',', ':'))
import os
print('wrote out/rebound.json  %.2f MB' % (os.path.getsize('out/rebound.json')/1e6))
