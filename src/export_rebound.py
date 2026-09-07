"""Export for the rebound / topobathy viewer.

  - the ICE-6G_C and ICE-7G_NA uplift lattices at their native 1 degree
  - the modern topo-bathy, decimated to 0.05 deg
  - a high-resolution transect from the Chicago outlet to North Bay, with the
    modern ground, the palaeo ground and the modelled water plane at each step
"""
import numpy as np, geopandas as gpd, json, gc
import model as M

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

G6 = ('data/ice6g_pts/ICE6G_Data_points.gdb', 'ICE6G_datapoint_%05d')
G7 = ('data/ice7g/ICE7G_Data_points_all.gdb', 'ICE7G_datapoints_%05d')

lat6 = {}; lat7 = {}
for t in STEPS:
    if t == 0:
        continue
    lo, la, g6 = lattice(*G6, t); _, _, g7 = lattice(*G7, t)
    lat6[t] = np.round(g6, 1).tolist(); lat7[t] = np.round(g7, 1).tolist()
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
    lon1d=LON1D, lat1d=LAT1D, ice6=lat6, ice7=lat7,
    bg=dict(w=dm.shape[1], h=dm.shape[0], bbox=[M.W, M.S, M.E, M.N],
            z=np.round(dm).astype('int16').ravel().tolist()),
    transect=dict(lon=np.round(tlon, 3).tolist(), lat=np.round(tlat, 3).tolist(),
                  dist=np.round(dist, 1).tolist(),
                  modern=np.round(M.dem[trow, tcol], 1).tolist(),
                  marks=[dict(lon=a, lat=b, name=c) for a, b, c in MARKS]),
    frames=frames), open('out/rebound.json', 'w'), separators=(',', ':'))
import os
print('wrote out/rebound.json  %.2f MB' % (os.path.getsize('out/rebound.json')/1e6))
