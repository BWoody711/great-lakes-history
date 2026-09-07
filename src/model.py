"""
Paleo-hydrography of the Great Lakes basins.

DEM        : ETOPO 2022 15" (NOAA NCEI), includes the NOAA Great Lakes bathymetry
GIA        : ICE-6G_C (VM5a) land deformation, Godbout, Brouard & Roy, PANGAEA 947536
Ice margins: NADI-1 optimal isochrones, Dalton et al. 2023, QSR 321, 108345

Per 500-yr timestep:
  paleo surface = modern topo-bathy - uplift since t
  ice = wall; domain edge = open drain
  grayscale reconstruction by erosion -> filled surface
  water = filled - surface; lakes labelled; level = filled value; spill = lowest
  non-lake neighbour, i.e. the sill the lake actually overflows at.
"""
import numpy as np, geopandas as gpd, json, os, sys, time
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import gaussian_filter1d
from scipy import ndimage as ndi
from skimage.morphology import reconstruction
from rasterio.features import rasterize, shapes
from rasterio.transform import from_origin
from shapely.geometry import shape
from shapely.ops import unary_union

W, E, S, N = -95.0, -70.0, 38.0, 53.0
RES = 1/120.0
MODEL = 'ICE7G'
import os as _os
MARGIN = _os.environ.get('MARGIN', 'OPTIMAL')   # 'OPTIMAL', 'MAX' or 'MIN'      # 'ICE6G' or 'ICE7G'
GDB = {'ICE6G': ('data/raw/ice6g_pts/ICE6G_Data_points.gdb', 'ICE6G_datapoint_%05d'),
       'ICE7G': ('data/raw/ice7g/ICE7G_Data_points_all.gdb', 'ICE7G_datapoints_%05d')}[MODEL]
NADI = 'data/raw/nadi1/NADI-1 shapefiles Dalton et al. QSR'

dem = np.load('data/dem30_min.npy')
ny, nx = dem.shape
TR = from_origin(W, N, RES, RES)
lon_c = W + (np.arange(nx)+0.5)*RES
lat_c = N - (np.arange(ny)+0.5)*RES
AREA_ROW = ((RES*111.32)*(RES*111.32*np.cos(np.radians(lat_c)))).astype('f4')

MIN_AREA_KM2 = 120.0
SIMPLIFY = 0.010

# Every 500-yr solved timestep, oldest to 0. Matches STEPS in export.py /
# export_rebound.py; kept here too because the GIA cube below is built once
# for the whole run, not per call.
STEPS_T = list(range(14500, 500, -500)) + [0]


def _gia_grid_raw(t_yr):
    g = gpd.read_file(GDB[0], layer=GDB[1] % t_yr)
    lon = np.where(g['lon'].values > 180, g['lon'].values-360, g['lon'].values)
    lat = g['lat'].values
    deep = g['Topo_0Ka'].values < -4000
    esl = -float(np.median(g['Topo_diff'].values[deep]))
    m = (lon >= W-5) & (lon <= E+5) & (lat >= S-5) & (lat <= N+5)
    lo, la, gi = lon[m], lat[m], g['GIA'].values[m].astype('f8')
    ulo, ula = np.unique(lo), np.unique(la)
    grid = np.full((len(ula), len(ulo)), np.nan)
    grid[np.searchsorted(ula, la), np.searchsorted(ulo, lo)] = gi
    assert not np.isnan(grid).any(), 'gappy GIA lattice'
    spl = RectBivariateSpline(ula, ulo, grid, kx=3, ky=3, s=0)
    # lat_c runs north->south; evaluate ascending then flip back
    z = spl(lat_c[::-1], lon_c)[::-1]
    return z.astype('f4'), esl


_gia_cube_cache = None


def _gia_cube():
    """Uplift-since-t for every solved timestep, smoothed along the time axis.

    ICE-6G/ICE-7G give GIA on a 1 deg lattice (see README); RectBivariateSpline
    interpolates that down to the 30" grid, but at a small, shallow basin (a
    few tens of km across, a few tens of m of sill-to-floor relief -- Lake
    Simcoe is the case that surfaced this) a handful of metres of spline
    wobble between control points is enough to open or close its sill. Real
    postglacial uplift is a smooth, slowly-relaxing function of time at any
    fixed point; it does not vanish for a few timesteps and come back. Fitting
    each timestep's spline independently, from a different scattered point set
    each time, has no reason to respect that -- so any single-timestep jump
    that isn't echoed by its neighbours is spline noise, not signal, and is
    filtered out here rather than left to flip small basins on and off.
    This does not add resolution the 1 deg lattice doesn't have; it only
    removes noise the per-timestep fit adds on top of it.
    """
    global _gia_cube_cache
    if _gia_cube_cache is not None:
        return _gia_cube_cache
    cache_path = f'data/gia_cube_{MODEL}.npz'
    if os.path.exists(cache_path):
        d = np.load(cache_path)
        _gia_cube_cache = (d['gia'], d['esl'])
        return _gia_cube_cache
    print(f'building smoothed GIA cube for {MODEL}, {len(STEPS_T)} timesteps...', flush=True)
    grids, esls = [], []
    for t in STEPS_T:
        if t == 0:
            grids.append(np.zeros_like(dem)); esls.append(0.0)
        else:
            g, e = _gia_grid_raw(t)
            grids.append(g); esls.append(e)
            print(f'  {t}', flush=True)
    gia = np.stack(grids).astype('f4')
    esl = np.array(esls, dtype='f8')
    # sigma=1 step (~500-1000 yr, since the final gap is 1000 yr): a light
    # touch that damps single-timestep spikes without smearing the multi-ka
    # trend (deglaciation, forebulge collapse) the model depends on elsewhere.
    gia = gaussian_filter1d(gia, sigma=1.0, axis=0, mode='nearest')
    esl = gaussian_filter1d(esl, sigma=1.0, axis=0, mode='nearest')
    np.savez(cache_path, gia=gia, esl=esl)
    _gia_cube_cache = (gia, esl)
    return _gia_cube_cache


def gia_grid(t_yr):
    gia, esl = _gia_cube()
    idx = STEPS_T.index(t_yr)
    return gia[idx], float(esl[idx])


def ice_mask(label):
    f = f'{NADI}/{label}_cal_{MARGIN}_NADI-1_Dalton_etal_QSR.shp'
    g = gpd.read_file(f, bbox=(W-1, S-1, E+1, N+1))
    if len(g) == 0:
        return np.zeros((ny, nx), bool), None
    geom = unary_union(list(g.geometry.values))
    r = rasterize([(geom, 1)], out_shape=(ny, nx), transform=TR, fill=0, dtype='uint8')
    return r.astype(bool), geom


def fill_depressions(surf):
    seed = np.full(surf.shape, surf.max(), dtype='f4')
    seed[0, :] = surf[0, :]; seed[-1, :] = surf[-1, :]
    seed[:, 0] = surf[:, 0]; seed[:, -1] = surf[:, -1]
    return reconstruction(seed, surf, method='erosion').astype('f4')


def descend(fill, r, c, maxsteps=400):
    """Steepest descent on the filled surface; returns the lowest elevation reached.
    A rim cell on a rising shore cannot descend and returns its own elevation."""
    ny_, nx_ = fill.shape
    cur = float(fill[r, c])
    for _ in range(maxsteps):
        r0, r1 = max(r-1, 0), min(r+2, ny_)
        c0, c1 = max(c-1, 0), min(c+2, nx_)
        win = fill[r0:r1, c0:c1]
        dr, dc = np.unravel_index(int(np.argmin(win)), win.shape)
        nr, nc = r0+dr, c0+dc
        if float(fill[nr, nc]) >= cur - 1e-4:
            break
        r, c, cur = nr, nc, float(fill[nr, nc])
        if r in (0, ny_-1) or c in (0, nx_-1):
            break
    return cur


def polygonise(mask):
    out = []
    for geom, val in shapes(mask.astype('uint8'), mask=mask, transform=TR):
        g = shape(geom)
        a = g.area*(111.32**2)*np.cos(np.radians(g.centroid.y))
        if a < MIN_AREA_KM2:
            continue
        gs = g.simplify(SIMPLIFY, preserve_topology=True)
        if gs.is_empty or gs.area == 0:
            continue
        out.append((gs, a))
    return out


def run(t_yr, label):
    t0 = time.time()
    if t_yr > 0:
        gia, esl = gia_grid(t_yr)
    else:
        gia, esl = np.zeros_like(dem), 0.0
    ice, icegeom = ice_mask(label) if label else (np.zeros((ny, nx), bool), None)
    surf = (dem - gia).astype('f4')
    del gia

    work = surf.copy(); work[ice] = 9000.0
    filled = fill_depressions(work)
    water = ((filled - work) > 0.05) & (~ice)
    del work

    slab, _ = ndi.label(surf <= esl)
    edge = set(slab[:, -1].tolist()) | set(slab[-1, :].tolist()); edge.discard(0)
    marine = np.isin(slab, sorted(edge)) if edge else np.zeros_like(water)
    marine &= ~ice          # ice-covered ground is not open water
    del slab
    # A basin that holds a lake standing above sea level is fresh, whatever its
    # floor elevation. Only ground with no lake surface over it can be sea, so
    # the marine incursion is trimmed to what is left, not the other way round.
    marine &= ~water

    llab, nl = ndi.label(water)
    sizes = ndi.sum(np.broadcast_to(AREA_ROW[:, None], llab.shape).astype('f4'),
                    llab, np.arange(1, nl+1))
    keep = (np.where(sizes >= MIN_AREA_KM2)[0] + 1).astype(np.int32)
    big = np.isin(llab, keep)

    lakes = []
    objs = ndi.find_objects(llab)
    order_by_area = keep[np.argsort(sizes[keep-1])[::-1]]
    for rank, lid in enumerate(order_by_area):
        sl = objs[lid-1]
        r0 = max(sl[0].start-1, 0); r1 = min(sl[0].stop+1, ny)
        c0 = max(sl[1].start-1, 0); c1 = min(sl[1].stop+1, nx)
        s2 = (slice(r0, r1), slice(c0, c1))
        sub = (llab[s2] == lid)
        sel_full = (llab == lid)
        level = float(np.median(filled[s2][sub]))
        dry = (filled[s2] - surf[s2]) <= 0.01          # not inside any depression
        ring = (ndi.binary_dilation(sub, np.ones((3, 3), bool))
                & ~sub & ~ice[s2] & dry)
        spill = None
        if ring.any() and rank < 6:
            # The outlet is where the lake touches ground that genuinely gets
            # away: cells below the lake level that are connected to the domain
            # edge. A low rim cell in a severed channel stub, or one of many
            # tied cells on a flat shore, fails this test.
            esc = filled < (level - 0.02)
            seed = np.zeros_like(esc)
            seed[0, :] = esc[0, :]; seed[-1, :] = esc[-1, :]
            seed[:, 0] = esc[:, 0]; seed[:, -1] = esc[:, -1]
            E8 = np.ones((3, 3), bool)
            away = ndi.binary_propagation(seed, mask=esc, structure=E8)
            # walk back up any flat spillway reach, but never through the lake
            route = ndi.binary_propagation(away, mask=((filled <= level+0.02) & ~sel_full),
                                           structure=E8)
            cand = ring & ndi.binary_dilation(route[s2], E8)
            if not cand.any():
                cand = ring & (surf[s2] <= surf[s2][ring].min() + 0.02)
            rr, cc = np.where(cand)
            k = int(np.argmin(surf[s2][rr, cc]))
            spill = [round(float(lon_c[c0+cc[k]]), 4),
                     round(float(lat_c[r0+rr[k]]), 4),
                     round(float(surf[s2][rr[k], cc[k]]), 1)]
        lakes.append(dict(id=int(lid), area=float(sizes[lid-1]),
                          level=round(level, 1), spill=spill))
    lakes.sort(key=lambda d: -d['area'])
    print(f'  {label or "0ka":>7s}  esl {esl:6.1f}  ice {ice.mean()*100:4.1f}%  '
          f'lakes {len(lakes):3d}  largest {lakes[0]["area"] if lakes else 0:8.0f} km2'
          f'  [{time.time()-t0:.0f}s]', flush=True)
    return dict(t=t_yr, esl=esl, icegeom=icegeom, lakes=lakes, llab=llab,
                big=big, marine=marine, filled=filled, surf=surf, ice=ice)


if __name__ == '__main__':
    r = run(int(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else None)
    for lk in r['lakes'][:14]:
        s = lk['spill']
        print(f"    {lk['area']:9.0f} km2   level {lk['level']:7.1f} m   "
              + (f"spill {s[0]:8.3f},{s[1]:7.3f} @ {s[2]:7.1f} m" if s else "no spill"))
