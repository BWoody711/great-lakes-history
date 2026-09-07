"""
Modern Great Lakes polygons, three ways.

  1. GSHHG full resolution (Wessel & Smith), split at the connecting straits
  2. Natural Earth 1:10m physical lakes, already one feature per lake
  3. Derived here by flooding the NOAA ETOPO 2022 15" grid to each lake datum

All in EPSG:4326. Areas are computed on an equal-area projection (EPSG:6933),
not from degrees.
"""
import numpy as np, geopandas as gpd, os
from shapely.geometry import LineString, Polygon, MultiPolygon, box
from shapely.ops import unary_union
from scipy import ndimage as ndi
from rasterio.features import shapes
from rasterio.transform import from_origin
from shapely.geometry import shape as shp

OUT = '/mnt/user-data/outputs/great_lakes_shapefiles'
os.makedirs(OUT, exist_ok=True)
EA = 'EPSG:6933'

# lake seeds and their chart datums (m, IGLD 1985)
LAKES = [('Superior',      -87.50, 47.60, 183.5),
         ('Michigan',      -86.80, 44.00, 176.0),
         ('Huron',         -82.30, 44.60, 176.0),
         ('St. Clair',     -82.70, 42.45, 174.9),
         ('Erie',          -81.40, 42.10, 174.0),
         ('Ontario',       -77.50, 43.60,  74.2)]

# narrow cuts across the channels that join the basins, so the lakes separate
CUTS = [LineString([(-84.76, 45.68), (-84.76, 45.98)]),      # Straits of Mackinac
        LineString([(-84.52, 46.40), (-84.22, 46.62)]),      # St Marys R., Sault
        LineString([(-83.28, 42.14), (-82.94, 42.14)]),      # Detroit River
        LineString([(-82.60, 42.72), (-82.30, 42.72)]),      # St Clair River
        LineString([(-79.12, 42.86), (-78.86, 43.10)])]      # Niagara River
CUT = unary_union([c.buffer(0.004) for c in CUTS])


def split_and_name(geoms):
    parts = []
    for g in geoms:
        d = g.difference(CUT)
        parts += [d] if d.geom_type == 'Polygon' else list(d.geoms)
    rows = []
    for name, lo, la, lvl in LAKES:
        hit = [p for p in parts if p.contains(__import__('shapely').geometry.Point(lo, la))]
        if not hit:
            print('  no polygon found for', name); continue
        g = max(hit, key=lambda p: p.area)
        rows.append(dict(name=name, datum_m=lvl, geometry=g))
    return rows


def finish(rows, path, source, note):
    gdf = gpd.GeoDataFrame(rows, crs='EPSG:4326')
    gdf['area_km2'] = (gdf.to_crs(EA).area/1e6).round(0)
    gdf['nvert'] = gdf.geometry.apply(
        lambda g: sum(len(p.exterior.coords) for p in
                      ([g] if g.geom_type == 'Polygon' else g.geoms)))
    gdf['source'] = source
    gdf = gdf[['name', 'datum_m', 'area_km2', 'nvert', 'source', 'geometry']]
    gdf.to_file(path)
    print(f'\n{os.path.basename(path)}   ({note})')
    print(gdf.drop(columns='geometry').to_string(index=False))
    return gdf


# ---- 1. GSHHG full resolution -------------------------------------------
g = gpd.read_file('data/gshhg/GSHHS_shp/f/GSHHS_f_L2.shp', bbox=(-93, 41, -75, 50))
g = g[g.area > 0.05]
a = finish(split_and_name(list(g.geometry.values)),
           f'{OUT}/great_lakes_gshhg_full.shp',
           'GSHHG 2.3.7 full res (Wessel & Smith), LGPL',
           'highest detail')

# ---- 2. Natural Earth ----------------------------------------------------
ne = gpd.read_file('data/ne_lakes/ne_10m_lakes.shp')
ne = ne[ne['name_alt'].astype(str).str.contains('Great Lakes', na=False)]
rows = []
for name, lo, la, lvl in LAKES:
    from shapely.geometry import Point
    hit = ne[ne.geometry.contains(Point(lo, la))]
    if len(hit):
        rows.append(dict(name=name, datum_m=lvl, geometry=hit.iloc[0].geometry))
b = finish(rows, f'{OUT}/great_lakes_naturalearth.shp',
           'Natural Earth 1:10m physical lakes, public domain', 'generalised')

# ---- 3. Flooded from the NOAA ETOPO 2022 15" grid ------------------------
W, E, S, N, RES = -95.0, -70.0, 38.0, 53.0, 1/240.0
dem = np.load('data/dem15.npy')
ny, nx = dem.shape
TR = from_origin(W, N, RES, RES)
lon_c = W + (np.arange(nx)+0.5)*RES
lat_c = N - (np.arange(ny)+0.5)*RES
rows = []
for name, lo, la, lvl in LAKES:
    wet = dem <= lvl
    lab, _ = ndi.label(wet)
    lid = lab[int((N-la)/RES), int((lo-W)/RES)]
    m = (lab == lid)
    del lab, wet
    polys = [shp(gm) for gm, _ in shapes(m.astype('uint8'), mask=m, transform=TR)]
    gm = max(polys, key=lambda p: p.area)
    gm = gm.simplify(0.0015, preserve_topology=True).difference(CUT)
    if gm.geom_type == 'MultiPolygon':
        from shapely.geometry import Point
        cand = [p for p in gm.geoms if p.contains(Point(lo, la))]
        gm = cand[0] if cand else max(gm.geoms, key=lambda p: p.area)
    rows.append(dict(name=name, datum_m=lvl, geometry=gm))
    del m
c = finish(rows, f'{OUT}/great_lakes_etopo15s.shp',
           'flooded from ETOPO 2022 15" (NOAA NCEI) at chart datum', '~450 m cells')

print('\npublished areas for comparison:')
for n, v in [('Superior', 82100), ('Michigan', 57800), ('Huron', 59600),
             ('St. Clair', 1114), ('Erie', 25700), ('Ontario', 18960)]:
    print(f'  {n:12s} {v:8,d} km2')
