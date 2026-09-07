"""
Great Lakes shoreline with shoreline-type attributes, v2.

Changes from v1:
  * base geometry is the GSHHG full-resolution lake outlines, not the ETOPO
    raster boundaries.  The raster boundaries stair-step at the cell size, so
    every segment orientation was quantised to multiples of 45 degrees, which
    is useless if you care about shore-normal or drift direction.
  * the Canadian shore is typed from Environment and Climate Change Canada's
    Great Lakes Nearshore Waters Assessment, which NOAA's ESI does not cover.

Two typologies, kept in separate fields because they are not the same thing:
  ESI   - substrate and form (sand, gravel, bedrock, wetland, armored, ...)
  ECCC  - hardened vs natural, plus a wave-energy class for the reach
"""
import geopandas as gpd, pandas as pd, numpy as np, glob, struct, subprocess
from shapely.geometry import LineString, MultiLineString
from shapely.ops import substring
from io import StringIO

GLA = 'EPSG:3175'
SEG = 500.0
MAXD_ESI = 2000.0
MAXD_ECCC = 3000.0
OUT = '/mnt/user-data/outputs/great_lakes_shapefiles'
D = '/home/claude/data'

ESI_ATLASES = [
    ('Lake Michigan',       2025, f'{D}/esi/Lake_Michigan_ESI_2025_GDB/**/*.gdb'),
    ('Lake Ontario',        2023, f'{D}/esi/Lake_Ontario_ESI_2023_GDB/**/*.gdb'),
    ('Lake Erie',           2022, f'{D}/esi_erie/**/*.gdb'),
    ('St Marys R.',         2021, f'{D}/esi/Great_Lakes_St_Marys_River_2021_ESI_GDB/**/*.gdb'),
    ('St Lawrence R.',      2021, f'{D}/esi/Great_Lakes_St_Lawrence_River_2021_ESI_GDB/**/*.gdb'),
    ('Straits of Mackinac', 2019, f'{D}/esi/GL_Straits_of_Mackinac_2019_GDB/**/*.gdb'),
    ('St Clair/Detroit',    2019, f'{D}/esi/GL_StClair_Detroit_River_System_2019_GDB/**/*.gdb'),
]
HURON_MDB = (f'{D}/esi/LakeHuron_1994_GDB/LakeHuron_1994_GDB/LakeHuron/LakeHuronESI.mdb')
SCAT = f'{D}/scat/ShorelineClassification_ON_OpenDataCatalogue.gdb'
SCAT_MAP = {
    'Bedrock Cliff/Vertical': 'bedrock', 'Bedrock Platform': 'bedrock',
    'Pebble/Cobble Beach or Bank': 'gravel', 'Boulder Beach or Bank': 'gravel',
    'Mixed Sediment Beach or Bank': 'mixed sand and gravel',
    'Sand Beach or Bank': 'sand', 'Sediment Cliff': 'eroding bank (unconsol.)',
    'Marsh': 'wetland', 'Mud Tidal Flat': 'flats',
    'Man-Made Solid': 'armored', 'Man-Made Permeable': 'armored',
    'Not Classified': None}
MAXD_SCAT = 2000.0

ECCC = [
    ('Erie',     f'{D}/eccc/erie/LakeErieNearshoreWatersAssessment_z17N.gdb',
     'LE_CoastalProcesses_Shoreline_Hardening'),
    ('Ontario',  f'{D}/eccc/ontario/LakeOntarioNearshoreWatersAssessment_z17N.gdb',
     'LO_CoastalProcesses_Shoreline_Hardening'),
    ('Huron',    f'{D}/eccc/huron/LakeHuronNearshoreWatersAssessment_z17N.gdb',
     'LH_CoastalProcesses_Shoreline_Hardening'),
    ('Superior', f'{D}/eccc/superior/LakeSuperiorNearshoreWatersAssessment_z16N.gdb',
     'LS_CoastalProcesses_Shoreline_Hardening'),
]


def substrate(code):
    if not isinstance(code, str) or not code:
        return None
    c = code.split('/')[0].strip().upper()
    num = ''.join(ch for ch in c if ch.isdigit())
    let = ''.join(ch for ch in c if ch.isalpha())
    if not num:
        return None
    n = int(num)
    return {1: 'armored' if let == 'B' else 'bedrock',
            2: 'bedrock',
            3: 'eroding bank (unconsol.)' if let == 'B' else 'sand',
            4: 'sand',
            5: 'mixed sand and gravel',
            6: 'armored' if let == 'B' else 'gravel',
            7: 'flats',
            8: 'armored' if let in ('B', 'C') else ('veg. bank or bluff' if let == 'F' else 'bedrock'),
            9: 'veg. bank or bluff' if let == 'B' else 'flats',
            10: 'wetland'}.get(n)


def read_huron():
    csv = subprocess.run(['mdb-export', '-b', 'hex', HURON_MDB, 'esil_arcs'],
                         capture_output=True, text=True).stdout
    df = pd.read_csv(StringIO(csv))
    geoms, keep = [], []
    for i, blob in enumerate(df['Shape'].fillna('')):
        try:
            b = bytes.fromhex(blob)
            if len(b) < 44 or struct.unpack('<i', b[:4])[0] != 3:
                continue
            nparts, npts = struct.unpack('<ii', b[36:44])
            off = 44 + 4*nparts
            parts = list(struct.unpack(f'<{nparts}i', b[44:off])) + [npts]
            xy = np.frombuffer(b[off:off+16*npts], dtype='<f8').reshape(npts, 2)
            ls = [LineString(xy[parts[k]:parts[k+1]]) for k in range(nparts)
                  if parts[k+1]-parts[k] >= 2]
            if not ls:
                continue
            geoms.append(ls[0] if len(ls) == 1 else MultiLineString(ls))
            keep.append(i)
        except Exception:
            continue
    g = gpd.GeoDataFrame(df.iloc[keep].copy(), geometry=geoms, crs='EPSG:4326')
    g['atlas'], g['year'] = 'Lake Huron', 1994
    g['esi'] = g['ESI'].astype(str)
    g['exposure'] = g['ENVIR'].astype(str)
    g['shoretype'] = None
    print(f'  Lake Huron          1994  {len(g):6d} arcs (parsed from .mdb)')
    return g[['esi', 'shoretype', 'exposure', 'atlas', 'year', 'geometry']]


def load_esil(name, year, pattern):
    p = glob.glob(pattern, recursive=True)[0]
    g = gpd.read_file(p, layer='ESIL').to_crs('EPSG:4326')
    cols = {c.upper(): c for c in g.columns}
    g['esi'] = g[cols['ESI']].astype(str)
    g['shoretype'] = g[cols['LANDWARD_SHORETYPE']] if 'LANDWARD_SHORETYPE' in cols else None
    g['exposure'] = g[cols['ENVIR']].astype(str) if 'ENVIR' in cols else None
    g['atlas'], g['year'] = name, year
    print(f'  {name:20s}{year}  {len(g):6d} segments')
    return g[['esi', 'shoretype', 'exposure', 'atlas', 'year', 'geometry']]


print('NOAA ESI atlases')
esi = pd.concat([load_esil(*a) for a in ESI_ATLASES] + [read_huron()], ignore_index=True)
esi = gpd.GeoDataFrame(esi, geometry='geometry', crs='EPSG:4326').to_crs(GLA)
lut = (esi.dropna(subset=['shoretype']).groupby('esi')['shoretype']
          .agg(lambda s: s.value_counts().index[0]))
esi['shoretype'] = esi['shoretype'].fillna(esi['esi'].map(lut))
esi['substrate'] = esi['esi'].map(substrate)

print('\nECCC Great Lakes Nearshore Waters Assessment')
ec = []
for lake, path, layer in ECCC:
    g = gpd.read_file(path, layer=layer).to_crs(GLA)
    g['eccc_hard'] = g['ShType'].astype(str)
    g['eccc_unit'] = g['UnitType'].astype(str)
    g['eccc_year'] = pd.to_numeric(g.get('AssessmentYear'), errors='coerce')
    print(f'  {lake:20s}{int(g["eccc_year"].median()):4d}  {len(g):6d} reaches, '
          f'{g.length.sum()/1000:6.0f} km, '
          f'{(g["eccc_hard"]=="Hardened").mean()*100:.0f}% hardened')
    ec.append(g[['eccc_hard', 'eccc_unit', 'eccc_year', 'geometry']])
ec = gpd.GeoDataFrame(pd.concat(ec, ignore_index=True), geometry='geometry', crs=GLA)

print('\nECCC Shoreline Classification (SCAT), Ontario')
sc = gpd.read_file(SCAT, layer='O14Oceans_ShorelineClass_ON').to_crs(GLA)
sc['scat_cls'] = sc['SCAT_Class_EN'].astype(str)
sc['scat_sub'] = sc['scat_cls'].map(SCAT_MAP)
print(f'  {len(sc):,} segments province-wide, {sc.length.sum()/1000:,.0f} km; '
      f'{sc["scat_sub"].notna().sum():,} carry a substrate class')

# ---------------------------------------------------- base geometry: GSHHG
lakes = gpd.read_file(f'{OUT}/great_lakes_gshhg_full.shp').to_crs(GLA)
rows = []
for _, lk in lakes.iterrows():
    polys = [lk.geometry] if lk.geometry.geom_type == 'Polygon' else list(lk.geometry.geoms)
    for poly in polys:
        for ring in [poly.exterior] + list(poly.interiors):
            n = max(1, int(round(ring.length/SEG)))
            step = ring.length/n
            for k in range(n):
                s = substring(ring, k*step, (k+1)*step)
                if s.length > 1:
                    rows.append(dict(lake=lk['name'], geometry=s))
shore = gpd.GeoDataFrame(rows, crs=GLA)
shore['seg_km'] = shore.length/1000

# segment bearing, which is the point of using vector geometry
def bearing(g):
    a, b = g.coords[0], g.coords[-1]
    return (np.degrees(np.arctan2(b[0]-a[0], b[1]-a[1]))) % 180.0
shore['bearing'] = shore.geometry.apply(bearing).round(1)
print(f'\nGSHHG shoreline segmented: {len(shore):,} segments, {shore.seg_km.sum():,.0f} km')

mid = shore.copy(); mid['geometry'] = shore.geometry.interpolate(0.5, normalized=True)

# ---------------------------------------------------------------- ESI join
for c in ['esi', 'shoretype', 'substrate', 'exposure', 'atlas', 'year', 'esi_m']:
    shore[c] = np.nan if c in ('year', 'esi_m') else None
for name, year, _ in ESI_ATLASES + [('Lake Huron', 1994, None)]:
    todo = shore.index[shore['atlas'].isna()]
    if len(todo) == 0:
        break
    src = esi[esi['atlas'] == name]
    if not len(src):
        continue
    j = gpd.sjoin_nearest(mid.loc[todo, ['geometry']], src, how='left',
                          max_distance=MAXD_ESI, distance_col='esi_m')
    j = j[~j.index.duplicated(keep='first')]
    hit = j['atlas'].notna()
    for c in ['esi', 'shoretype', 'substrate', 'exposure', 'atlas', 'year', 'esi_m']:
        shore.loc[j.index[hit], c] = j.loc[hit, c].values

# --------------------------------------------------------------- SCAT join
sc_ok = sc[sc['scat_sub'].notna()]
j = gpd.sjoin_nearest(mid[['geometry']], sc_ok[['scat_cls', 'scat_sub', 'geometry']],
                      how='left', max_distance=MAXD_SCAT, distance_col='scat_m')
j = j[~j.index.duplicated(keep='first')]
for c in ['scat_cls', 'scat_sub', 'scat_m']:
    shore[c] = j[c].values

# --------------------------------------------------------------- ECCC join
j = gpd.sjoin_nearest(mid[['geometry']], ec, how='left',
                      max_distance=MAXD_ECCC, distance_col='eccc_m')
j = j[~j.index.duplicated(keep='first')]
for c in ['eccc_hard', 'eccc_unit', 'eccc_year', 'eccc_m']:
    shore[c] = j[c].values

# ------------------------------------------------- combined substrate field
shore['sub_src'] = np.where(shore['substrate'].notna(), 'NOAA ESI', None)
f2 = shore['substrate'].isna() & shore['scat_sub'].notna()
shore.loc[f2, 'substrate'] = shore.loc[f2, 'scat_sub']
shore.loc[f2, 'sub_src'] = 'ECCC SCAT'
f3 = shore['substrate'].isna() & (shore['eccc_hard'] == 'Hardened')
shore.loc[f3, 'substrate'] = 'armored'
shore.loc[f3, 'sub_src'] = 'ECCC (hardened only)'
# ECCC "Natural" says only that it is not hardened, so no substrate is claimed

# ------------------------------------------- fused drift class
# One general attribute, built for a single question: can this shore take part
# in alongshore drift?  Derived from the source codes, not from the collapsed
# substrate field, because a couple of ESI ranks collapse ambiguously.
ESI_DRIFT = {
    '1A': 'bedrock', '1C': 'bedrock', '2A': 'bedrock', '2': 'bedrock', '8A': 'bedrock',
    '1B': 'armoured', '6B': 'armoured', '8B': 'armoured', '8C': 'armoured',
    '3A': 'sediment', '3B': 'sediment', '3': 'sediment', '4': 'sediment',
    '5': 'sediment', '6A': 'sediment', '8F': 'sediment',
    '7': 'low energy', '9A': 'low energy', '9B': 'low energy',
    '10A': 'low energy', '10B': 'low energy', '10C': 'low energy',
    '10D': 'low energy', '10': 'low energy',
}
SCAT_DRIFT = {
    'Bedrock Cliff/Vertical': 'bedrock', 'Bedrock Platform': 'bedrock',
    'Pebble/Cobble Beach or Bank': 'sediment', 'Boulder Beach or Bank': 'sediment',
    'Mixed Sediment Beach or Bank': 'sediment', 'Sand Beach or Bank': 'sediment',
    'Sediment Cliff': 'sediment',
    'Marsh': 'low energy', 'Mud Tidal Flat': 'low energy',
    'Man-Made Solid': 'armoured', 'Man-Made Permeable': 'armoured',
}


def drift_class(row):
    e = row['esi']
    if isinstance(e, str) and e:
        c = ESI_DRIFT.get(e.split('/')[0].strip().upper())
        if c:
            return c, 'NOAA ESI'
    sc = row['scat_cls']
    if isinstance(sc, str) and sc in SCAT_DRIFT:
        return SCAT_DRIFT[sc], 'ECCC SCAT'
    if row['eccc_hard'] == 'Hardened':
        return 'armoured', 'ECCC hardening'
    return None, None


dc = shore.apply(drift_class, axis=1, result_type='expand')
shore['drift_cls'] = dc[0]
shore['drift_src'] = dc[1]
# does the shore hold mobile clastic sediment that can move alongshore?
shore['mobile'] = np.where(shore['drift_cls'].isna(), None,
                           np.where(shore['drift_cls'] == 'sediment', 'Y', 'N'))

got = shore['substrate'].notna()
any_info = got | shore['eccc_hard'].notna()
print(f'\nsubstrate assigned : {shore.loc[got,"seg_km"].sum():,.0f} km '
      f'({100*shore.loc[got,"seg_km"].sum()/shore.seg_km.sum():.0f}%)')
print(f'any type information: {shore.loc[any_info,"seg_km"].sum():,.0f} km '
      f'({100*shore.loc[any_info,"seg_km"].sum()/shore.seg_km.sum():.0f}%)')

shore = shore.to_crs('EPSG:4326')
shore[['lake', 'esi', 'shoretype', 'substrate', 'sub_src', 'exposure', 'atlas',
       'year', 'scat_cls', 'drift_cls', 'drift_src', 'mobile',
       'eccc_hard', 'eccc_unit', 'eccc_year', 'bearing',
       'seg_km', 'esi_m', 'scat_m', 'eccc_m', 'geometry']].to_file(f'{OUT}/great_lakes_shoreline_typed.shp')

d = (shore.groupby(['lake', shore['drift_cls'].fillna('unknown')])['seg_km']
          .sum().unstack(fill_value=0).round(0))
order = [c for c in ['bedrock', 'sediment', 'armoured', 'low energy', 'unknown'] if c in d]
d = d[order]
d.to_csv(f'{OUT}/great_lakes_drift_class_summary.csv')
print('\nfused drift class, km by lake')
print(d.to_string())
print('\ntotals: ' + ', '.join(f'{k} {v:,.0f} km ({100*v/d.values.sum():.0f}%)'
                                for k, v in d.sum().items()))

t = (shore[got].groupby(['lake', 'substrate'])['seg_km'].sum().unstack(fill_value=0))
t['NATURAL, type unknown'] = shore[(~got) & (shore['eccc_hard'] == 'Natural')] \
    .groupby('lake')['seg_km'].sum()
t['NO DATA'] = shore[~any_info].groupby('lake')['seg_km'].sum()
t = t.fillna(0).round(0)
t.to_csv(f'{OUT}/great_lakes_shoreline_type_summary.csv')
print()
print(t.to_string())
