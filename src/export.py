import numpy as np, json, gc, time
from scipy import ndimage as ndi
from shapely.geometry import shape, mapping, box
from shapely.ops import unary_union
from rasterio.features import shapes
import model as M

STEPS = [(t, f'{t/1000:g}ka') for t in range(14500, 500, -500)] + [(0, None)]
CLIP = box(M.W, M.S, M.E, M.N)
LAKE_SIMP, ICE_SIMP = 0.012, 0.030
MIN_SHOW = 400.0


def rings(geom, simp, prec=3):
    g = geom.simplify(simp, preserve_topology=True)
    if g.is_empty:
        return []
    polys = [g] if g.geom_type == 'Polygon' else list(g.geoms)
    out = []
    for p in polys:
        if p.geom_type != 'Polygon':
            continue
        ext = [[round(x, prec), round(y, prec)] for x, y in p.exterior.coords]
        if len(ext) < 4:
            continue
        holes = []
        for r in p.interiors:
            h = [[round(x, prec), round(y, prec)] for x, y in r.coords]
            if len(h) >= 4 and abs(shape({'type': 'Polygon', 'coordinates': [h]}).area) > 3e-3:
                holes.append(h)
        out.append([ext]+holes)
    return out


def mask_to_geom(mask, min_km2):
    parts = []
    for geom, _ in shapes(mask.astype('uint8'), mask=mask, transform=M.TR):
        g = shape(geom)
        a = g.area*(111.32**2)*np.cos(np.radians(g.centroid.y))
        if a >= min_km2:
            parts.append(g)
    return unary_union(parts) if parts else None


# modern basin masks, used to name the palaeolakes
base = M.run(0, None)
basins = {}
for nm, (lo, la) in {'Superior': (-87.5, 47.6), 'Michigan': (-86.8, 44.0),
                     'Huron': (-82.3, 44.6), 'Georgian Bay': (-80.9, 45.4),
                     'Erie': (-81.4, 42.1), 'Ontario': (-77.5, 43.6)}.items():
    r = int((M.N-la)/M.RES); c = int((lo-M.W)/M.RES)
    lid = base['llab'][r, c]
    basins[nm] = (base['llab'] == lid) if lid else np.zeros(base['llab'].shape, bool)
# split Michigan-Huron-Georgian, which is one lake today, by longitude/latitude cuts
mh = basins['Michigan']
LONG, LATG = np.meshgrid(M.lon_c, M.lat_c)
basins['Michigan'] = mh & (LONG < -84.75)
basins['Huron'] = mh & (LONG >= -84.75) & ~((LONG > -81.4) & (LATG < 45.6))
basins['Georgian Bay'] = mh & (LONG > -81.4) & (LATG < 45.9)
del base, mh
gc.collect()

frames = []
for t, label in STEPS:
    r = M.run(t, label)
    surf, filled, llab = r['surf'], r['filled'], r['llab']
    feats = []
    for lk in r['lakes']:
        if lk['area'] < MIN_SHOW:
            continue
        sel = (llab == lk['id'])
        g = mask_to_geom(sel, MIN_SHOW*0.35)
        if g is None:
            continue
        names = [nm for nm, bm in basins.items()
                 if (sel & bm).sum() > 0.18*bm.sum()]
        feats.append(dict(r=rings(g, LAKE_SIMP), lv=lk['level'],
                          a=round(lk['area']), b=names, sp=lk['spill']))
        del sel
    icegeom = r['icegeom'].intersection(CLIP) if r['icegeom'] is not None else None
    mar = mask_to_geom(r['marine'], 400.0) if r['marine'].any() else None
    frames.append(dict(t=t, esl=round(r['esl'], 1),
                       lakes=feats,
                       ice=rings(icegeom, ICE_SIMP) if icegeom else [],
                       sea=rings(mar, ICE_SIMP) if mar else []))
    del r, surf, filled, llab
    gc.collect()

# modern shoreline reference, from the 0 ka solution
json.dump(dict(res_arcsec=30, bbox=[M.W, M.S, M.E, M.N], frames=frames),
          open('out/paleo_%s.json' % M.MARGIN.lower(), 'w'), separators=(',', ':'))
import os
print('wrote out/paleo_'+M.MARGIN.lower()+'.json  %.2f MB' % (os.path.getsize('out/paleo_%s.json' % M.MARGIN.lower())/1e6))
