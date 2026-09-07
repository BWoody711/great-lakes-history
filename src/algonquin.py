"""
The Lake Algonquin water plane, imposed rather than solved.

The Algonquin strandline is horizontal at about 184 m in the zone south of the
hinge line, and rises north-eastward from there. That tilt is differential
uplift, so if the modern strandline elevation at x is

    z_now(x) = 184 + [ U(x,t) - U(ref,t) ]

with U the uplift since time t, then in the palaeo frame used by the model,
surf = dem - U, the same water surface is simply flat:

    z_then(x) = z_now(x) - U(x,t) = 184 - U(ref,t)

So the known lake is one horizontal flood of the palaeo surface, at a level set
by the uplift at the reference locality alone. Port Huron is the reference.
This is the geologically known lake drawn on the model's own topography, which
is what makes it comparable with the lake the model solves for.
"""
import numpy as np, json, gc
from scipy import ndimage as ndi
import model as M
from export import mask_to_geom, rings, LAKE_SIMP
from shapely.geometry import Polygon
from rasterio.features import rasterize

ALG_M = 184.0                                  # mapped strandline, hinge-free zone
REF = (-82.42, 43.00)                          # Port Huron
STEPS = [13000, 12500, 12000, 11500]           # Algonquin interval; main phase 12.6-11.5
SEEDS = [(-82.3, 44.6), (-86.8, 44.0), (-80.9, 45.3), (-84.5, 45.6)]
# The lake's own outlets. Water at the Algonquin plane leaves through these, so
# the flood must not be allowed to walk out through them into the basins below;
# plugging the sill is what confines the plane to the lake it belongs to.
SILLS = [(-82.42, 43.00, 'Port Huron'), (-87.80, 41.72, 'Chicago'),
         (-78.74, 44.54, 'Kirkfield'), (-79.45, 46.32, 'North Bay'),
         (-78.20, 44.30, 'Trent'), (-84.35, 46.50, 'Sault')]
PLUG_KM = 12.0


def sill_plug():
    m = np.zeros((M.ny, M.nx), bool)
    for lo, la, _ in SILLS:
        dlat = PLUG_KM/111.32
        dlon = PLUG_KM/(111.32*np.cos(np.radians(la)))
        r0 = int((M.N-(la+dlat))/M.RES); r1 = int((M.N-(la-dlat))/M.RES)
        c0 = int(((lo-dlon)-M.W)/M.RES); c1 = int(((lo+dlon)-M.W)/M.RES)
        m[max(r0, 0):r1, max(c0, 0):c1] = True
    return m


PLUG = sill_plug()

# The plane on its own floods east through the deeply depressed Ontario and
# Ottawa basins, which is the same escape that suppresses Algonquin in the
# solved run. The known lake did not occupy those basins, so the overlay is
# confined to the upper Great Lakes basin. This clip fixes which basins the
# lake may occupy; the shoreline position inside it is still set entirely by
# the topography and the water plane.
BASIN = Polygon([(-93.2, 46.2), (-92.6, 49.6), (-88.0, 49.7), (-84.0, 48.3),
                 (-81.0, 47.3), (-79.3, 46.7), (-78.5, 45.3), (-79.1, 43.9),
                 (-80.4, 43.2), (-82.1, 42.6), (-82.9, 42.7), (-83.5, 43.7),
                 (-85.0, 43.1), (-86.0, 41.4), (-88.0, 41.4), (-88.7, 42.7),
                 (-90.6, 44.1), (-92.1, 45.2)])
INSIDE = rasterize([(BASIN, 1)], out_shape=(M.ny, M.nx), transform=M.TR,
                   fill=0, dtype='uint8').astype(bool)

out = {}
for margin in ('OPTIMAL', 'MAX'):
    M.MARGIN = margin
    fr = []
    for t in STEPS:
        gia, esl = M.gia_grid(t)
        ice, _ = M.ice_mask(f'{t/1000:g}ka')
        surf = (M.dem - gia).astype('f4')
        uref = float(gia[int((M.N-REF[1])/M.RES), int((REF[0]-M.W)/M.RES)])
        level = ALG_M - uref
        wet = (surf <= level) & ~ice & ~PLUG & INSIDE
        lab, n = ndi.label(wet)
        keep = set()
        for lo, la in SEEDS:
            v = lab[int((M.N-la)/M.RES), int((lo-M.W)/M.RES)]
            if v:
                keep.add(int(v))
        if keep:
            mask = np.isin(lab, sorted(keep))
            area = float((np.broadcast_to(M.AREA_ROW[:, None], mask.shape)*mask).sum())
            geom = mask_to_geom(mask, 400.0)
            r = rings(geom, LAKE_SIMP) if geom else []
        else:
            area, r = 0.0, []
        fr.append(dict(t=t, level=round(level, 1), uref=round(uref, 1),
                       area=round(area), r=r))
        print(f'{margin:8s} {t/1000:5g} ka  U(Port Huron) {uref:7.1f} m  '
              f'plane {level:6.1f} m  area {area:9,.0f} km2', flush=True)
        del gia, ice, surf, wet, lab
        gc.collect()
    out[margin.lower()] = fr

json.dump(dict(strandline=ALG_M, ref=REF, frames=out),
          open('out/algonquin.json', 'w'), separators=(',', ':'))
import os
print('wrote out/algonquin.json  %.0f kB' % (os.path.getsize('out/algonquin.json')/1e3))
