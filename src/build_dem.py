import numpy as np, netCDF4, os

W, E, S, N = -95.0, -70.0, 38.0, 53.0
RES = 1/240.0                       # 15 arc-seconds
nx = int(round((E-W)/RES)); ny = int(round((N-S)/RES))
print("target grid", nx, "x", ny)

# north-up array, row 0 = north edge
dem = np.full((ny, nx), np.nan, dtype=np.float32)
lon_c = W + (np.arange(nx)+0.5)*RES
lat_c = N - (np.arange(ny)+0.5)*RES

for tile in ["N45W105","N45W090","N45W075","N60W105","N60W090","N60W075"]:
    f = f"data/etopo_{tile}.nc"
    d = netCDF4.Dataset(f)
    tlat = d['lat'][:].astype('f8'); tlon = d['lon'][:].astype('f8')
    # index ranges of the target grid that fall inside this tile
    ci = np.where((lon_c >= tlon.min()-RES/2) & (lon_c <= tlon.max()+RES/2))[0]
    ri = np.where((lat_c >= tlat.min()-RES/2) & (lat_c <= tlat.max()+RES/2))[0]
    if len(ci)==0 or len(ri)==0:
        d.close(); continue
    # nearest source indices (grids are co-registered so this is exact)
    si_x = np.round((lon_c[ci]-tlon[0])/RES).astype(int)
    si_y = np.round((lat_c[ri]-tlat[0])/RES).astype(int)      # tlat ascending
    ok_x = (si_x>=0)&(si_x<len(tlon)); ok_y = (si_y>=0)&(si_y<len(tlat))
    ci, si_x = ci[ok_x], si_x[ok_x]; ri, si_y = ri[ok_y], si_y[ok_y]
    block = d['z'][si_y.min():si_y.max()+1, si_x.min():si_x.max()+1]
    block = np.asarray(block, dtype=np.float32)
    block = block[si_y-si_y.min(), :][:, si_x-si_x.min()]
    dem[ri.min():ri.max()+1, ci.min():ci.max()+1] = block
    print("  merged", tile, block.shape, "z range %.0f..%.0f" % (np.nanmin(block), np.nanmax(block)))
    d.close()

print("nan cells:", int(np.isnan(dem).sum()))
np.save("data/dem15.npy", dem)

# sanity probes
def at(lon, lat):
    r = int((N-lat)/RES); c = int((lon-W)/RES); return float(dem[r,c])
for name, lon, lat in [("Superior mid-basin", -87.5, 47.6), ("Huron deep (Goderich)", -82.2, 44.3),
                       ("Michigan deep", -86.8, 44.2), ("Erie central", -81.6, 42.0),
                       ("Ontario deep", -76.9, 43.6), ("Chicago sill area", -87.8, 41.75),
                       ("Port Huron", -82.42, 43.00), ("North Bay", -79.45, 46.32),
                       ("Rome NY", -75.45, 43.21), ("Fort Wayne", -85.13, 41.08)]:
    print(f"  {name:24s} {at(lon,lat):8.1f} m")
