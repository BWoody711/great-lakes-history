"""
Observational calibration of the GIA field over the Great Lakes.

ICE-6G_C / ICE-7G_NA are global inversions. Their Great Lakes tilt is a long
way from their data: neither was fitted to lake-level gauges (Peltier, pers.
comm. in Mainville & Craymer 2005), and the region sits on the shoulder of the
Laurentide load where the answer is most sensitive to mantle viscosity and
lithospheric thickness. Three independent observational studies say the raw
model tilt should not be taken at face value here:

  Mainville, A. & Craymer, M.R. (2005) Present-day tilting of the Great Lakes
    region based on water level gauges. GSA Bull. 117, 1070-1080.
    55 gauges, 1860-2000 monthly means, relative vertical velocities good to
    +-1 cm/century within a lake (+-6 between lakes). Their headline negative
    result: ICE-4G contour gradients are "too small" over the Great Lakes;
    ICE-3G fits better. ICE-5G/6G/7G inherit ICE-4G's VM2-family viscosity, so
    the same gradient problem is expected to survive into ICE-7G -- and it
    does; see the report printed by __main__.

  Sella, G.F. et al. (2007) Observation of glacial isostatic adjustment in
    "stable" North America with GPS. GRL 34, L02306.
    362 GPS sites. ~10 mm/yr uplift at Hudson Bay, subsidence of 1-2 mm/yr
    south of the Great Lakes, hinge line running through the lakes. This is an
    absolute frame, which the gauges are not, so it is what fixes the datum.

  Brierley-Green, C. (2023) Glacial isostatic adjustment modelling for crustal
    motion in North America. MSc thesis, Univ. of Victoria.
    Grid search over mantle viscosity and lithospheric thickness against North
    American GNSS: best vertical fit at 100 km lithosphere, and even then a
    1.30 mm/yr vertical RMS residual remains against a spherically symmetric
    Earth. That residual is the point: a laterally homogeneous forward model
    cannot be tuned to fit this region exactly, so a smooth empirical
    correction fitted to the observations is the honest way to use them.

WHAT THIS MODULE DOES

  1. Stitches the four per-lake gauge solutions into one relative frame using
     the inter-lake ties Mainville & Craymer used for their Figure 7
     (Gros Cap = Thessalon; Toronto = Collingwood - 6; Buffalo = Port Weller).
     Reproduces their published Rossport-Calumet maximum of ~57 cm/century,
     which is the check that the stitch is right.
  2. Differences that against the base model's present-day uplift rate.
  3. Fits a smooth cubic trend surface to the residual, with a free datum
     offset per lake carrying the +-6 cm/century inter-lake tie uncertainty as
     a prior. Cubic is chosen by leave-one-out CV; quartic overfits.
  4. Propagates the present-day rate correction back through time as

         dU(x,t) = dv(x) . T . (1 - exp(-t/T))

     with T the late-stage relaxation time (~4 ka, measured from the base
     model's own rate decay -- see __main__).

WHY THAT FORM, AND WHAT IT DELIBERATELY DOES NOT DO

  dU(x,0) = 0, so modern topography is untouched -- it is data, not model.
  d(dU)/dt at t=0 = dv(x), so the corrected field reproduces the observed
  present-day tilt by construction. And dU saturates at dv.T rather than
  growing, which is the important choice: a present-day rate discrepancy is
  NOT evidence that the same fractional error held through deglaciation. The
  deglacial uplift is set by the ice load history and is independently
  constrained by raised strandlines; the gauges constrain only the slow tail.
  Propagating the residual back with the physically-correct-looking but
  divergent dv.T.(exp(t/T)-1) would put tens of metres of unconstrained
  correction at 14.5 ka and swamp exactly the signal the shorelines pin down.
  So the correction is bounded to the ~4 ka window the gauge data can actually
  speak for -- a couple of metres of differential tilt across the basins,
  which is the Nipissing/Algonquin outlet window and nothing older.

  Outside the gauge network the correction is tapered to zero. There is no
  data there, and the base model already matches Sella's ~10 mm/yr at Hudson
  Bay, so there is nothing to fix and every reason not to let a cubic
  extrapolate.
"""
import numpy as np, geopandas as gpd, csv

GAUGE_CSV = 'data/obs/mainville_craymer_2005_table3.csv'
GDB = {'ICE6G': ('data/raw/ice6g_pts/ICE6G_Data_points.gdb', 'ICE6G_datapoint_%05d'),
       'ICE7G': ('data/raw/ice7g/ICE7G_Data_points_all.gdb', 'ICE7G_datapoints_%05d')}

RATE_STEP = 500          # yr; present-day rate taken as GIA(500)/500
TAU       = 4000.0       # yr; late-stage relaxation time (measured, see __main__)
DEGREE    = 3            # trend-surface degree, chosen by LOO CV
TIE_SIGMA = 6.0          # cm/century, Mainville & Craymer inter-lake tie uncertainty
TAPER_IN  = 100.0        # km from nearest gauge: correction at full strength
TAPER_OUT = 350.0        # km from nearest gauge: correction is zero
LAKES     = ['superior', 'huronmich', 'erie', 'ontario']


def gauges():
    """The 55 gauges, stitched into one relative frame (Lakeport = 0)."""
    rows = [r for r in csv.DictReader(
        l for l in open(GAUGE_CSV) if not l.startswith('#'))]
    v = {r['gauge']: float(r['vel']) for r in rows}
    off = {'huronmich': 0.0, 'superior': v['Thessalon'] - v['Gros Cap']}
    off['ontario'] = (v['Collingwood'] - 6.0) - v['Toronto']
    off['erie'] = (v['Port Weller'] + off['ontario']) - v['Buffalo Harbor']
    keep = [r for r in rows if r['rejected'] == '0']
    return dict(
        name=[r['gauge'] for r in keep],
        lake=np.array([r['lake'] for r in keep]),
        lat=np.array([float(r['lat']) for r in keep]),
        lon=np.array([float(r['lon']) for r in keep]),
        # floor the reported sigma: the outlet gauges are 0 by construction
        sigma=np.array([max(float(r['sigma']), 0.3) for r in keep]),
        vel=np.array([float(r['vel']) + off[r['lake']] for r in keep]),
        offsets=off,
        span=v['Rossport'] + off['superior'] - v['Calumet Harbor'])


def _spline(model, t):
    from scipy.interpolate import RectBivariateSpline
    path, fmt = GDB[model]
    g = gpd.read_file(path, layer=fmt % t)
    lon = np.where(g['lon'].values > 180, g['lon'].values-360, g['lon'].values)
    lat = g['lat'].values
    m = (lon >= -105) & (lon <= -60) & (lat >= 32) & (lat <= 60)
    lo, la, gi = lon[m], lat[m], g['GIA'].values[m].astype('f8')
    ulo, ula = np.unique(lo), np.unique(la)
    grid = np.full((len(ula), len(ulo)), np.nan)
    grid[np.searchsorted(ula, la), np.searchsorted(ulo, lo)] = gi
    assert not np.isnan(grid).any(), 'gappy GIA lattice'
    return RectBivariateSpline(ula, ulo, grid, kx=3, ky=3, s=0)


def model_rate(model, lat, lon, t=RATE_STEP):
    """Base-model present-day uplift rate at scattered points, cm/century."""
    spl = _spline(model, t)
    r = np.array([spl(a, o)[0, 0]
                  for a, o in zip(np.atleast_1d(lat), np.atleast_1d(lon))])
    return r * (100.0 * 100.0 / t)      # m over t yr -> cm/century


# ----------------------------------------------------------------- trend fit
def _terms(deg, lo, la, lon0, lat0):
    """Polynomial basis, scaled by 5 deg so the normal equations stay tame."""
    X = [np.ones_like(lo)]
    for n in range(1, deg+1):
        for i in range(n+1):
            X.append(((lo-lon0)/5.0)**(n-i) * ((la-lat0)/5.0)**i)
    return np.array(X).T


def _design(deg, lo, la, lk, lon0, lat0):
    P = _terms(deg, lo, la, lon0, lat0)
    L = np.array([[1.0*(x == k) for k in LAKES] for x in lk])
    return np.hstack([P, L]), P.shape[1]


def _solve(deg, lo, la, lk, y, w, lon0, lat0):
    X, npoly = _design(deg, lo, la, lk, lon0, lat0)
    A = np.vstack([X*np.sqrt(w)[:, None],
                   np.hstack([np.zeros((4, npoly)), np.eye(4)/TIE_SIGMA])])
    b = np.r_[y*np.sqrt(w), np.zeros(4)]
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    return c, npoly


def fit(model=None, deg=DEGREE, base_rate=None, G=None):
    """Fit the residual trend surface. Returns everything needed to evaluate
    it, plus the diagnostics __main__ reports."""
    G = G if G is not None else gauges()
    rate = base_rate if base_rate is not None else model_rate(model, G['lat'], G['lon'])
    resid = G['vel'] - rate
    w = 1.0/G['sigma']**2
    lon0, lat0 = G['lon'].mean(), G['lat'].mean()
    c, npoly = _solve(deg, G['lon'], G['lat'], G['lake'], resid, w, lon0, lat0)
    X, _ = _design(deg, G['lon'], G['lat'], G['lake'], lon0, lat0)
    pred_full = X @ c                                   # incl. per-lake datum
    surf_at_gauge = _terms(deg, G['lon'], G['lat'], lon0, lat0) @ c[:npoly]
    # The gauges carry no absolute datum, so only the *pattern* is ours to
    # apply; demean it and leave the base model regional mean rate alone.
    # Sella Hudson Bay rate and hinge-line latitude both check out against the
    # base model, so that mean is the part worth keeping.
    mean = float(np.average(surf_at_gauge, weights=w))
    dv = surf_at_gauge - mean
    # The "before" to quote is the datum-only fit -- per-lake offsets and
    # nothing else. Residual about a single global mean would be dominated by
    # the arbitrary datum the gauges do not carry, and would flatter the
    # correction by counting a constant it never claimed to fix.
    if deg == 0:
        wrms_before = float(np.sqrt(np.average((resid-pred_full)**2, weights=w)))
    else:
        c0, _ = _solve(0, G['lon'], G['lat'], G['lake'], resid, w, lon0, lat0)
        X0, _ = _design(0, G['lon'], G['lat'], G['lake'], lon0, lat0)
        wrms_before = float(np.sqrt(np.average((resid-X0 @ c0)**2, weights=w)))
    return dict(model=model, deg=deg, coef=c[:npoly], lon0=lon0, lat0=lat0,
                mean=mean, gauge_lat=G['lat'], gauge_lon=G['lon'],
                lo=float(dv.min()), hi=float(dv.max()),
                resid=resid, pred=pred_full, rate=rate, w=w, G=G,
                wrms_before=wrms_before,
                wrms_after=float(np.sqrt(np.average((resid-pred_full)**2, weights=w))))


def surface(F, lon, lat):
    """Evaluate the fitted correction dv (cm/century) on a lon/lat mesh:
    demeaned, saturated to the range seen at the gauges, and tapered out.

    A cubic diverges outside the region it was fitted in, so two guards are
    needed rather than one. tanh saturation bounds the correction to the range
    the gauges actually demonstrate -- smoothly, since a hard clip would leave
    a kink in the uplift field and hence a ridge in the palaeo-surface -- and
    the distance taper then takes it to zero where there are no gauges at all.
    """
    lo = np.asarray(lon, dtype='f8'); la = np.asarray(lat, dtype='f8')
    v = _terms(F['deg'], lo.ravel(), la.ravel(), F['lon0'], F['lat0']) @ F['coef']
    v = v.reshape(lo.shape) - F['mean']
    mid, half = 0.5*(F['hi']+F['lo']), 0.5*(F['hi']-F['lo'])
    v = mid + half*np.tanh((v-mid)/half)
    return v * _taper(F, lo, la)


def _taper(F, lon, lat):
    """cos^2 rolloff on distance to the nearest gauge."""
    la0 = np.radians(0.5*(F['gauge_lat'].min()+F['gauge_lat'].max()))
    gx = F['gauge_lon']*111.32*np.cos(la0); gy = F['gauge_lat']*111.32
    x = lon*111.32*np.cos(la0); y = lat*111.32
    d = np.full(x.shape, np.inf)
    for a, b in zip(gx, gy):
        np.minimum(d, np.hypot(x-a, y-b), out=d)
    u = np.clip((d-TAPER_IN)/(TAPER_OUT-TAPER_IN), 0.0, 1.0)
    return np.cos(0.5*np.pi*u)**2


def rate_field(model, lon, lat, F=None):
    """The correction dv on a lon/lat mesh, in m/yr, ready to scale by
    time_factor(). Separating the space and time parts matters: the full dU
    cube is 29 x 1800 x 3000 floats, and it never has to exist."""
    F = F if F is not None else fit(model)
    return (surface(F, lon, lat) / 1e4).astype('f4')      # cm/century -> m/yr


def time_factor(t, tau=TAU):
    """tau.(1 - exp(-t/tau)): 0 at t=0, slope 1 there, saturating at tau."""
    return tau*(1.0 - np.exp(-t/tau))


def uplift_correction(F, lon, lat, steps, tau=TAU):
    """dU(x,t) in metres for each t in `steps`, on the given lon/lat mesh.

    dv . tau . (1 - exp(-t/tau)): exact present-day rate, zero at t=0, and
    bounded at dv.tau so the shoreline-constrained deglacial uplift is left
    alone. See the module docstring for why it is bounded and not divergent.
    """
    dv = surface(F, lon, lat) / 1e4              # cm/century -> m/yr
    return np.stack([dv*time_factor(t, tau) for t in steps]).astype('f4')


# --------------------------------------------------------------- diagnostics
def relaxation_time(model, sites=None, window=6000):
    """Measure T by fitting ln(rate) against t over the last `window` yr.
    Sites near the hinge are dropped: their rate is small and changes sign, so
    the log fit there is meaningless, not informative."""
    sites = sites or [('N Superior', 48.8, -87.5), ('Michipicoten', 47.96, -84.90),
                      ('North Bay', 46.32, -79.45), ('Parry Sound', 45.34, -80.04),
                      ('Lakeport', 43.15, -82.5), ('Chicago', 41.7, -87.5),
                      ('Cleveland', 41.5, -81.6)]
    ts = list(range(0, window+1, 500))
    spl = {t: (None if t == 0 else _spline(model, t)) for t in ts}
    U = lambda t, a, o: 0.0 if t == 0 else spl[t](a, o)[0, 0]
    out = []
    for nm, a, o in sites:
        mid = np.array([(x+y)/2 for x, y in zip(ts[:-1], ts[1:])], dtype='f8')
        r = np.array([(U(y, a, o)-U(x, a, o))/(y-x)*1000
                      for x, y in zip(ts[:-1], ts[1:])])
        if not (np.all(r > 0) or np.all(r < 0)):
            out.append((nm, r[0], None, None)); continue
        sl, ic = np.polyfit(mid, np.log(np.abs(r)), 1)
        r2 = 1 - np.var(np.log(np.abs(r))-(ic+sl*mid))/np.var(np.log(np.abs(r)))
        out.append((nm, r[0], 1/sl, r2))
    return out


def _loo(deg, G, resid, w, lon0, lat0):
    e = np.empty(len(resid))
    for i in range(len(resid)):
        k = np.r_[np.arange(i), np.arange(i+1, len(resid))]
        c, _ = _solve(deg, G['lon'][k], G['lat'][k], G['lake'][k],
                      resid[k], w[k], lon0, lat0)
        X, _ = _design(deg, G['lon'][i:i+1], G['lat'][i:i+1],
                       G['lake'][i:i+1], lon0, lat0)
        e[i] = resid[i] - (X @ c)[0]
    return float(np.sqrt(np.average(e**2, weights=w)))


def report(model='ICE7G'):
    G = gauges()
    print(f'=== GIA calibration report: {model} ===\n')
    print('-- gauge stitch (Mainville & Craymer 2005, Fig. 7 ties) --')
    for k, v in G['offsets'].items():
        print(f'   {k:10s} datum offset {v:+6.1f} cm/century')
    print(f'   Rossport-Calumet span {G["span"]:.1f} cm/century '
          f'(published: ~57)  {"OK" if abs(G["span"]-57) < 1 else "MISMATCH"}')
    print(f'   {len(G["name"])} gauges used, 3 rejected by the authors\n')

    print('-- Sella et al. 2007 absolute checks on the base model --')
    for nm, a, o in [('Hudson Bay (58N,86W)', 58.0, -86.0),
                     ('N Lake Superior', 48.8, -87.5),
                     ('S Lake Michigan', 41.7, -87.5),
                     ('S of lakes (40N,84W)', 40.0, -84.0)]:
        print(f'   {nm:22s} {model} {model_rate(model, a, o)[0]/10:+6.2f} mm/yr')
    print('   GPS: ~+10 Hudson Bay, -1 to -2 south of the lakes, hinge through'
          ' the lakes\n')

    print('-- late-stage relaxation time from the base model rate decay --')
    Ts = []
    for nm, r0, T, r2 in relaxation_time(model):
        if T is None or T < 0:
            print(f'   {nm:14s} rate {r0:+6.2f} mm/yr   (near hinge, skipped)')
        else:
            print(f'   {nm:14s} rate {r0:+6.2f} mm/yr   T = {T:5.0f} yr  r2={r2:.3f}')
            Ts.append(T)
    print(f'   median T = {np.median(Ts):.0f} yr;  module uses TAU = {TAU:.0f} yr\n')

    print('-- per-lake tilt gradient, observed vs base model --')
    rate = model_rate(model, G['lat'], G['lon'])
    for lk in LAKES:
        s = G['lake'] == lk
        w = 1/G['sigma'][s]**2
        o = G['vel'][s]-np.average(G['vel'][s], weights=w)
        m = rate[s]-np.average(rate[s], weights=w)
        k = np.sum(w*o*m)/np.sum(w*m*m)
        print(f'   {lk:10s} n={s.sum():2d}  observed range {o.max()-o.min():5.1f}'
              f'  model range {m.max()-m.min():5.1f}  obs/model {k:5.2f}')
    print('   obs/model < 1: model tilts too steeply.  > 1: too flat.\n')

    print('-- trend-surface degree, leave-one-out CV (cm/century) --')
    resid = G['vel']-rate; w = 1/G['sigma']**2
    lon0, lat0 = G['lon'].mean(), G['lat'].mean()
    for deg in (0, 1, 2, 3, 4):
        F = fit(model, deg=deg, base_rate=rate, G=G)
        tag = '  <- datum only, no pattern' if deg == 0 else \
              ('  <- chosen' if deg == DEGREE else '')
        print(f'   deg {deg}  fit wRMS {F["wrms_after"]:5.2f}   '
              f'LOO wRMS {_loo(deg, G, resid, w, lon0, lat0):5.2f}{tag}')

    F = fit(model, base_rate=rate, G=G)
    print(f'\n-- fitted correction (degree {DEGREE}) --')
    print(f'   dv over the gauge network: {F["lo"]:+.1f} to {F["hi"]:+.1f} cm/century')
    print(f'   bounded dU at t >> TAU:    {F["lo"]*TAU/1e4:+.2f} to '
          f'{F["hi"]*TAU/1e4:+.2f} m')
    print(f'   residual wRMS {F["wrms_before"]:.2f} (datum only) -> '
          f'{F["wrms_after"]:.2f} cm/century (gauge precision ~1)')
    worst = np.argsort(-np.abs(F['resid']-F['pred']))[:5]
    print('   largest remaining misfits:')
    for i in worst:
        print(f'     {G["name"][i]:16s} {G["lake"][i]:10s} '
              f'obs {G["vel"][i]:+6.1f}  model {F["rate"][i]:+6.1f}  '
              f'resid after fit {F["resid"][i]-F["pred"][i]:+5.1f}')


if __name__ == '__main__':
    import sys
    report(sys.argv[1] if len(sys.argv) > 1 else 'ICE7G')
