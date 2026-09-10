# Great Lakes History

Palaeogeography of the Great Lakes from 14.5 ka to the present, solved rather than
drawn, plus a modern shoreline dataset carrying substrate and drift-class attributes.

Every shoreline in the reconstruction is the output of flooding a real
topo-bathymetric grid. Nothing is traced from published maps.

## Method

For each 500-year timestep:

1. **Palaeo surface** = modern topo-bathymetry minus the isostatic uplift that has
   occurred since that date.
2. **Ice** from the NADI-1 isochrone for that date, rasterised as an impassable wall.
3. **Depression fill** by grayscale reconstruction with the domain edge as the only
   drain. Water depth is the filled surface minus the ground, so each lake's level is
   the elevation of the lowest sill on its own rim — found, not assumed.
4. **Outlet** located where the lake touches ground that is below its level and
   connected to the domain edge, walking back up any flat spillway reach.
5. **Marine flooding** handled separately from the eustatic curve read off the deep
   ocean cells of the same rebound file, then trimmed by the lakes: a basin under a
   water plane above sea level is fresh whatever its floor elevation.

Domain 95°W–70°W, 38°N–53°N at 30 arc-seconds, 5.4 million cells.

## Control run

Nothing is tuned. At 0 ka the lakes emerge from flooding the modern grid:

| | modelled | actual |
|---|---|---|
| Michigan–Huron | 119,151 km², 174.4 m | 117,400 km², 176.0 m |
| Superior | 83,026 km², 180.5 m | 82,100 km², 183.5 m |
| Erie | 26,233 km², 170.1 m | 25,700 km², 174.0 m |
| Ontario | 19,829 km², 73.5 m | 18,960 km², 74.2 m |

All four modern outlets are located by the model: Port Huron (−82.42, 43.00), the
Sault (−84.46, 46.48), Niagara (−78.91, 42.92), Thousand Islands (−75.80, 44.50).

Levels run 1–4 m low because 30″ cells smooth the outlet channels, and because the
DEM had to be aggregated by minimum rather than mean — averaging closed the St. Clair,
Niagara and St. Marys channels outright, which are narrower than one cell.

## What it reproduces, and what it does not

Reproduces:

- A proglacial lake ponded across the Erie, Saginaw and Huron basins at ~218 m at
  14 ka, close to the mapped Lake Warren shoreline.
- The Champlain Sea, unprompted, in the St. Lawrence and Ottawa lowlands.
- The Nipissing highstand: confluent upper lakes at **183.8 m spilling at Port Huron**
  at 5 ka, against a mapped strandline of about 184 m.
- Superior separating from Michigan–Huron at the Sault between 2.5 and 2 ka
  (between 2 and 1.5 ka before the gauge calibration).

Does not:

- **Main Lake Algonquin.** Under the maximum ice margin the model gets a 168,000 km²
  confluent upper-lakes lake at 166 m draining through Kirkfield at 11.5 ka — right
  lake, right outlet, 18 m low. Main Algonquin at 184 m on Port Huron never forms,
  and the reason is structural: it existed because the Port Huron readvance *re-closed*
  Kirkfield after it had opened, and NADI-1 shrinks monotonically at every step in all
  three of its variants. No choice of margin envelope can produce it.
- **Erie's Holocene history**, which comes out backwards — falling 200→170 m instead of
  rising 145→174 m. That is an overdriven peripheral forebulge, not a routing error.
- **Sill geometry of the past.** The model floods today's channels. Niagara has cut
  ~11 km of gorge since 12.4 ka; Chicago and St. Clair are dredged.
- **Small, shallow basins below the resolution of the rebound field.** Lake Simcoe
  solves as intermittently dry for several thousand years around the mid-Holocene,
  which the archaeological record does not support — the basin has almost certainly
  held water continuously since deglaciation. ICE-7G/ICE-6G are 1° global lattices;
  a bedrock sill a few kilometres across, with only a few tens of metres of relief, is
  well below that. `model.py` now smooths the interpolated uplift field across
  timesteps at each grid cell (`_gia_cube`, a light Gaussian filter in time) on the
  reasoning that real postglacial uplift is smooth in time, so a single-timestep spline
  departure from its neighbours is fit noise, not signal. That change was verified
  against the control run and the Nipissing highstand — both unchanged to within
  rounding — but it did **not** rescue Simcoe: the same 12–3 ka dry gap comes out
  whether or not the smoothing is applied. That null result is itself informative — it
  means the problem there is not per-timestep interpolation jitter but a sustained,
  smooth-in-time misestimate of the basin's sill-to-floor differential, i.e. the 1°
  lattice genuinely cannot resolve a basin this small. No amount of temporal smoothing
  fixes a spatial resolution limit, and pushing the filter harder to force Simcoe right
  would mean smearing the real multi-millennial trend (Nipissing, Champlain Sea timing)
  to fix one basin by hand, which is the tuning this project specifically avoids.

Switching the rebound model from ICE-6G_C to ICE-7G_NA fixed the Nipissing highstand
and left the Algonquin failure untouched, which is what localises that failure to the
ice-margin dataset rather than to the rebound field.

## Calibrating the rebound field

ICE-6G and ICE-7G are global inversions, and neither was ever fitted to Great Lakes
lake-level gauges. Over these basins the modelled tilt is an extrapolation, and it is
measurably wrong. `src/gia_calib.py` corrects it against three observational datasets
and is on by default; `CALIB=0` restores the raw base model.

Mainville & Craymer (2005) give relative vertical velocities at 55 gauges from 140
years of monthly water levels, good to ±1 cm/century within a lake. Stitched into one
frame with their own inter-lake ties — which reproduces their published 57 cm/century
Rossport–Calumet maximum, the check that the stitch is right — they show ICE-7G tilting
Michigan–Huron about 30% too steeply and Erie less than half steeply enough:

| basin | gauges | observed range | ICE-7G range | observed / model |
|---|---:|---:|---:|---:|
| Superior | 10 | 52.8 | 46.5 | 1.06 |
| Michigan–Huron | 20 | 41.4 | 32.4 | 0.77 |
| Erie | 13 | 14.2 | 4.8 | 2.27 |
| Ontario | 9 | 22.5 | 12.0 | 1.28 |

Ranges in cm/century about each basin's own datum. This is the same complaint Mainville
& Craymer made against ICE-4G, whose VM2-family viscosity ICE-5G/6G/7G all inherit.

The correction is a cubic trend surface (chosen by leave-one-out CV; quartic overfits)
fitted to the residual with a free datum offset per lake carrying their ±6 cm/century
inter-lake tie uncertainty as a prior. It brings the weighted residual from 4.24 to
1.80 cm/century against a gauge precision near 1.

The gauges carry no absolute datum, so only the pattern is applied and the base model's
regional mean rate is kept. Sella et al. (2007) is what justifies keeping it: 362 GPS
sites put Hudson Bay uplift near 10 mm/yr where ICE-7G gives 10.5, and the hinge line
through the lakes where ICE-7G puts it at 43–44°N. Their one disagreement is the far
field, where GPS shows 1–2 mm/yr of subsidence south of the lakes against ICE-7G's
0.3–1.0. The calibration passes that as an independent test it was not fitted to:
modelled subsidence at Chicago moves from −0.4 to −1.5 mm/yr, into the observed band.

Propagating the present-day rate correction back in time uses
`ΔU(t) = Δv · τ · (1 − exp(−t/τ))`, with τ ≈ 4 ka measured from the decay of ICE-7G's
own uplift rate. It is zero at the present, so modern topography is untouched, and it
**saturates** rather than growing. The single-mode alternative `Δv · τ · (exp(t/τ) − 1)`
looks more physical and is wrong here: it would put tens of metres of unconstrained
correction at 14.5 ka and swamp the deglacial uplift, which the raised strandlines
already constrain. The gauges measure the slow tail and are bounded to it. Outside the
gauge network the correction tapers to zero.

What it changed: Lake Superior separates from Michigan–Huron at the Sault between
2.5 and 2 ka rather than between 2 and 1.5 ka — one timestep earlier, and both inside
the range the late-Holocene record allows, so a shift in timing rather than a fix; the
Michigan basin sits 6–7 m lower against the North Bay sill at 11 ka. Everywhere else
the largest lake's level moves 0.2–2.9 m and its spill point does not move at all.
Brierley-Green (2023) is the reason to treat this as an empirical patch rather than a
better Earth model: the best fit over mantle viscosity and lithospheric thickness
against North American GNSS still leaves 1.30 mm/yr of vertical residual, so lateral
heterogeneity is doing work no spherically symmetric model can absorb.

`python src/gia_calib.py [ICE7G|ICE6G]` prints the full validation report — the stitch
check, the Sella comparisons, the measured relaxation times, the per-basin gradients,
the CV table and the largest remaining misfits.

## Layout

    index.html                     the deployed page: palaeolakes, ice, Algonquin,
                                   uplift/topo-bathymetry/model-diff rasters and the
                                   rebound transect as togglable layers on one map,
                                   fetching data/*.json at runtime
    assets/
      style.css, app.js            styling and logic for index.html
    viewers/
      great-lakes-modelled.html    standalone single-file snapshot of the palaeolake
                                   view, all data embedded; kept for offline sharing
      great-lakes-rebound.html     standalone single-file snapshot of the rebound view,
                                   all data embedded; kept for offline sharing
    data/
      paleo_optimal.json           solved lakes, 29 timesteps, NADI-1 optimal margin
      paleo_max.json               same, maximum margin
      algonquin.json               the mapped 184 m Algonquin plane on the model's
                                   palaeo-topography, 13–11.5 ka
      rebound.json                 uplift lattices, background grid, transect
      raw/                         fetched inputs, not committed (fetch_data.sh)
      gia_cube_ICE6G.npz,
      gia_cube_ICE7G.npz           cached, time-smoothed uplift lattices; not
                                   committed, rebuilt on first run of model.py
    gis/
      great_lakes.gpkg             all vector layers in one file (preferred)
      great_lakes_shoreline_typed.shp   500 m segments with substrate and drift class
      great_lakes_gshhg_full.shp   lake polygons, GSHHG, best angles
      great_lakes_etopo15s.shp     lake polygons, flooded from ETOPO, best areas
      great_lakes_naturalearth.shp lake polygons, generalised, public domain
      DATA_README.txt              full field descriptions and provenance
    src/
      build_dem.py                 mosaic the ETOPO tiles
      model.py                     the flooding model
      gia_calib.py                 gauge/GPS calibration of the rebound field
      export.py                    run all timesteps, write GeoJSON
      export_rebound.py            rebound lattices and transect
      algonquin.py                 the imposed Algonquin water plane
      make_shapefiles.py           modern lake polygons, three sources
      join_shoretype4.py           shoreline type join and fused drift class
      fetch_data.sh                download every input dataset

## Shoreline dataset

`gis/great_lakes_shoreline_typed` cuts the GSHHG outlines into 500 m segments and tags
each with substrate from NOAA ESI (US) and ECCC SCAT (Canada), plus a fused
`drift_cls` built for one question: can this shore take part in alongshore drift?

| | bedrock | sediment | armoured | low energy | unknown |
|---|---|---|---|---|---|
| Erie | 83 | 395 | 583 | 127 | 34 |
| Huron | 549 | 896 | 267 | 555 | 102 |
| Michigan | 14 | 818 | 465 | 521 | 5 |
| Ontario | 127 | 377 | 438 | 130 | 198 |
| Superior | 420 | 526 | 67 | 50 | 1,243 |

km per class. The remaining unknown is almost all the US shore of Lake Superior, whose
only ESI atlas is 1994 and sits in an undefined local oblique mercator that cannot be
georeferenced without inventing a fit.

GSHHG is used as the base geometry because raster boundaries quantise orientation:
44.5% of ETOPO-derived segment bearings fall within 5° of a 45° multiple, against
22.1% for GSHHG, and 22% is what an unbiased distribution gives.

## Sources

| layer | source |
|---|---|
| Topography and bathymetry | ETOPO 2022 15″, NOAA NCEI, carrying the NOAA Great Lakes bathymetric grids |
| Ice margins | NADI-1, Dalton et al. 2023, *Quaternary Science Reviews* 321, 108345 (Zenodo 8161764) |
| Isostatic rebound | ICE-7G_NA (VM7), Roy & Peltier 2017–2018, via Godbout, Brouard & Roy, PANGAEA 947536 |
| Present-day tilt (calibration) | Mainville & Craymer 2005, *GSA Bulletin* 117, 1070–1080 — 55 Great Lakes water-level gauges |
| Present-day uplift (calibration) | Sella et al. 2007, *Geophys. Res. Lett.* 34, L02306 — 362 North American GPS sites |
| Earth rheology (calibration) | Brierley-Green 2023, MSc thesis, University of Victoria — GNSS grid search over viscosity and lithospheric thickness |
| Shoreline geometry | GSHHG 2.3.7 (Wessel & Smith); Natural Earth 1:10m |
| Shoreline type, US | NOAA Environmental Sensitivity Index atlases, 1994–2025 |
| Shoreline type, Canada | ECCC Shoreline Classification (SCAT); ECCC Great Lakes Nearshore Waters Assessment |

`src/fetch_data.sh` downloads all of them. Raw inputs are not committed — they run to
roughly 700 MB.

## Reproducing

    bash src/fetch_data.sh          # ~700 MB of downloads
    python src/build_dem.py
    python src/gia_calib.py         # optional: print the calibration report
    MARGIN=OPTIMAL python src/export.py
    MARGIN=MAX     python src/export.py
    python src/export_rebound.py
    python src/algonquin.py
    python src/make_shapefiles.py
    python src/join_shoretype4.py

Needs numpy, scipy, scikit-image, rasterio, geopandas, shapely, pyogrio, netCDF4, and
`mdbtools` plus `e00compr` for the two 1994-vintage ESI atlases. Peak memory is about
1.5 GB at 30″; at 15″ it needs roughly 4 GB and gives sharper outlet channels.

`build_dem.py` mosaics the ETOPO tiles to `data/dem15.npy` and also writes the 30″
grid `model.py` actually reads, `data/dem30_min.npy`, by 2×2 block-minimum decimation
(minimum, not mean, for the reason given in the control run above). The first call to
`model.py` per rebound model also builds and caches `data/gia_cube_{ICE6G,ICE7G}.npz`
— every solved timestep's uplift lattice, smoothed across time — which both margin
runs of `export.py` then share; expect that one-time build to take a few minutes.

## Deploying

`index.html` at the repository root is a static page with no build step: it fetches
`data/*.json` at load time, so it only needs to be served over HTTP, not opened as a
`file://` URL. On GitHub Pages: repository **Settings → Pages → Source → Deploy from a
branch**, branch `main`, folder `/ (root)`. `.nojekyll` is committed so Pages serves the
`data/` and `assets/` directories as-is. To preview locally: `python -m http.server` from
the repository root, then open `http://localhost:8000/`.

## Caveats worth carrying forward

- ICE-7G is a 1° global model. The lattice visible under the rebound fields is the real
  resolution of the data; a bedrock sill a few kilometres across is far below it — see
  Lake Simcoe under "What it reproduces, and what it does not" for the concrete case
  and what was tried against it. The gauge calibration corrects the tilt over these
  basins; it does not add resolution, and it cannot.
- The calibration is fitted to *modern* rates and bounded to the last few thousand
  years. It should not be read as a claim about deglacial-age uplift, which it is
  deliberately built not to touch, and its 1.80 cm/century residual is still above the
  ~1 cm/century the gauges resolve.
- ESI atlas vintages range from 1994 to 2025. Lake Huron's US shore is typed from
  30-year-old mapping.
- ECCC nearshore reaches average 17 km on Superior against 1.3 km on Erie, so the
  hardening attribute is coarse in the north.
- GSHHG lake polygons derive from CIA WDBII and run 3–8% large on area; use the ETOPO
  polygons when area matters and GSHHG when orientation matters.
