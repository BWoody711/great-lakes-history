Modern Great Lakes polygons
===========================
Six features each: Superior, Michigan, Huron (incl. Georgian Bay), St. Clair,
Erie, Ontario.  EPSG:4326.  Attributes: name, datum_m (chart datum, IGLD 1985),
area_km2 (computed on EPSG:6933, not from degrees), nvert, source.

Three independent versions, because none is best at everything:

  great_lakes_etopo15s.shp
      Flooded from the NOAA NCEI ETOPO 2022 15" grid at each lake's chart
      datum, then vectorised.  Cell size ~450 m.  Most vertices and the
      closest areas; consistent with the palaeogeography reconstruction,
      which uses the same grid.  St. Clair runs large because the
      surrounding flats flood at this cell size.

  great_lakes_naturalearth.shp
      Natural Earth 1:10m physical lakes.  Public domain, no attribution
      required.  Areas are excellent, vertex density is low (~3 km spacing),
      so it is a good basemap layer and a poor analysis layer.

  great_lakes_gshhg_full.shp
      GSHHG 2.3.7 full resolution (Wessel & Smith).  Included because it is
      the usual high-resolution answer, but note its lake polygons derive
      from CIA WDBII rather than World Vector Shorelines, and here they run
      3-8% large: Superior 85,195 km2 and Huron 64,460 km2 against published
      82,100 and 59,600.  Licensed LGPL, unlike the other two.

Area check against published values (km2):

  lake        published     etopo15s   naturalearth   gshhg_full
  Superior       82,100       82,221         82,310       85,195
  Michigan       57,800       57,868         57,540       58,539
  Huron          59,600       59,447         59,937       64,460
  St. Clair       1,114        1,424          1,148        1,223
  Erie           25,700       26,052         25,736       25,713
  Ontario        18,960       19,420         19,465       19,716

All three sources join Superior, Michigan and Huron into single polygons,
since the basins are connected at the Straits of Mackinac and the St. Marys
River.  They have been split here with narrow cuts across the Straits of
Mackinac, the St. Marys, Detroit, St. Clair and Niagara rivers, so each lake
is one feature.  If you want the hydrologically confluent Michigan-Huron as
a single body, dissolve those two.

Built by make_shapefiles.py.


Shoreline type
==============
great_lakes_shoreline_typed.shp
    GSHHG lake outlines cut into 500 m segments (18,453 segments, 9,226 km),
    tagged from three sources.  EPSG:4326.

    Base geometry is GSHHG, not the ETOPO raster boundaries.  A raster
    boundary stair-steps at the cell size, so segment orientation is quantised
    to multiples of 45 degrees.  Measured on this data: 44.5% of ETOPO-derived
    bearings fall within 5 degrees of a 45 degree multiple, against 22.1% for
    GSHHG - and 22% is what an unbiased distribution gives.  Peak-to-median in
    a 36-bin bearing histogram is 5.4 for ETOPO, 1.7 for GSHHG.  Use this one
    for shore-normals, fetch or drift direction.
    The cost: GSHHG is generalised.  9,226 km here against 20,335 km for the
    ETOPO version - partly raster length inflation, mostly GSHHG dropping
    small islands and the fine structure of Georgian Bay and the North
    Channel.  Both base layers are in this folder.

    Fields:
      lake, bearing (0-180, 0 = north), seg_km
      drift_cls  THE FUSED FIELD.  One general class, built for a single
                 question: can this shore take part in alongshore drift?
      drift_src  which dataset that class came from
      mobile     Y if drift_cls is sediment, N otherwise, null if unknown
      substrate  physical class, merged from the three sources below
      sub_src    which source it came from
      esi        NOAA ESI rank (3B, 6A, 8C, ...)
      shoretype  full ESI description
      exposure   ESI ENVIR code (E / S / L)
      atlas,year ESI atlas and its year
      scat_cls   ECCC SCAT class, verbatim
      eccc_hard  ECCC nearshore assessment: Hardened / Natural
      eccc_unit  ECCC wave-energy class for the reach
      eccc_year  ECCC assessment year
      esi_m, scat_m, eccc_m   distance to the source feature each attribute
                 came from (caps 2000 / 2000 / 3000 m)

drift_cls
---------
Four classes plus null.  Assigned from the source codes directly, not from
the collapsed "substrate" field, because two ESI ranks collapse ambiguously
there (8F sheltered vegetated bluff and 9B vegetated low bank both land in
"veg. bank or bluff" but behave differently).

  bedrock     Rock shore.  Effectively no mobile sediment, and a fixed
              boundary that terminates or compartmentalises a drift cell.
              ESI 1A, 1C, 2, 2A, 8A;  SCAT Bedrock Cliff/Vertical,
              Bedrock Platform.
  sediment    Unconsolidated clastic shore, the drift-capable class: beaches,
              banks, bars and eroding bluffs, whether they are transporting
              or supplying.  ESI 3A, 3B, 4, 5, 6A, 8F;  SCAT Sand / Pebble-
              Cobble / Boulder / Mixed Sediment Beach or Bank, Sediment Cliff.
  armoured    Man-made hard shore.  Drift is interrupted or the shore is
              fixed; treat as a barrier, not as a source.  ESI 1B, 6B, 8B,
              8C;  SCAT Man-Made Solid, Man-Made Permeable;  and, as a last
              resort, ECCC nearshore ShType = Hardened.
  low energy  Marsh, swamp, flats, vegetated low bank.  Sheltered, negligible
              wave-driven transport.  ESI 7, 9A, 9B, 10A-10D;  SCAT Marsh,
              Mud Tidal Flat.
  null        No source covers the segment, or SCAT records Not Classified,
              or ECCC records only Natural, which does not distinguish sand
              from bedrock.

Totals: bedrock 1,194 km (13%), sediment 3,014 km (33%), armoured 1,961 km
(21%), low energy 1,473 km (16%), unknown 1,582 km (17%).

Two judgement calls worth knowing before you use it.  Boulder Beach or Bank
is in "sediment" because it is clastic, but boulders are largely immobile in
this wave climate - filter on scat_cls if that matters to you.  And ESI 8A,
"Sheltered Scarps (Bedrock/Mud/Clay)", is put in bedrock although the code
itself is ambiguous about the substrate; it is only 29 km on Erie, but check
it if it falls in your study area.

Sources, in the priority order used to fill "substrate" and "drift_cls":

  1. NOAA ESI, US shore, 4,292 km.  Substrate and form.  Atlases, newest wins:
     Michigan 2025, Ontario 2023, Erie 2022, St Marys R. 2021, St Lawrence R.
     2021, Straits of Mackinac 2019, St Clair/Detroit 2019, Huron 1994.

  2. ECCC Shoreline Classification (SCAT), Canadian shore, 3,250 km.  Part of
     the Environmental Emergencies shoreline database, 12,343 segments
     province-wide, of which 7,339 carry a class.  Classes: Bedrock
     Cliff/Vertical, Bedrock Platform, Pebble/Cobble Beach or Bank, Boulder
     Beach or Bank, Mixed Sediment Beach or Bank, Sand Beach or Bank, Sediment
     Cliff, Marsh, Mud Tidal Flat, Man-Made Solid, Man-Made Permeable, Not
     Classified.  Mapped into the same substrate vocabulary; the verbatim
     class is kept in scat_cls.  data-donnees.ec.gc.ca, path
     /sites/emergencies/shoreline-segmentation-with-shoreline-cleanup-
     assessment-technique-scat-classification/ontario-shoreline-classification

  3. ECCC Great Lakes Nearshore Waters Assessment, 101 km.  Only fills gaps
     where it reports Hardened, since Natural does not distinguish sand from
     bedrock.  Its real value here is eccc_hard and eccc_unit as independent
     attributes: Erie 2018 (1,075 km, 51% hardened), Ontario 2019 (1,573 km,
     51%), Huron 2021 (4,091 km, 49%), Superior 2020 (1,628 km, 34%).  Note
     its reaches are long - Superior averages 17 km per reach against 1.3 km
     for Erie.

Both ECCC portals are JavaScript catalogues with no plain directory listing.
Their APIs are /api/path_contents?path= and /api/file?path= on
data-donnees.ec.gc.ca and data-donnees.az.ec.gc.ca.  See join_shoretype3.py.

Coverage: substrate on 7,643 km (83%); any type information on 8,016 km (87%).
The 1,154 km with nothing is the US shore of Lake Superior, whose only ESI
atlas is 1994 and sits in an undefined local oblique mercator that cannot be
georeferenced without inventing a fit.

great_lakes_drift_class_summary.csv      km by lake and fused drift class
great_lakes_shoreline_type_summary.csv   km by lake and substrate
lake_erie_esi_shoretype_2022.shp         raw Lake Erie ESIL, unjoined
