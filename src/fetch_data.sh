#!/usr/bin/env bash
# Download every input dataset. About 700 MB. Run from the repository root.
set -euo pipefail
mkdir -p data/raw && cd data/raw

echo "== ETOPO 2022 15 arc-second, NOAA NCEI"
B=https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/15s/15s_surface_elev_netcdf
for T in N45W105 N45W090 N45W075 N60W105 N60W090 N60W075; do
  curl -L -o "etopo_$T.nc" "$B/ETOPO_2022_v1_15s_${T}_surface.nc"
done

echo "== NADI-1 ice margin isochrones, Dalton et al. 2023"
curl -L -o nadi1.zip \
  "https://zenodo.org/records/8161764/files/NADI-1%20shapefiles%20Dalton%20et%20al.%20QSR.zip?download=1"
unzip -q -o nadi1.zip -d nadi1

echo "== ICE-7G_NA and ICE-6G_C deformation points, PANGAEA 947536"
# Peltier's own server and the PMIP4 mirror both time out; this is the repackaging
curl -L -o ice7g.zip "https://download.pangaea.de/dataset/947536/files/ICE7G_Data_points_all.zip"
curl -L -o ice6g.zip "https://download.pangaea.de/dataset/947536/files/ICE6G_Data_points.zip"
unzip -q -o ice7g.zip -d ice7g && unzip -q -o ice6g.zip -d ice6g_pts

echo "== GSHHG 2.3.7 full resolution shorelines (LGPL)"
curl -L -o gshhg.zip "https://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip"
unzip -q -o gshhg.zip "GSHHS_shp/f/*" -d gshhg

echo "== Natural Earth 1:10m lakes (public domain)"
curl -L -o ne_lakes.zip "https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_lakes.zip"
unzip -q -o ne_lakes.zip -d ne_lakes

echo "== NOAA Environmental Sensitivity Index atlases"
B=https://response.restoration.noaa.gov/sites/default/files/esimaps/gisdata
mkdir -p esi && cd esi
for F in Lake_Michigan_ESI_2025_GDB Lake_Ontario_ESI_2023_GDB LakeHuron_1994_GDB \
         GL_Straits_of_Mackinac_2019_GDB GL_StClair_Detroit_River_System_2019_GDB \
         Great_Lakes_St_Marys_River_2021_ESI_GDB Great_Lakes_St_Lawrence_River_2021_ESI_GDB \
         LakeSuperior_1994_Source; do
  curl -L -o "$F.zip" "$B/$F.zip" && unzip -q -o "$F.zip" -d "$F"
done
curl -L -o Lake_Erie_ESI_2022_GDB.zip "$B/Lake_Erie_ESI_2022_GDB.zip"
unzip -q -o Lake_Erie_ESI_2022_GDB.zip -d ../esi_erie
cd ..

echo "== ECCC datasets"
# Both ECCC portals are JavaScript catalogues with no plain directory listing.
# Their file API is /api/file?path=<url-encoded absolute path>.
python3 - <<'PY'
import urllib.parse, urllib.request, os
def get(host, path, out):
    u = f"https://{host}/api/file?path={urllib.parse.quote('/'+path.lstrip('/'), safe='')}"
    urllib.request.urlretrieve(u, out)
    print(f"  {out}  {os.path.getsize(out)/1e6:.1f} MB")

R = "sites/scientificknowledge/great-lakes-nearshore-waters-assessment"
os.makedirs("eccc", exist_ok=True)
for k, p in {
  "erie":     f"{R}/lake-erie-nearshore-waters-assessment/en/LakeErieNearshoreWatersAssessment_z17N.gdb.zip",
  "ontario":  f"{R}/lake-ontario-nearshore-waters-assessment/en/LakeOntarioNearshoreWatersAssessment_z17N.gdb.zip",
  "huron":    f"{R}/lake-huron-nearshore-waters-assessment/English/LakeHuronNearshoreWatersAssessment_z17N.gdb.zip",
  "superior": f"{R}/lake-superior-nearshore-waters-assessment/en/LakeSuperiorNearshoreWatersAssessment_z16N.gdb.zip",
}.items():
    get("data-donnees.az.ec.gc.ca", p, f"eccc/{k}.gdb.zip")

S = ("sites/emergencies/shoreline-segmentation-with-shoreline-cleanup-assessment-"
     "technique-scat-classification/ontario-shoreline-classification/"
     "ShorelineClassification_ON_OpenDataCatalogue.gdb.zip")
os.makedirs("scat", exist_ok=True)
get("data-donnees.ec.gc.ca", S, "scat/ShorelineClassification_ON_OpenDataCatalogue.gdb.zip")
PY
for f in eccc/*.gdb.zip; do unzip -q -o "$f" -d "eccc/$(basename "${f%.gdb.zip}")"; done
unzip -q -o scat/*.gdb.zip -d scat/

echo
echo "Done. The 1994 atlases also need:  apt-get install mdbtools e00compr"
