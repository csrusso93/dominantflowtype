# -*- coding: utf-8 -*-
"""Optional reference-DEM search / download via the OpenTopography REST API."""
from __future__ import annotations


def opentopo_search(bbox_wgs84, api_key):
    """List OpenTopography point-cloud datasets intersecting bbox (best-effort)."""
    import requests
    s, w, n, e = bbox_wgs84
    try:
        url = "https://portal.opentopography.org/API/otCatalog"
        params = dict(productFormat="PointCloud", minx=w, miny=s, maxx=e,
                      maxy=n, detail="false", outputFormat="json",
                      include_federated="true")
        r = requests.get(url, params=params, timeout=60)
        cat = r.json()
        sets = cat.get("Datasets", [])
        print(f"  [OpenTopography] {len(sets)} dataset(s) intersect AOI:")
        for d in sets:
            ds = d.get("Dataset", {})
            print(f"     - {ds.get('name')}  ({ds.get('alternateName')}) "
                  f"[{ds.get('temporalCoverage')}]")
        return sets
    except Exception as ex:
        print(f"  [OpenTopography] catalog search failed: {ex}")
        return []


def opentopo_download(bbox_wgs84, api_key, demtype, out_path):
    """Download a global/US DEM for the AOI via OpenTopography REST API."""
    import requests
    s, w, n, e = bbox_wgs84
    try:
        if demtype.upper().startswith("USGS"):
            url = "https://portal.opentopography.org/API/usgsdem"
            params = dict(datasetName=demtype, south=s, north=n, west=w,
                          east=e, outputFormat="GTiff", API_Key=api_key)
        else:
            url = "https://portal.opentopography.org/API/globaldem"
            params = dict(demtype=demtype, south=s, north=n, west=w, east=e,
                          outputFormat="GTiff", API_Key=api_key)
        r = requests.get(url, params=params, timeout=180)
        if r.status_code == 200 and r.content[:2] in (b"II", b"MM"):
            with open(out_path, "wb") as fh:
                fh.write(r.content)
            print(f"  [OpenTopography] saved {demtype} DEM -> {out_path}")
            return out_path
        print(f"  [OpenTopography] download failed "
              f"(HTTP {r.status_code}): {r.text[:200]}")
    except Exception as ex:
        print(f"  [OpenTopography] download error: {ex}")
    return None
