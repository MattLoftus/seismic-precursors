"""ANSS ComCat catalog fetch with full event records.

Returns a pandas DataFrame with columns:
    time (datetime64[ns, UTC]), latitude, longitude, depth_km, magnitude, eventid

Two entry points:
- `fetch_comcat_catalog`: circular (lat, lon, radius_km).
- `fetch_comcat_catalog_bbox`: rectangular (min/max lat/lon) for region polygons
   per `papers/pre_registration.md`.

Both chunk the request by year to stay within ComCat's ~20k-per-query limit and
to give visible progress for long pulls. Both support an optional CSV cache.
"""
from __future__ import annotations

import datetime as dt
import sys
import time
from pathlib import Path

import pandas as pd


def fetch_comcat_catalog(
    lat: float,
    lon: float,
    radius_km: float,
    start: dt.datetime,
    end: dt.datetime,
    m_min: float,
    chunk_years: int = 2,
    cache_path: Path | None = None,
    log: callable | None = print,
) -> pd.DataFrame:
    """Pull ANSS ComCat events around (lat, lon) within radius_km, [start, end), M >= m_min.

    Returns a DataFrame sorted by time ascending. If `cache_path` is provided and
    exists, loads from that CSV instead of re-fetching.
    """
    if cache_path is not None and Path(cache_path).is_file():
        if log:
            log(f"[data] loading cached catalog from {cache_path}")
        # ComCat times are mixed-format; use ISO8601 to handle both microsecond
        # and no-microsecond rows.
        df = pd.read_csv(cache_path)
        df["time"] = pd.to_datetime(df["time"], utc=True, format="ISO8601")
        return df

    from libcomcat.search import search

    if log:
        log(
            f"[data] ComCat: M>={m_min}, r={radius_km:.0f}km @ ({lat:.3f},{lon:.3f}), "
            f"{start.date()}..{end.date()}"
        )

    rows: list[dict] = []
    cur = start
    chunk = 0
    while cur < end:
        nxt = dt.datetime(min(cur.year + chunk_years, end.year), 1, 1)
        if nxt <= cur:
            nxt = end
        chunk += 1
        for attempt in (1, 2, 3):
            try:
                events = search(
                    starttime=cur,
                    endtime=nxt,
                    minmagnitude=m_min,
                    latitude=lat,
                    longitude=lon,
                    maxradiuskm=radius_km,
                )
                break
            except Exception as e:
                if attempt == 3:
                    raise
                wait = 5 * attempt
                if log:
                    log(
                        f"[data]   chunk {chunk} {cur.year}-{nxt.year} attempt {attempt}/3 "
                        f"failed ({type(e).__name__}); retrying in {wait}s"
                    )
                time.sleep(wait)

        for ev in events:
            if ev.magnitude is None:
                continue
            rows.append(
                {
                    "time": ev.time,
                    "latitude": float(ev.latitude),
                    "longitude": float(ev.longitude),
                    "depth_km": float(ev.depth) if ev.depth is not None else float("nan"),
                    "magnitude": float(ev.magnitude),
                    "eventid": str(ev.id),
                }
            )
        if log:
            log(
                f"[data]   chunk {chunk:>2d} {cur.year}-{nxt.year}: "
                f"{len(events):>5d} events  (running total: {len(rows)})"
            )
        cur = nxt

    df = pd.DataFrame(rows)
    if len(df):
        df["time"] = pd.to_datetime(df["time"], utc=True)
        df = df.sort_values("time").reset_index(drop=True)

    if cache_path is not None:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(cache_path, index=False)
        if log:
            log(f"[data] cached {len(df)} events -> {cache_path}")
    if log:
        log(
            f"[data] done -- {len(df)} events, "
            f"M [{df['magnitude'].min():.2f}, {df['magnitude'].max():.2f}]"
            if len(df)
            else "[data] done -- 0 events"
        )

    return df


def fetch_comcat_catalog_bbox(
    lat_min: float,
    lat_max: float,
    lon_min: float,
    lon_max: float,
    start: dt.datetime,
    end: dt.datetime,
    m_min: float,
    chunk_years: int = 1,
    cache_path: Path | None = None,
    log: callable | None = print,
) -> pd.DataFrame:
    """Pull ANSS ComCat events within a rectangular bounding box.

    Returns a DataFrame sorted by time ascending (same schema as
    `fetch_comcat_catalog`). Uses min/max lat/lon directly (no haversine).
    """
    if cache_path is not None and Path(cache_path).is_file():
        if log:
            log(f"[data] loading cached catalog from {cache_path}")
        df = pd.read_csv(cache_path)
        df["time"] = pd.to_datetime(df["time"], utc=True, format="ISO8601")
        return df

    from libcomcat.search import search

    if log:
        log(
            f"[data] ComCat bbox: M>={m_min}, "
            f"lat=[{lat_min:.2f},{lat_max:.2f}] lon=[{lon_min:.2f},{lon_max:.2f}], "
            f"{start.date()}..{end.date()}"
        )

    rows: list[dict] = []
    cur = start
    chunk = 0
    while cur < end:
        nxt = dt.datetime(min(cur.year + chunk_years, end.year), 1, 1)
        if nxt <= cur:
            nxt = end
        chunk += 1
        for attempt in (1, 2, 3):
            try:
                events = search(
                    starttime=cur,
                    endtime=nxt,
                    minmagnitude=m_min,
                    minlatitude=lat_min, maxlatitude=lat_max,
                    minlongitude=lon_min, maxlongitude=lon_max,
                )
                break
            except Exception as e:
                if attempt == 3:
                    raise
                wait = 5 * attempt
                if log:
                    log(
                        f"[data]   chunk {chunk} {cur.year}-{nxt.year} attempt {attempt}/3 "
                        f"failed ({type(e).__name__}); retrying in {wait}s"
                    )
                time.sleep(wait)
        for ev in events:
            if ev.magnitude is None:
                continue
            rows.append({
                "time": ev.time,
                "latitude": float(ev.latitude),
                "longitude": float(ev.longitude),
                "depth_km": float(ev.depth) if ev.depth is not None else float("nan"),
                "magnitude": float(ev.magnitude),
                "eventid": str(ev.id),
            })
        if log:
            log(
                f"[data]   chunk {chunk:>2d} {cur.year}-{nxt.year}: "
                f"{len(events):>6d} events  (running total: {len(rows)})"
            )
        cur = nxt

    df = pd.DataFrame(rows)
    if len(df):
        df["time"] = pd.to_datetime(df["time"], utc=True)
        df = df.sort_values("time").reset_index(drop=True)

    if cache_path is not None:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(cache_path, index=False)
        if log:
            log(f"[data] cached {len(df)} events -> {cache_path}")
    if log:
        log(
            f"[data] done -- {len(df)} events, "
            f"M [{df['magnitude'].min():.2f}, {df['magnitude'].max():.2f}]"
            if len(df)
            else "[data] done -- 0 events"
        )
    return df
