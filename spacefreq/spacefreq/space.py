"""Operations to calculate the landscape area occupied by a variant."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
from pyproj import Geod
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon
from shapely.geometry import polygon as shapely_polygon


_GEOD = Geod(ellps="WGS84")


def get_carrier_locations(
    genotype_calls: np.ndarray,
    metadata: pd.DataFrame,
    allele: int,
) -> np.ndarray:
    """Return the unique (x, y) coordinates of samples carrying *allele*.

    A sample is considered a carrier if it has at least one copy of the
    allele (heterozygous or homozygous alternate).  At multi-allelic sites
    only the exact allele index is matched, so carriers of a different
    alternate are never counted.

    Parameters
    ----------
    genotype_calls:
        2-D array of shape (n_samples, ploidy) for a single variant.
    metadata:
        DataFrame with columns ``sampleID``, ``x`` (longitude), ``y`` (latitude).
        Row order must correspond to the sample order in *genotype_calls*.
    allele:
        Allele index to query (1 = first alternate, 2 = second alternate, …).

    Returns
    -------
    np.ndarray of shape (n_unique_locations, 2) with columns [x, y].
    """
    # mask: True for samples that carry at least one copy of the allele
    mask = np.sum(genotype_calls == allele, axis=1) >= 1
    carriers = metadata[mask]
    if len(carriers) == 0:
        return np.empty((0, 2))
    locs = np.unique(carriers[["x", "y"]].to_numpy(), axis=0)
    return locs


def calculate_area(
    coords: np.ndarray,
    min_locs: int = 3,
    transect: float = 1.0,
    sample_area: float = 1.0,
) -> float:
    """Calculate the landscape area (km²) occupied by a set of coordinates.

    Dispatch rules (n = number of unique locations):
    - n < min_locs  → NaN
    - n == 1        → *sample_area*
    - n == 2        → geodetic distance (km) × *transect*
    - n >= 3        → area of convex hull (km²) via WGS84 geodesy

    Parameters
    ----------
    coords:
        Array of shape (n, 2) with columns [x (longitude), y (latitude)].
    min_locs:
        Minimum number of unique locations required to report an area.
    transect:
        Width (km) of the transect used when n == 2.
    sample_area:
        Default area (km²) for a single sampling location.

    Returns
    -------
    float area in km², or NaN when n < min_locs.
    """
    n = len(coords)

    if n < min_locs:
        return math.nan

    if n == 1:
        return float(sample_area)

    if n == 2:
        return _geodetic_distance_km(coords) * transect

    return _convex_hull_area_km2(coords)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _geodetic_distance_km(coords: np.ndarray) -> float:
    """Geodetic distance in km between two (lon, lat) points."""
    lon1, lat1 = coords[0, 0], coords[0, 1]
    lon2, lat2 = coords[1, 0], coords[1, 1]
    _, _, dist_m = _GEOD.inv(lon1, lat1, lon2, lat2)
    return abs(dist_m) / 1000.0


def _convex_hull_area_km2(coords: np.ndarray) -> float:
    """Area of the convex hull over (lon, lat) points in km²."""
    hull = ConvexHull(coords)
    hull_coords = coords[hull.vertices]
    # Close the polygon ring
    poly = Polygon(list(zip(hull_coords[:, 0], hull_coords[:, 1])))
    poly = shapely_polygon.orient(poly)
    area_m2, _ = _GEOD.geometry_area_perimeter(poly)
    return abs(area_m2) / 1_000_000.0
