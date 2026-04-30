"""Region polygons and station assignments per `papers/pre_registration.md`.

This module is the single source of truth for which lat/lon bounding box and
which broadband stations belong to each region. Every Round-C+ experiment
script imports `REGIONS` from here rather than duplicating literals.

Pre-reg SHA: a4f1c6fa093ecea8bd00077fb8dbcab008c70eaa
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Station:
    """FDSN station identifier."""
    network: str
    station: str

    def code(self) -> str:
        return f"{self.network}.{self.station}"


@dataclass(frozen=True)
class Region:
    """A pre-registered region: bounding box + stations + tectonic regime.

    `fdsn_client` is the ObsPy FDSN service short-name for the region's primary
    waveform host. Defaults to "IRIS" (global IU/II network); California
    overrides to "NCEDC" (Berkeley/NCSN archive); Italy overrides to "INGV"
    (Italian national archive). This field was added in Session 12 after exp10
    discovered BK.PKD waveforms are not in EarthScope/IRIS for most of the
    2000–2024 period — they are NCEDC-hosted.
    """
    name: str
    lat_min: float
    lat_max: float
    lon_min: float
    lon_max: float
    tectonic_regime: str
    primary: Station
    backups: tuple[Station, ...] = field(default_factory=tuple)
    is_test: bool = False
    fdsn_client: str = "IRIS"

    def contains(self, lat: float, lon: float) -> bool:
        return (self.lat_min <= lat <= self.lat_max and
                self.lon_min <= lon <= self.lon_max)

    def all_stations(self) -> list[Station]:
        return [self.primary, *self.backups]


# === Training regions (6) ===
CALIFORNIA = Region(
    name="California",
    lat_min=32.0, lat_max=42.0, lon_min=-125.0, lon_max=-114.0,
    tectonic_regime="strike-slip + thrust margins",
    primary=Station("BK", "PKD"),
    backups=(Station("BK", "CMB"), Station("IU", "COR")),
    fdsn_client="NCEDC",
)
CASCADIA = Region(
    name="Cascadia",
    lat_min=40.0, lat_max=50.0, lon_min=-130.0, lon_max=-121.0,
    tectonic_regime="subduction (Juan de Fuca + Pacific)",
    primary=Station("UW", "LON"),
    backups=(Station("IU", "COR"), Station("UW", "RATT")),
)
JAPAN = Region(
    name="Japan",
    lat_min=30.0, lat_max=46.0, lon_min=128.0, lon_max=148.0,
    tectonic_regime="subduction (Pacific + Philippine)",
    primary=Station("II", "MAJO"),
    backups=(Station("IU", "MAJO"), Station("II", "INU")),
)
CHILE = Region(
    name="Chile",
    lat_min=-45.0, lat_max=-18.0, lon_min=-76.0, lon_max=-68.0,
    tectonic_regime="subduction (Nazca)",
    primary=Station("IU", "LCO"),
    backups=(Station("C1", "GO01"), Station("G", "PEL")),
)
TURKEY = Region(
    name="Turkey",
    lat_min=36.0, lat_max=42.0, lon_min=26.0, lon_max=44.0,
    tectonic_regime="strike-slip (NAF + EAF) + thrust",
    primary=Station("IU", "ANTO"),
    backups=(Station("KO", "ANTO"), Station("II", "KIV")),
)
ITALY = Region(
    name="Italy",
    lat_min=36.0, lat_max=47.0, lon_min=6.0, lon_max=19.0,
    tectonic_regime="continental thrust + extension",
    # Pre-reg v1 §3 listed IV.AQU as primary. Session 13 verification found
    # IV.AQU has zero data availability via INGV, ORFEUS, or IRIS — the
    # post-2009 network reorganization moved L'Aquila Observatory's broadband
    # data under MN.AQU (MedNet, INGV-affiliated, IRIS-archived). MN.AQU is
    # the same physical station; the network-code update is documented as
    # an operational implementation of pre-reg v1 §3's "95% availability
    # fallback rule" rather than a station swap.
    primary=Station("MN", "AQU"),
    backups=(Station("IV", "AQU"), Station("IV", "MGAB"), Station("G", "SSB")),
    fdsn_client="IRIS",
)

# === Test regions (held out) (2) ===
MEXICO = Region(
    name="Mexico",
    lat_min=14.0, lat_max=32.0, lon_min=-118.0, lon_max=-86.0,
    tectonic_regime="subduction (Cocos) + Gulf of California rift",
    primary=Station("IU", "TEIG"),
    backups=(Station("IU", "TUC"), Station("G", "UNM")),
    is_test=True,
)
ALASKA = Region(
    name="Alaska",
    lat_min=51.0, lat_max=72.0, lon_min=-180.0, lon_max=-130.0,
    tectonic_regime="subduction (Pacific) + transform",
    primary=Station("IU", "COLA"),
    backups=(Station("IU", "KDAK"), Station("AT", "KIAG")),
    is_test=True,
)


REGIONS: dict[str, Region] = {
    r.name: r for r in (
        CALIFORNIA, CASCADIA, JAPAN, CHILE, TURKEY, ITALY,   # training
        MEXICO, ALASKA,                                       # test
    )
}

TRAINING_REGIONS: list[Region] = [r for r in REGIONS.values() if not r.is_test]
TEST_REGIONS: list[Region] = [r for r in REGIONS.values() if r.is_test]


def get_region(name: str) -> Region:
    """Look up a region by name (case-insensitive)."""
    for k, v in REGIONS.items():
        if k.lower() == name.lower():
            return v
    raise KeyError(f"unknown region '{name}'; known: {list(REGIONS.keys())}")
