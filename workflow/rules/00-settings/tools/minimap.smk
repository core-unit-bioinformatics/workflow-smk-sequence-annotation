
import pathlib

RUN_MINIMAP_REGIONDB = config.get("run_minimap_regiondb", False)

MINIMAP_REGION_DB = config.get("minimap_region_db", None)

if MINIMAP_REGION_DB is not None:
    assert isinstance(MINIMAP_REGION_DB, dict)
    MINIMAP_REGION_DB_NAMES = sorted(MINIMAP_REGION_DB.keys())
