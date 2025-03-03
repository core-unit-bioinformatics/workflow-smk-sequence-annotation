
import pathlib

RUN_MINIMAP_REGIONDB = config.get("run_minimap_regiondb", False)

MINIMAP_REGION_DB = config.get("minimap_region_db", None)
MINIMAP_REGION_DB_NAMES = []

if MINIMAP_REGION_DB is not None:
    assert isinstance(MINIMAP_REGION_DB, dict)
    MINIMAP_REGION_DB_NAMES = sorted(MINIMAP_REGION_DB.keys())


RUN_MINIMAP_LABELREF = config.get("run_minimap_labelref", False)

MINIMAP_LABELED_REFERENCES = config.get("minimap_labeled_references", None)
MINIMAP_LABELED_REFERENCE_NAMES = []
if MINIMAP_LABELED_REFERENCES is not None:
    assert isinstance(MINIMAP_LABELED_REFERENCES, dict)
    MINIMAP_LABELED_REFERENCE_NAMES = sorted(MINIMAP_LABELED_REFERENCES.keys())

    for label_name, label_data in MINIMAP_LABELED_REFERENCES.items():
        assert "sequence" in label_data
        assert "labels" in label_data
