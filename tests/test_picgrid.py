"""Shared PICGRID serialization and precision regression tests."""

import importlib.machinery
import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_core():
    """! @brief Load the conductor module. @return Loaded module. """
    path = ROOT / "picurv_cli" / "core.py"
    loader = importlib.machinery.SourceFileLoader("picurv_core_picgrid_tests", str(path))
    spec = importlib.util.spec_from_loader("picurv_core_picgrid_tests", loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def load_grid_generator():
    """! @brief Load the canonical grid generator. @return Loaded module. """
    path = ROOT / "generators" / "grid.gen"
    loader = importlib.machinery.SourceFileLoader("picurv_grid_picgrid_tests", str(path))
    spec = importlib.util.spec_from_loader("picurv_grid_picgrid_tests", loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def test_picgrid_writers_share_round_trip_safe_precision(tmp_path):
    """! @brief Guard coordinate precision in every PICGRID writer. @param[in] tmp_path Temporary path. """
    core = load_core()
    grid = load_grid_generator()
    settings = {"im": 4, "jm": 4, "km": 4, "xMins": 0.0, "xMaxs": 1.0,
                "yMins": 0.0, "yMaxs": 1.0, "zMins": 0.0, "zMaxs": 1.0,
                "rxs": 1.0, "rys": 1.0, "rzs": 1.0}
    generated, normalized = tmp_path/"generated.grid", tmp_path/"normalized.grid"
    core.generate_picgrid_from_programmatic_settings(settings, str(generated), 3.0)
    core.validate_and_nondimensionalize_picgrid(str(generated), str(normalized), 7.0, expected_nblk=1)
    direct = tmp_path / "direct.grid"
    coordinate = np.array([[[0.12345678901234567]]])
    grid.export_grid_to_picgrid(str(direct), coordinate, coordinate, coordinate)

    legacy = tmp_path / "legacy.grid"
    legacy.write_text(
        "1\n1 1 1\n0.12345678901234567\n0.23456789012345678\n0.34567890123456789\n",
        encoding="utf-8",
    )
    converted = tmp_path / "converted.grid"
    grid.convert_legacy_1d_to_picgrid(str(legacy), str(converted))

    for path in (generated, normalized, direct, converted):
        coordinate = path.read_text(encoding="utf-8").splitlines()[3].split()[0]
        assert len(coordinate.split("e")[0].replace(".", "").replace("-", "")) >= 18
