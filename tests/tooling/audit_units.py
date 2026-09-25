#!/usr/bin/env python3
"""!
@file audit_units.py
@brief Enforce the physical-units rule and the published units indexes.

Every configuration input is physical and is converted to solver units exactly once;
every catalogued field records the dimension post-processing scales it by. This checks
that no solver-facing input lacks a recorded dimension, and that both indexes on page 19
match what the code records: the input index against `INPUT_QUANTITIES` and its sibling
tables in `picurv_cli/core.py`, the field index against the `FIELD_DIM_*` argument of
every entry in `src/field_catalog.c` and `src/particle_field_catalog.c`.

It does not check that a conversion is performed where the table says; the non-unit
ingress tests and the units-equivalence smoke run check that.
"""

from __future__ import annotations

import importlib.machinery
import importlib.util
import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
CORE = REPO_ROOT / "picurv_cli" / "core.py"
IC_GENERATOR = REPO_ROOT / "generators" / "ic.gen"
EULERIAN_CATALOG = REPO_ROOT / "src" / "field_catalog.c"
PARTICLE_CATALOG = REPO_ROOT / "src" / "particle_field_catalog.c"
PAGE = REPO_ROOT / "docs" / "pages" / "19_Nondimensionalization.md"

#: Schemas whose keys configure the solver. Cluster, study, and workspace files hold
#: scheduling and bookkeeping; study overrides address case keys, which are covered.
SOLVER_INPUT_SCHEMAS = {
    "case": "_CASE_SCHEMA",
    "solver": "_SOLVER_SCHEMA",
    "monitor": "_MONITOR_SCHEMA",
    "post": "_POST_SCHEMA",
}

#: Free-form mappings whose keys have their own dimension tables.
DEDICATED_TABLE_SURFACES = {
    ("case", "boundary_conditions", "[]", "params"),
    ("case", "boundary_conditions", "[]", "[]", "params"),
    ("case", "properties", "initial_conditions", "params"),
}

#: Built-in initial-condition generators, which are not Python providers.
BUILTIN_IC_GENERATORS = {"zero", "constant", "streamwise_constant", "poiseuille"}

ENTRY_RE = re.compile(r"\b(?:FIELD_ENTRY|FIELD_COORDINATE_ENTRY|PARTICLE_FIELD_ENTRY)\(")
DIMENSION_RE = re.compile(r"FIELD_DIM_(\w+)")
NAME_RE = re.compile(r'"([^"]+)"')
INDEX_ROW_RE = re.compile(r"^\|\s*`([^`]+)`\s*\|\s*(\w+)\s*\|\s*`?(\w+)`?\s*\|", re.M)


def load_module(name: str, path: Path):
    """!
    @brief Load a repository script as an importable module.
    @param[in] name Module name to register.
    @param[in] path Source path of the script.
    @return Loaded module object.
    """
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def schema_leaves(schema: dict):
    """!
    @brief Yield every key path a schema accepts that is not itself a nested mapping.
    @details A path mapped to `None` accepts arbitrary keys and is yielded as a leaf: its
             contents are covered as a whole or by a dedicated table. A list's children
             are keyed under `path + ("[]",)`, so any path prefixing an entry is interior.
    @param[in] schema Key schema mapping parent paths to allowed key sets.
    @return Generator of key-path tuples.
    """
    interiors = {
        parent[:length]
        for parent, keys in schema.items() if keys is not None
        for length in range(1, len(parent) + 1)
    }
    for parent, keys in schema.items():
        for key in keys or ():
            path = parent + (key,)
            if path not in interiors:
                yield path


def uncovered_inputs(core) -> list:
    """!
    @brief Name every solver-facing input with no recorded physical dimension.
    @param[in] core Loaded conductor core.
    @return Dotted names of uncovered inputs; empty when complete.
    """
    problems = []
    for role, schema_name in SOLVER_INPUT_SCHEMAS.items():
        for path in schema_leaves(getattr(core, schema_name)):
            full = (role,) + path
            if full in DEDICATED_TABLE_SURFACES:
                continue
            try:
                core.input_quantity(full)
            except KeyError:
                problems.append(".".join(full))

    accepted = set(core._DEPRECATED_BC_PARAM_ALIASES)
    for spec in core.BC_HANDLER_SPECS.values():
        accepted |= spec["required_params"] | spec["optional_params"]
    problems.extend(f"boundary_conditions[].params.{key}"
                    for key in sorted(accepted - set(core.BC_PARAM_QUANTITIES)))

    generators = BUILTIN_IC_GENERATORS | core._PYTHON_INITIAL_CONDITION_PROVIDERS
    problems.extend(f"initial_conditions.generator {name}"
                    for name in sorted(generators ^ set(core.IC_PARAM_QUANTITIES)))
    ic_gen = load_module("audit_units_ic_gen", IC_GENERATOR)
    provider_params = {
        "spectral_random_velocity": core.SPECTRAL_RANDOM_VELOCITY_PARAMS,
        "channel_spectral_velocity": ic_gen.WALL_SPECTRAL_PARAMS,
        "duct_spectral_velocity": ic_gen.WALL_SPECTRAL_PARAMS,
    }
    for generator, keys in provider_params.items():
        roots = {name.split(".", 1)[0] for name in core.IC_PARAM_QUANTITIES.get(generator, {})}
        problems.extend(f"initial_conditions.params.{key} ({generator})" for key in sorted(set(keys) - roots))
    return problems


def dimension_name(core, dimension) -> str:
    """!
    @brief Name a dimension triple by its constant in the conductor core.
    @param[in] core Loaded conductor core.
    @param[in] dimension `(length, velocity, density)` exponents.
    @return Constant name, e.g. `VELOCITY`.
    """
    for name in ("DIMENSIONLESS", "LENGTH", "VELOCITY", "TIME", "WAVENUMBER", "VOLUME_FLUX",
                 "DIFFUSIVITY", "PRESSURE", "DENSITY", "DYNAMIC_VISCOSITY"):
        if getattr(core, name) == dimension:
            return name
    raise ValueError(f"dimension {dimension} has no named constant")


def physical_inputs(core) -> dict:
    """!
    @brief Every input that carries a physical dimension, as the input index lists it.
    @details Dimensionless and non-quantity inputs are omitted: there is nothing to
             convert. Payload inputs whose dimension follows their content are listed
             with dimension `PAYLOAD`.
    @param[in] core Loaded conductor core.
    @return Mapping of documented input name to `(dimension name, conversion site)`.
    """
    rows = {}

    def record(name, quantity):
        """!
        @brief Add one input to the index when it carries something to convert.
        @param[in] name Documented input name.
        @param[in] quantity `(dimension, site)` entry.
        """
        dimension, site = quantity
        if site in ("", "passthrough"):
            return
        if dimension is not None and dimension == core.DIMENSIONLESS:
            return
        rows[name] = ("PAYLOAD" if dimension is None else dimension_name(core, dimension), site)

    for path, quantity in core.INPUT_QUANTITIES.items():
        record(f"{path[0]}.yml: " + ".".join(path[1:]), quantity)
    for key, quantity in core.BC_PARAM_QUANTITIES.items():
        record(f"case.yml: boundary_conditions[].params.{key}", quantity)
    for generator, table in core.IC_PARAM_QUANTITIES.items():
        for key, quantity in table.items():
            record(f"case.yml: initial_conditions.params.{key} ({generator})", quantity)
    return rows


def compiled_field_dimensions() -> dict:
    """!
    @brief Field name to `(catalog, FIELD_DIM_* suffix)` for every compiled entry.
    @return Mapping of canonical field names to their catalog and dimension.
    @throws ValueError when an entry carries no dimension argument.
    """
    fields = {}
    for catalog, path in (("Eulerian", EULERIAN_CATALOG), ("particle", PARTICLE_CATALOG)):
        text = path.read_text(encoding="utf-8")
        body = text[text.index("gFieldCatalog" if catalog == "Eulerian" else "gParticleFieldCatalog"):]
        for match in ENTRY_RE.finditer(body):
            depth, index = 1, match.end()
            while depth:
                depth += {"(": 1, ")": -1}.get(body[index], 0)
                index += 1
            entry = body[match.end():index - 1]
            name = NAME_RE.search(entry).group(1)
            dimension = DIMENSION_RE.search(entry)
            if dimension is None:
                raise ValueError(f"{catalog} field '{name}' has no FIELD_DIM_* argument")
            fields[(catalog, name)] = dimension.group(1)
    return fields


def page_section(heading: str) -> str:
    """!
    @brief One `@section` body of page 19.
    @param[in] heading Section anchor to extract.
    @return Section text up to the next section.
    """
    text = PAGE.read_text(encoding="utf-8")
    start = text.index(f"@section {heading}")
    following = text.find("@section", start + 1)
    return text[start: following if following != -1 else len(text)]


def documented_rows(heading: str) -> dict:
    """!
    @brief Rows of the first table in a section: first column to second and third.
    @param[in] heading Section anchor.
    @return Mapping of backticked first-column text to the next two cells.
    """
    return {(first, second): third
            for first, second, third in INDEX_ROW_RE.findall(page_section(heading))}


def main() -> int:
    """!
    @brief Fail when an input lacks a dimension or an index disagrees with the code.
    @return Process status code.
    """
    core = load_module("audit_units_core", CORE)
    problems = [f"no recorded dimension: {name}" for name in uncovered_inputs(core)]

    inputs = physical_inputs(core)
    documented_inputs = {name: (dimension, site)
                         for (name, dimension), site in documented_rows("p19_inputs_sec").items()}
    for name in sorted(set(inputs) - set(documented_inputs)):
        problems.append(f"input index: '{name}' is physical but not listed")
    for name in sorted(set(documented_inputs) - set(inputs)):
        problems.append(f"input index: '{name}' is listed but records nothing to convert")
    for name in sorted(set(inputs) & set(documented_inputs)):
        if inputs[name] != documented_inputs[name]:
            problems.append(f"input index: '{name}' records {inputs[name]}, page says {documented_inputs[name]}")

    fields = compiled_field_dimensions()
    documented_fields = {(catalog, name): dimension
                         for (name, catalog), dimension in documented_rows("p19_fields_sec").items()}
    for key in sorted(set(fields) - set(documented_fields)):
        problems.append(f"field index: {key[0]} field '{key[1]}' is not listed")
    for key in sorted(set(documented_fields) - set(fields)):
        problems.append(f"field index: {key[0]} field '{key[1]}' is not in the compiled catalog")
    for key in sorted(set(fields) & set(documented_fields)):
        if fields[key] != documented_fields[key]:
            problems.append(f"field index: {key[1]} is FIELD_DIM_{fields[key]}, page says {documented_fields[key]}")

    if problems:
        print("Units rule or its published indexes do not match the code:", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        print("\nRecord a new input in picurv_cli/core.py INPUT_QUANTITIES (or BC_PARAM_QUANTITIES /\n"
              "IC_PARAM_QUANTITIES), give a new field a FIELD_DIM_* argument, and update the\n"
              "indexes in docs/pages/19_Nondimensionalization.md.", file=sys.stderr)
        return 1
    print(f"Units audit passed: every solver-facing input records a dimension; {len(inputs)} "
          f"physical inputs and {len(fields)} catalogued fields match page 19. Whether each "
          "conversion is performed is checked by the ingress tests and smoke run, not here.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
