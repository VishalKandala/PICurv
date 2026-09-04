"""!
@file picurv_cli/storage/__init__.py
@brief PICurv run, study, and workspace storage lifecycle.

The package is split by responsibility - models, transport, safety, inventory,
packaging, compatibility, catalog, operations - and re-exports the whole surface,
so importers and tests address `picurv_cli.storage` exactly as before.
"""

from .models import *  # noqa: F401,F403
from .transport import *  # noqa: F401,F403
from .safety import *  # noqa: F401,F403
from .inventory import *  # noqa: F401,F403
from .packaging import *  # noqa: F401,F403
from .compatibility import *  # noqa: F401,F403
from .catalog import *  # noqa: F401,F403
from .operations import *  # noqa: F401,F403

from . import models  # noqa: F401
from . import transport  # noqa: F401
from . import safety  # noqa: F401
from . import inventory  # noqa: F401
from . import packaging  # noqa: F401
from . import compatibility  # noqa: F401
from . import catalog  # noqa: F401
from . import operations  # noqa: F401

# Private helpers the conductor and the test suite address directly.
from .models import (
    _PARALLEL_GZIP_VALUES,
    _atomic_write_json,
    _find_upwards,
    _human_bytes,
    _parse_tags,
    _read_json,
    _sha256_file,
    _state_path,
    _utc_now,
)
from .transport import (
    _blob_remote,
    _chunk_remote_path,
    _load_remote_manifest,
    _object_remote,
    _read_remote_bytes,
    _remote_blob_present,
    _remote_join,
    _remote_sha256,
    _run_rclone,
    _upload_verified,
)
from .safety import (
    _assert_archive_safe,
    _collect_job_ids,
    _lock_owner_active,
    _require_free_space,
    _slurm_activity,
)
from .inventory import (
    _artifact_runtime_roots,
    _checkpoint_step_from_bundle,
    _artifact_component_layout,
    _checkpoint_steps,
    _classify_component,
    _classify_workspace_component,
    _completed_study_members,
    _discover_dependencies,
    _discover_external_paths,
    _find_incomplete_checkpoints,
    _inventory_fingerprint,
    _path_is_within,
    _resolve_configured_path,
    _select_completed_case_ids,
    _validate_case_id,
    _walk_archive_entries,
)
from .packaging import (
    _build_chunk_specs,
    _chunk_extension,
    _compression_size_range,
    _entry_retained_by_policy,
    _extract_chunk,
    _merge_tree,
    _resolve_offload_policy,
    _safe_component_name,
    _select_compression,
    _validate_tar_members,
    _write_tar_chunk,
)
from .compatibility import (
    _capture_parameter_summary,
    _capture_run_assets,
    _capture_study_context,
    _capture_workspace_assets,
    _config_fingerprints,
    _git_provenance,
    _rebase_restored_text_paths,
    _restore_study_context,
)
from .catalog import (
    _find_reusable_archive,
)
from .operations import (
    _download_archive_components,
    _expand_checkpoint_selection,
    _prune_archived_payload,
    _register_existing_archive,
    _render_plan,
    _render_status,
    _resolve_archive_id_from_args,
)
