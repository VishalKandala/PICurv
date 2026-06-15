"""Argument parser construction and command dispatch for PICurv."""

from .core import *


# ==============================================================================
# MAIN COMMAND-LINE INTERFACE PARSER
# ==============================================================================

def _add_run_parser(subparsers):
    """!
    @brief Attach `run` parser with staged execution and dry-run support.
    @param[in] subparsers Argument passed to `_add_run_parser()`.
    @return Value returned by `_add_run_parser()`.
    """
    p_run = subparsers.add_parser(
        "run",
        help="Execute a simulation workflow (solve and/or post-process).",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Execute solver and/or post-processing stages.\n\n"
            "Notes:\n"
            "  - --num-procs applies to solver stage launches.\n"
            "  - Post-processing defaults to a single MPI rank/task.\n"
            "  - With --solve, --continue resumes the existing run directory in-place.\n"
            "  - With --post-process, --continue resumes the same recipe from the first unfinished step\n"
            "    and caps the launch to the highest fully available contiguous source frontier.\n\n"
            "Diagnostics:\n"
            "  - PETSc and runtime memory diagnostics live under monitor.yml -> diagnostics.\n"
            "  - Use --dry-run to inspect resolved PETSc flags and expected log destinations.\n\n"
            "Examples:\n"
            "  picurv run --solve -n 8 --case case.yml --solver solver.yml --monitor monitor.yml\n"
            "  picurv run --solve --restart-from runs/old_run --case case.yml --solver solver.yml --monitor monitor.yml\n"
            "  picurv run --solve --continue --run-dir runs/my_run --case case.yml --solver solver.yml --monitor monitor.yml\n"
            "  picurv run --post-process --run-dir runs/my_run --post post.yml\n"
            "  picurv run --post-process --continue --run-dir runs/my_run --post post.yml\n"
            "  picurv run --solve --case case.yml --solver solver.yml --monitor monitor.yml --dry-run"
        ),
        epilog="Next: run `picurv validate ...` first for config-only checks.",
    )
    run_group = p_run.add_argument_group("stages")
    run_group.add_argument("--solve", action="store_true", help="Execute the solver stage (creates a new run directory).")
    run_group.add_argument("--post-process", action="store_true", help="Execute the post-processing stage on a run directory.")

    solver_group = p_run.add_argument_group("solver inputs (required for --solve)")
    solver_group.add_argument("--case", help="Path to the case definition file (e.g., case.yml).")
    solver_group.add_argument("--solver", help="Path to the solver settings profile (e.g., solver.yml).")
    solver_group.add_argument("--monitor", help="Path to the monitoring, diagnostics, and I/O profile (e.g., monitor.yml).")
    solver_group.add_argument(
        "--restart-from",
        help="Path to an existing run directory to restart from.\n"
             "Creates a new run directory and copies/references checkpoint data from the source.",
    )
    run_group.add_argument(
        "--continue",
        action="store_true",
        dest="continue_run",
        help="Resume an existing run directory in-place. Requires --run-dir.\n"
             "With --solve, appends to existing solver output/logs.\n"
             "With --post-process, resumes the same recipe from the first unfinished step\n"
             "and skips already-complete work inside the current live source frontier.",
    )

    post_group = p_run.add_argument_group("post-processor inputs (required for --post-process)")
    post_group.add_argument("--run-dir", help="Path to an existing run directory.\n(Used with --post-process or --continue).")
    post_group.add_argument("--post", help="Path to the post-processing recipe file (e.g., post.yml).")

    p_run.add_argument(
        "-n",
        "--num-procs",
        type=int,
        default=1,
        help="Number of MPI processes for the solver stage. Post-processing defaults to 1 rank.",
    )
    p_run.add_argument("--cluster", help="Path to cluster.yml for Slurm execution mode.")
    p_run.add_argument("--scheduler", help="Explicit scheduler selector (currently 'slurm').")
    p_run.add_argument("--no-submit", action="store_true", help="Stage run artifacts without starting local execution or Slurm submission.")
    p_run.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve and print planned commands/artifacts, including diagnostic flags and log paths, without writing files.",
    )
    p_run.add_argument(
        "--format",
        dest="output_format",
        choices=["text", "json"],
        default="text",
        help="Output format for --dry-run (default: text).",
    )
    return p_run


def _add_sweep_parser(subparsers):
    """!
    @brief Perform add sweep parser.
    @param[in] subparsers Argument passed to `_add_sweep_parser()`.
    @return Value returned by `_add_sweep_parser()`.
    """
    p_sweep = subparsers.add_parser(
        "sweep",
        help="Launch or continue a Slurm-based parameter sweep/study.",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Launch studies from study.yml + cluster.yml, whether the study uses\n"
            "cross-product parameter combinations or explicit coupled parameter_sets, continue\n"
            "partially-completed studies, or re-aggregate metrics.\n\n"
            "Examples:\n"
            "  picurv sweep --study study.yml --cluster cluster.yml\n"
            "  picurv sweep --study study.yml --cluster cluster.yml --no-submit\n"
            "  picurv sweep --continue --study-dir studies/<id>\n"
            "  picurv sweep --continue --study-dir studies/<id> --cluster cluster_more_time.yml\n"
            "  picurv sweep --reaggregate --study-dir studies/<id>"
        ),
        epilog=(
            "Study files support either `parameters` cross-product expansion or explicit `parameter_sets`.\n"
            "Next: inspect studies/<study_id>/results/metrics_table.csv for aggregated metrics."
        ),
    )
    p_sweep.add_argument("--study", help="Path to study.yml defining either a `parameters` cross-product expansion or explicit parameter_sets, plus metrics.")
    p_sweep.add_argument("--cluster", help="Path to cluster.yml defining Slurm resources.")
    p_sweep.add_argument("--no-submit", action="store_true", help="Generate all study artifacts without submitting jobs.")
    p_sweep.add_argument("--study-dir", help="Path to an existing study directory (for --continue/--reaggregate).")
    p_sweep.add_argument("--continue", action="store_true", dest="continue_study",
                         help="Continue a partially-completed study. Requires --study-dir.")
    p_sweep.add_argument("--reaggregate", action="store_true",
                         help="Re-run metrics aggregation on a completed study. Requires --study-dir.")
    return p_sweep


def _add_validate_parser(subparsers):
    """!
    @brief Perform add validate parser.
    @param[in] subparsers Argument passed to `_add_validate_parser()`.
    @return Value returned by `_add_validate_parser()`.
    """
    p_validate = subparsers.add_parser(
        "validate",
        help="Validate config files without launching solver/post.",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Validate one or more config roles. No solver/post execution and no run/study artifact writes.\n\n"
            "Examples:\n"
            "  picurv validate --case case.yml --solver solver.yml --monitor monitor.yml\n"
            "  picurv validate --post post.yml --cluster cluster.yml\n"
            "  picurv validate --study study.yml --cluster cluster.yml --strict"
        ),
        epilog="Next: run `picurv run --dry-run ...` to inspect resolved commands/artifacts.",
    )
    p_validate.add_argument("--case", help="Path to case.yml")
    p_validate.add_argument("--solver", help="Path to solver.yml")
    p_validate.add_argument("--monitor", help="Path to monitor.yml (logging, diagnostics, profiling, and I/O)")
    p_validate.add_argument("--restart-from", help="Path to a run directory to validate restart from.")
    p_validate.add_argument("--continue", action="store_true", dest="continue_run", help="Validate continue-in-place mode. Requires --run-dir.")
    p_validate.add_argument("--run-dir", help="Path to an existing run directory (for --continue validation).")
    p_validate.add_argument("--post", help="Path to post.yml")
    p_validate.add_argument("--cluster", help="Path to cluster.yml")
    p_validate.add_argument("--study", help="Path to study.yml")
    p_validate.add_argument("--strict", action="store_true", help="Enable additional strict checks for selected roles.")
    return p_validate

def _add_precompute_parser(subparsers):
    """!
    @brief Attach `precompute` parser for deterministic artifact generation.
    @param[in] subparsers Argument passed to `_add_precompute_parser()`.
    @return Configured parser.
    """
    p_precompute = subparsers.add_parser(
        "precompute",
        help="Generate deterministic case artifacts without launching the solver.",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Generate configured deterministic artifacts, such as grid_gen grids,\n"
            "generated prescribed-flow inlet profiles, and ic_gen initial conditions,\n"
            "into a run-like config layout.\n\n"
            "Examples:\n"
            "  picurv precompute --case case.yml\n"
            "  picurv precompute --case case.yml --output-dir precomputed/channel"
        ),
        epilog="Next: inspect generated artifacts or run solve with the configured generated modes directly.",
    )
    p_precompute.add_argument("--case", required=True, help="Path to case.yml containing grid/profile/IC generator settings.")
    p_precompute.add_argument(
        "--output-dir",
        help="Directory where precomputed artifacts should be written. Defaults to precomputed/<case-name>.",
    )
    return p_precompute


def _add_summarize_parser(subparsers):
    """!
    @brief Attach `summarize` parser for read-only run-health views.
    @param[in] subparsers Argument passed to `_add_summarize_parser()`.
    @return Value returned by `_add_summarize_parser()`.
    """
    p_summarize = subparsers.add_parser(
        "summarize",
        help="Summarize run configs/health and plot scalar log histories.",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Build read-only configuration overviews, run-health summaries, and scalar time-history plots.\n"
            "It does not modify solver output and works for active or completed runs.\n\n"
            "Examples:\n"
            "  picurv summarize --run-dir runs/my_run --overview\n"
            "  picurv summarize --run-dir runs/my_run --case --solver\n"
            "  picurv summarize --run-dir runs/my_run --latest\n"
            "  picurv summarize --run-dir runs/my_run --monitor --step 500\n"
            "  picurv summarize --run-dir runs/my_run --list-plot-series\n"
            "  picurv summarize --run-dir runs/my_run --plot momentum.residual_norm --last 100\n"
            "  picurv summarize --run-dir runs/my_run --latest --format json"
        ),
        epilog="Next: use `picurv run ...` to create runs or `picurv sweep ...` for multi-case studies.",
    )
    p_summarize.add_argument("--run-dir", required=True, help="Path to the run directory to inspect.")
    p_summarize.add_argument(
        "--overview",
        action="store_true",
        help="Summarize run metadata plus copied case, solver, and monitor configs without implicitly requesting health.",
    )
    p_summarize.add_argument("--case", action="store_true", help="Summarize the copied run-local case.yml.")
    p_summarize.add_argument("--solver", action="store_true", help="Summarize the copied run-local solver.yml.")
    p_summarize.add_argument("--monitor", action="store_true", help="Summarize the copied run-local monitor.yml.")
    plot_group = p_summarize.add_mutually_exclusive_group()
    plot_group.add_argument(
        "--list-plot-series",
        dest="list_plot_series",
        action="store_true",
        help="List scalar time histories available to --plot.",
    )
    plot_group.add_argument(
        "--list-series",
        dest="list_plot_series",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    plot_group.add_argument(
        "--plot",
        dest="plot_series",
        help="Plot one qualified scalar history, such as momentum.residual_norm (requires matplotlib).",
    )
    p_summarize.add_argument("--last", dest="last_n", type=int, help="Plot only the last N chronological records per line.")
    p_summarize.add_argument("--plot-output", help="Save the plot to this path instead of opening an interactive window.")
    p_summarize.add_argument("--linear-y", action="store_true", help="Force linear y-axis scaling instead of automatic residual/norm log scaling.")
    step_group = p_summarize.add_mutually_exclusive_group()
    step_group.add_argument("--step", type=int, help="Specific completed timestep to summarize.")
    step_group.add_argument(
        "--latest",
        action="store_true",
        help="Summarize the most recently appended completed step found in available artifacts (default behavior).",
    )
    step_group.add_argument(
        "--max-step",
        action="store_true",
        help="Summarize the numerically largest timestep found in available artifacts.",
    )
    p_summarize.add_argument(
        "--snapshot-rows",
        type=int,
        default=5,
        help="Number of sampled particle snapshot rows to preview when solver stream output contains them.",
    )
    p_summarize.add_argument(
        "--format",
        dest="output_format",
        choices=["text", "json"],
        default="text",
        help="Output format (default: text).",
    )
    return p_summarize


def _add_cancel_parser(subparsers):
    """!
    @brief Perform add cancel parser.
    @param[in] subparsers Argument passed to `_add_cancel_parser()`.
    @return Value returned by `_add_cancel_parser()`.
    """
    p_cancel = subparsers.add_parser(
        "cancel",
        help="Cancel Slurm-submitted jobs for an existing run directory.",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Look up scheduler/submission.json inside an existing run directory and cancel\n"
            "the recorded Slurm job(s) without requiring manual job-id lookup.\n\n"
            "Examples:\n"
            "  picurv cancel --run-dir runs/my_run\n"
            "  picurv cancel --run-dir runs/my_run --stage solve\n"
            "  picurv cancel --run-dir runs/my_run --stage solve --graceful\n"
            "  picurv cancel --run-dir runs/my_run --dry-run\n\n"
            "Default cancellation is a hard Slurm cancel (`scancel <job_id>`). With\n"
            "`--graceful`, solver jobs receive SIGUSR1 so PICurv can write the latest\n"
            "safe off-cadence step at the next runtime checkpoint before exiting.\n"
            "Post-process jobs still use ordinary hard cancellation."
        ),
        epilog="Next: use `picurv summarize --run-dir ...` to inspect whatever output the run already produced.",
    )
    p_cancel.add_argument("--run-dir", required=True, help="Path to the run directory whose Slurm job(s) should be canceled.")
    p_cancel.add_argument(
        "--stage",
        choices=["all", "solve", "post-process"],
        default="all",
        help="Which recorded stage job(s) to cancel (default: all).",
    )
    p_cancel.add_argument(
        "--dry-run",
        action="store_true",
        help="Show which `scancel` command(s) would run without actually canceling anything.",
    )
    p_cancel.add_argument(
        "--graceful",
        action="store_true",
        help=(
            "For solver jobs, request a clean runtime shutdown by sending SIGUSR1 instead of "
            "hard-canceling immediately. The solver writes the latest safe off-cadence step at "
            "the next checkpoint. Non-solver stages still use ordinary scancel."
        ),
    )
    return p_cancel


def _add_submit_parser(subparsers):
    """!
    @brief Perform add submit parser.
    @param[in] subparsers Argument passed to `_add_submit_parser()`.
    @return Value returned by `_add_submit_parser()`.
    """
    p_submit = subparsers.add_parser(
        "submit",
        help="Execute or submit previously staged artifacts from an existing run or study directory.",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Consume an existing artifact set created by picurv --no-submit and\n"
            "execute/submit it later without regenerating configs or scripts.\n\n"
            "Examples:\n"
            "  picurv submit --run-dir runs/my_run\n"
            "  picurv submit --run-dir runs/my_run --stage solve\n"
            "  picurv submit --study-dir studies/my_study --dry-run"
        ),
        epilog="Next: use `picurv summarize --run-dir ...` or `picurv cancel --run-dir ...` after submission as needed.",
    )
    target_group = p_submit.add_mutually_exclusive_group(required=True)
    target_group.add_argument("--run-dir", help="Path to a staged run directory created by `picurv run ... --no-submit`.")
    target_group.add_argument("--study-dir", help="Path to a staged study directory created by `picurv sweep --cluster ... --no-submit`.")
    p_submit.add_argument(
        "--stage",
        choices=["all", "solve", "post-process"],
        default="all",
        help="Which staged job(s) to submit (default: all).",
    )
    p_submit.add_argument(
        "--force",
        action="store_true",
        help="Allow re-submitting a stage already marked as submitted in scheduler/submission.json.",
    )
    p_submit.add_argument(
        "--dry-run",
        action="store_true",
        help="Show which local command(s) or `sbatch` command(s) would run without starting anything.",
    )
    return p_submit


def _add_init_parser(subparsers):
    """!
    @brief Perform add init parser.
    @param[in] subparsers Argument passed to `_add_init_parser()`.
    @return Value returned by `_add_init_parser()`.
    """
    p_init = subparsers.add_parser(
        "init",
        help="Initialize a new case study directory from a template.",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Create a study directory from examples/<template_name>.\n\n"
            "Examples:\n"
            "  picurv init flat_channel --dest my_case\n"
            "  picurv init bent_channel --dest my_bent_case"
        ),
        epilog="Next: run `picurv validate --case ... --solver ... --monitor ...` before execution.",
    )
    p_init.add_argument("template_name", help="Name of the case template directory to copy (e.g., 'flat_channel').")
    p_init.add_argument(
        "--dest",
        dest="dest_name",
        help="Optional name for the new directory. Defaults to the template name.\nPath is relative to your current working directory.",
    )
    p_init.add_argument(
        "--source-root",
        help="Optional override for the PICurv source repository root.\nUseful when running from a copied case without metadata.",
    )
    p_init.add_argument(
        "--pin-binaries",
        action="store_true",
        default=False,
        help="Copy simulator and postprocessor into the case directory.\nUse this to freeze specific binary versions for reproducibility\nor to protect running jobs from concurrent rebuilds.",
    )
    return p_init


def _add_build_parser(subparsers):
    """!
    @brief Perform add build parser.
    @param[in] subparsers Argument passed to `_add_build_parser()`.
    @return Value returned by `_add_build_parser()`.
    """
    p_build = subparsers.add_parser(
        "build",
        help="Build project executables using the Makefile.",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Calls the project's Makefile directly through `make`.\n"
            "If you do not pass an explicit make target, `picurv build` runs `make all`.\n"
            "Any arguments provided after 'build' are passed directly to make.\n\n"
            "Examples:\n"
            "  picurv build\n"
            "  picurv build clean-project\n"
            "  picurv build SYSTEM=cluster\n"
            "  picurv build all SYSTEM=cluster\n"
            "  picurv build postprocessor\n"
            "  ./picurv build clean-project   # from an initialized case directory"
        ),
        epilog="Next: run `picurv --help` or `picurv run --help` for execution commands.",
    )
    p_build.add_argument(
        "--source-root",
        help="Optional override for the PICurv source repository root.",
    )
    p_build.add_argument(
        "--case-dir",
        help="Optional case directory used to resolve .picurv-origin.json when not running from that case.",
    )
    p_build.add_argument(
        "make_args",
        nargs=argparse.REMAINDER,
        help="Arguments to pass directly to the make command (e.g., 'clean-project').",
    )
    return p_build


def _add_sync_binaries_parser(subparsers):
    """!
    @brief Perform add sync binaries parser.
    @param[in] subparsers Argument passed to `_add_sync_binaries_parser()`.
    @return Value returned by `_add_sync_binaries_parser()`.
    """
    p_sync_binaries = subparsers.add_parser(
        "sync-binaries",
        help="Pin specific binary versions into a case directory (optional).",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Copy simulator/postprocessor from the source repo bin/ into a case directory.\n"
            "This is optional — binaries are normally resolved via PATH. Use this only\n"
            "when you need to pin a specific build into a case for reproducibility.\n\n"
            "Examples:\n"
            "  picurv sync-binaries --case-dir my_case\n"
            "  picurv sync-binaries --source-root /path/to/PICurv"
        ),
        epilog="Next: run `./picurv build ...` first if the source repo bin/ is out of date.",
    )
    p_sync_binaries.add_argument("--case-dir", help="Optional case directory to refresh. Defaults to the current case.")
    p_sync_binaries.add_argument("--source-root", help="Optional override for the PICurv source repository root.")
    return p_sync_binaries


def _add_sync_config_parser(subparsers):
    """!
    @brief Perform add sync config parser.
    @param[in] subparsers Argument passed to `_add_sync_config_parser()`.
    @return Value returned by `_add_sync_config_parser()`.
    """
    p_sync_config = subparsers.add_parser(
        "sync-config",
        help="Refresh template-managed files in a case directory from examples/<template>.",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Copy updated example template files into an existing case directory.\n"
            "Modified files are preserved unless --overwrite is used.\n"
            "--prune removes only files that were previously tracked as template-managed and\n"
            "have since been removed from the source template.\n\n"
            "Examples:\n"
            "  ./picurv sync-config\n"
            "  ./picurv sync-config --overwrite\n"
            "  ./picurv sync-config --prune\n"
            "  ./bin/picurv sync-config --case-dir my_case --template-name flat_channel"
        ),
        epilog="Next: run `picurv validate ...` after syncing configs.",
    )
    p_sync_config.add_argument("--case-dir", help="Optional case directory to refresh. Defaults to the current case.")
    p_sync_config.add_argument("--source-root", help="Optional override for the PICurv source repository root.")
    p_sync_config.add_argument(
        "--template-name",
        help="Optional template name override (e.g., flat_channel). Required when metadata is absent.",
    )
    p_sync_config.add_argument("--overwrite", action="store_true", help="Overwrite case files even when they differ from the template.")
    p_sync_config.add_argument(
        "--prune",
        action="store_true",
        help="Remove stale files that were previously tracked as template-managed but no longer exist in the source template.",
    )
    return p_sync_config


def _add_pull_source_parser(subparsers):
    """!
    @brief Perform add pull source parser.
    @param[in] subparsers Argument passed to `_add_pull_source_parser()`.
    @return Value returned by `_add_pull_source_parser()`.
    """
    p_pull = subparsers.add_parser(
        "pull-source",
        help="Refresh source branches from an initialized case directory.",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Update the source repository without leaving an initialized case directory.\n\n"
            "By default this refreshes every local branch that tracks an upstream,\n"
            "then restores the branch you started on.\n\n"
            "Examples:\n"
            "  ./picurv pull-source\n"
            "  ./picurv pull-source --current-branch-only\n"
            "  ./picurv pull-source --no-rebase\n"
            "  ./picurv pull-source --remote origin --branch main"
        ),
        epilog="Next: run `./picurv build` and `./picurv sync-binaries` if source changes require rebuilt executables.",
    )
    p_pull.add_argument("--case-dir", help="Optional case directory used to resolve .picurv-origin.json.")
    p_pull.add_argument("--source-root", help="Optional override for the PICurv source repository root.")
    p_pull.add_argument("--remote", help="Optional git remote name (e.g., origin).")
    p_pull.add_argument("--branch", help="Optional branch name. If provided without --remote, origin is assumed.")
    p_pull.add_argument(
        "--current-branch-only",
        action="store_true",
        help="Only pull the currently checked out branch instead of iterating across all local tracking branches.",
    )
    p_pull.add_argument("--no-rebase", action="store_true", help="Use plain git pull instead of git pull --rebase.")
    return p_pull


def _add_status_source_parser(subparsers):
    """!
    @brief Perform add status source parser.
    @param[in] subparsers Argument passed to `_add_status_source_parser()`.
    @return Value returned by `_add_status_source_parser()`.
    """
    p_status = subparsers.add_parser(
        "status-source",
        help="Report source/case drift for an initialized case directory.",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Inspect whether the source repo, copied binaries, and template-managed files have drifted\n"
            "from the current case directory.\n\n"
            "Examples:\n"
            "  ./picurv status-source\n"
            "  ./picurv status-source --format json\n"
            "  ./bin/picurv status-source --case-dir my_case"
        ),
        epilog="Next: use `pull-source`, `build`, `sync-binaries`, or `sync-config` based on the reported drift.",
    )
    p_status.add_argument("--case-dir", help="Optional case directory used to resolve .picurv-origin.json.")
    p_status.add_argument("--source-root", help="Optional override for the PICurv source repository root.")
    p_status.add_argument("--template-name", help="Optional template name override when metadata is absent.")
    p_status.add_argument(
        "--format",
        dest="output_format",
        choices=["text", "json"],
        default="text",
        help="Output format (default: text).",
    )
    return p_status


def build_main_parser():
    """!
    @brief Build and return the top-level CLI parser.
    @return Value returned by `build_main_parser()`.
    """
    parser = argparse.ArgumentParser(
        description="picurv: A comprehensive conductor for the PICurv simulation platform.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Examples:\n"
            "  picurv validate --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml\n"
            "  picurv run --solve --post-process --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml --dry-run\n"
            "  picurv run --solve --post-process --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml --no-submit\n"
            "  picurv precompute --case case.yml --output-dir precomputed/my_case\n"
            "  picurv run --post-process --continue --run-dir runs/my_run --post post.yml\n"
            "  picurv summarize --run-dir runs/my_run --latest\n"
            "  picurv summarize --run-dir runs/my_run --list-plot-series\n"
            "  picurv summarize --run-dir runs/my_run --plot momentum.residual_norm --last 100\n"
            "  picurv submit --run-dir runs/my_run\n"
            "  picurv cancel --run-dir runs/my_run\n"
            "  picurv cancel --run-dir runs/my_run --stage solve --graceful\n"
            "  picurv sweep --study study.yml --cluster cluster.yml\n\n"
            "Next commands:\n"
            "  - First run: picurv init ... -> picurv validate ... -> picurv run ...\n"
            "  - Config debugging: picurv validate ...\n"
            "  - Artifact generation: picurv precompute ...\n"
            "  - Launch planning: picurv run ... --dry-run\n"
            "  - Post-only catch-up: picurv run --post-process --continue --run-dir ... --post ...\n"
            "  - Deferred submission: picurv submit --run-dir ...\n"
            "  - Run inspection/plots: picurv summarize ...\n"
            "  - Run cancellation: picurv cancel --run-dir ...; use --graceful for solver final-output shutdown"
        ),
    )
    parser.add_argument("-v", "--version", action="version", version=f"picurv {PICURV_VERSION}")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")
    _add_run_parser(subparsers)
    _add_sweep_parser(subparsers)
    _add_validate_parser(subparsers)
    _add_precompute_parser(subparsers)
    _add_summarize_parser(subparsers)
    _add_submit_parser(subparsers)
    _add_cancel_parser(subparsers)
    _add_init_parser(subparsers)
    _add_build_parser(subparsers)
    _add_sync_binaries_parser(subparsers)
    _add_sync_config_parser(subparsers)
    _add_pull_source_parser(subparsers)
    _add_status_source_parser(subparsers)
    return parser


def dispatch_command(args):
    """!
    @brief Validate argument combinations and dispatch to command handlers.
    @param[in] args Command-line style argument list supplied to the function.
    """
    if args.command == "run":
        if not args.solve and not args.post_process:
            fail_cli_usage("At least one stage (--solve or --post-process) must be selected.")
        if args.solve and (not args.case or not args.solver or not args.monitor):
            fail_cli_usage("--solve requires --case, --solver, and --monitor.")
        if args.post_process and not args.post:
            fail_cli_usage("--post-process requires --post.")
        if args.scheduler and not args.cluster:
            fail_cli_usage("--scheduler requires --cluster in this version.")
        run_workflow(args)
        return
    if args.command == "sweep":
        if args.continue_study:
            if not args.study_dir:
                fail_cli_usage("--continue requires --study-dir.")
            sweep_continue_workflow(args)
            return
        if args.reaggregate:
            if not args.study_dir:
                fail_cli_usage("--reaggregate requires --study-dir.")
            sweep_reaggregate_workflow(args)
            return
        # Default: launch new study
        if not args.study:
            fail_cli_usage("--study is required when launching a new sweep.")
        if not args.cluster:
            fail_cli_usage("--cluster is required when launching a new sweep.")
        sweep_workflow(args)
        return
    if args.command == "validate":
        validate_workflow(args)
        return
    if args.command == "precompute":
        try:
            precompute_workflow(args)
        except ValueError as e:
            emit_structured_error(
                ERROR_CODE_CFG_INVALID_VALUE,
                key="precompute",
                file_path=getattr(args, "case", "-"),
                message=str(e),
            )
            sys.exit(1)
        return
    if args.command == "summarize":
        summarize_workflow(args)
        return
    if args.command == "submit":
        submit_staged_jobs(args)
        return
    if args.command == "cancel":
        cancel_run_jobs(args)
        return
    if args.command == "init":
        init_case(args)
        return
    if args.command == "build":
        build_project(args)
        return
    if args.command == "sync-binaries":
        sync_case_binaries_command(args)
        return
    if args.command == "sync-config":
        sync_case_config_command(args)
        return
    if args.command == "pull-source":
        pull_source_repo(args)
        return
    if args.command == "status-source":
        status_source_command(args)
        return
    fail_cli_usage(f"Unsupported command '{args.command}'.")
