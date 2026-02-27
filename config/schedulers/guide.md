# Scheduler Config Guide

This directory stores reusable scheduler profiles for `pic.flow` cluster runs and sweeps.

Use with:
- `./scripts/pic.flow run --cluster config/schedulers/slurm_default.yml ...`
- `./scripts/pic.flow sweep --cluster config/schedulers/slurm_default.yml ...`

Reference schema:
- `examples/master_template/master_cluster.yml`
