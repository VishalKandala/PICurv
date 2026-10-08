- `examples/periodic_test/driven_channel/tools/wall_normal_profile.py` works again on
  current checkpoints. It looks the statistics window up in `checkpoint.meta`, reads the
  six-component `Ucat_m2` payload correctly, reports the actual Reynolds shear stress,
  and folds the two walls. A new `--wall-model-csv` option takes u_tau from the wall
  model, and `--u-tau` supplies it directly. The tool does not apply to `driven_duct`,
  which has only one homogeneous direction, and that example's comment now says so.
