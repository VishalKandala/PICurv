- The Werner–Wengle wall model (`wall_function.model: werner`) now takes the wall stress
  from the pointwise profile at the reference point. It previously applied the
  cell-integrated relation to a point velocity, which over-predicted the friction
  velocity by 12–17% depending on the reference speed; on the wall-modelled
  Re_tau ≈ 1000 channel this put u_tau 13.5% above Lee & Moser (2015). Runs that used
  this model will change. The corrected model has not yet been compared against a
  reference flow; the predicted −1% for that channel is an estimate from its mean
  profile, not a rerun. A unit-test assertion that let NaN pass is also fixed.
