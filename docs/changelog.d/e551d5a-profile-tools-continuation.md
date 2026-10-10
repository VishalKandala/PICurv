- The driven-periodic profile tools (`wall_normal_profile.py`, `cross_section_profile.py`)
  read `driven_flow.csv` and `wall_model.csv` from continued runs. They skip the
  continuation marker and count a repeated step once, instead of crashing.
