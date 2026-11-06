
- Fixed: a spatially varying analytical field set wrong boundary values and dummy cells on `+` faces of a non-periodic axis, because it read a face centre from the dummy layer. Its interior values were right.
- Internal: inlets and analytical fields share `PicurvSourceFaceSlotsCreate()`, the analytical engine shares `FinalizeCellFieldsFromBoundaryValues()`, and the verification scalar functions moved to `verification_sources.c`. A curvilinear analytical boundary now takes the face-metric direction inlets use.
