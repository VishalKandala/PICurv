"""!
@file petsc_binary.py
@brief Reader for the PETSc binary Vec files PICurv checkpoints hold, shared by the generators.

Not a command. `profile.gen` and `spectra.gen` import it from beside themselves, as
`spectra.gen` and `ic.gen` import `periodic_spectral.py`.
"""

import os
import struct

#: PETSc's class id for a Vec, the first header word of a VecView binary file.
PETSC_VEC_CLASS_ID = 1211214


def read_petsc_vec_binary(path: str, numpy, index_bytes: int = 4, label: str = "field file"):
    """!
    @brief Read a real64 PETSc binary VecView payload.
    @details The header's class id and count are PetscInt of `index_bytes`, the width of the
             PETSc that wrote the file, which the file does not record. Either byte order is
             accepted, with or without the class id; the file length must match exactly.
    @param[in] path PETSc binary vector path.
    @param[in] numpy The caller's numpy module.
    @param[in] index_bytes Bytes of a PetscInt in the writing PETSc, 4 or 8.
    @param[in] label How the caller names the file in errors.
    @return NumPy array of scalar values.
    """
    if not os.path.isfile(path):
        raise ValueError(f"{label} not found: {path}")
    data = open(path, "rb").read()
    word, code = index_bytes, "i" if index_bytes == 4 else "q"
    if len(data) < 2 * word:
        raise ValueError(f"{label} '{path}' is too short to be a PETSc Vec binary.")
    candidates = []
    for endian in (">", "<"):
        _class_id, vec_size = struct.unpack(endian + 2 * code, data[:2 * word])
        if vec_size > 0 and len(data) == 2 * word + vec_size * 8:
            candidates.append((endian, vec_size, 2 * word))
        (vec_size_no_class,) = struct.unpack(endian + code, data[:word])
        if vec_size_no_class > 0 and len(data) == word + vec_size_no_class * 8:
            candidates.append((endian, vec_size_no_class, word))
    if not candidates:
        raise ValueError(
            f"{label} '{path}' does not match a PETSc real64 Vec binary with "
            f"{index_bytes}-byte indices; a file from the other index width does not fit."
        )
    endian, vec_size, offset = candidates[0]
    values = numpy.frombuffer(data, dtype=numpy.dtype(endian + "f8"), count=vec_size, offset=offset)
    return values.astype(float, copy=False)
