"""!
@file petsc_binary.py
@brief The PETSc binary Vec envelope of PICurv's checkpoints and staged fields: its reader,
       writer and validator, shared by the conductor and the generators.

Not a command. `core.py`, `ic.gen`, `profile.gen` and `spectra.gen` import it from `picurv_cli/`,
as they import `periodic_spectral.py`. The envelope is a class id and a scalar count, each a
PetscInt of the writing build's width, then real64 values, all big-endian.
"""

import os
import struct

#: PETSc's class id for a Vec, the first header word of a VecView binary file.
PETSC_VEC_CLASS_ID = 1211214


def _index_code(index_bytes: int) -> str:
    """!
    @brief The struct code of a PetscInt of one width.
    @param[in] index_bytes Bytes of a PetscInt, 4 or 8.
    @return `i` or `q`.
    """
    return "i" if index_bytes == 4 else "q"


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
    with open(path, "rb") as fin:
        data = fin.read()
    word, code = index_bytes, _index_code(index_bytes)
    if len(data) < 2 * word:
        raise ValueError(f"{label} '{path}' is too short to be a PETSc Vec binary.")
    candidates = []
    for endian in (">", "<"):
        (vec_size,) = struct.unpack_from(endian + code, data, word)
        if vec_size > 0 and len(data) == 2 * word + vec_size * 8:
            candidates.append((endian, vec_size, 2 * word))
        (vec_size_no_class,) = struct.unpack_from(endian + code, data, 0)
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


def write_petsc_vec_binary(path: str, payload: bytes, index_bytes: int = 4) -> None:
    """!
    @brief Write a PETSc binary Vec: the class id and the count, then the values.
    @details The class id and count are PetscInt of `index_bytes`, so the reading solver's PETSc,
             which does not detect the width, reads them as it wrote its own.
    @param[in] path Output path.
    @param[in] payload Values as big-endian real64 bytes.
    @param[in] index_bytes Bytes of a PetscInt in the reading PETSc, 4 or 8.
    """
    with open(path, "wb") as fout:
        fout.write(struct.pack(">" + 2 * _index_code(index_bytes), PETSC_VEC_CLASS_ID, len(payload) // 8))
        fout.write(payload)


def validate_petsc_vec_binary(path: str, index_bytes: int) -> dict:
    """!
    @brief Validate the envelope the solver's ReadFieldData loads: big-endian, with the class id.
    @details Stricter than read_petsc_vec_binary(), which also takes either byte order, with or
             without the class id: this guards a file the solver itself will read.
    @param[in] path PETSc binary vector path.
    @param[in] index_bytes PetscInt width of the reading build.
    @return Summary containing the absolute path, scalar count and index width.
    """
    header_bytes = 2 * index_bytes
    with open(path, "rb") as fin:
        header = fin.read(header_bytes)
        if len(header) != header_bytes:
            raise ValueError(f"PETSc Vec file is too short: {path}")
        class_id, scalar_count = struct.unpack(">" + 2 * _index_code(index_bytes), header)
        payload_bytes = len(fin.read())
    if class_id != PETSC_VEC_CLASS_ID or scalar_count < 0 or payload_bytes != scalar_count * 8:
        raise ValueError(
            f"Invalid PETSc Vec file {path} for a PETSc with {index_bytes}-byte indices; a file "
            f"written for the other index width (or truncated) does not fit its header.")
    return {"path": os.path.abspath(path), "scalar_count": scalar_count, "index_bytes": index_bytes}
