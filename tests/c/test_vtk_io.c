/**
 * @file test_vtk_io.c
 * @brief C unit tests for VTK output preparation and file writing helpers.
 */

#include "test_support.h"

#include "vtk_io.h"

#include <stdio.h>
#include <string.h>
/**
 * @brief Tests coordinate extraction with interior-grid subsampling.
 */

static PetscErrorCode TestPrepareOutputCoordinatesSubsamplesInteriorGrid(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscScalar *coords = NULL;
    PetscInt nx = 0, ny = 0, nz = 0, npoints = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));

    PetscCall(PrepareOutputCoordinates(user, &coords, &nx, &ny, &nz, &npoints));
    PetscCall(PicurvAssertIntEqual(3, nx, "PrepareOutputCoordinates should output nx=IM-1"));
    PetscCall(PicurvAssertIntEqual(3, ny, "PrepareOutputCoordinates should output ny=JM-1"));
    PetscCall(PicurvAssertIntEqual(3, nz, "PrepareOutputCoordinates should output nz=KM-1"));
    PetscCall(PicurvAssertIntEqual(27, npoints, "PrepareOutputCoordinates should output (IM-1)*(JM-1)*(KM-1) points"));

    if (simCtx->rank == 0) {
        PetscCall(PicurvAssertBool((PetscBool)(coords != NULL), "PrepareOutputCoordinates should allocate coordinates on rank 0"));
        PetscCall(PicurvAssertRealNear(0.0, PetscRealPart(coords[0]), 1.0e-12, "First coordinate x should be 0"));
        PetscCall(PicurvAssertRealNear(0.0, PetscRealPart(coords[1]), 1.0e-12, "First coordinate y should be 0"));
        PetscCall(PicurvAssertRealNear(0.0, PetscRealPart(coords[2]), 1.0e-12, "First coordinate z should be 0"));

        PetscCall(PicurvAssertRealNear(2.0, PetscRealPart(coords[3 * (npoints - 1) + 0]), 1.0e-12, "Last subsampled coordinate x should be IM-2"));
        PetscCall(PicurvAssertRealNear(2.0, PetscRealPart(coords[3 * (npoints - 1) + 1]), 1.0e-12, "Last subsampled coordinate y should be JM-2"));
        PetscCall(PicurvAssertRealNear(2.0, PetscRealPart(coords[3 * (npoints - 1) + 2]), 1.0e-12, "Last subsampled coordinate z should be KM-2"));

        PetscCall(PetscFree(coords));
    }

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests Eulerian scalar field extraction with subsampling.
 */

static PetscErrorCode TestPrepareOutputEulerianFieldDataSubsamplesScalar(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscScalar *field_out = NULL;
    PetscInt start = 0, end = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));

    PetscCall(VecGetOwnershipRange(user->P_nodal, &start, &end));
    for (PetscInt idx = start; idx < end; ++idx) {
        PetscCall(VecSetValue(user->P_nodal, idx, (PetscScalar)idx, INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(user->P_nodal));
    PetscCall(VecAssemblyEnd(user->P_nodal));

    PetscCall(PrepareOutputEulerianFieldData(user, user->P_nodal, 1, &field_out));
    if (simCtx->rank == 0) {
        PetscCall(PicurvAssertBool((PetscBool)(field_out != NULL),
                                   "PrepareOutputEulerianFieldData should allocate subsampled data on rank 0"));
        PetscCall(PicurvAssertRealNear(0.0, PetscRealPart(field_out[0]), 1.0e-12, "Subsampled first value should map to source index 0"));
        PetscCall(PicurvAssertRealNear(42.0, PetscRealPart(field_out[26]), 1.0e-12,
                                       "Subsampled last value should map to source index (2,2,2)=42 in 4x4x4"));
        PetscCall(PetscFree(field_out));
    }

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests particle-data extraction with output subsampling.
 */

static PetscErrorCode TestPrepareOutputParticleDataSubsampling(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    VTKMetaData meta;
    PetscInt n_total = 0;
    PetscReal (*position)[3] = NULL;
    PetscReal (*velocity)[3] = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PetscMemzero(&meta, sizeof(meta)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 5, "ske"));

    pps.particle_output_freq = 2;
    PetscCall(PetscStrncpy(pps.particle_fields, "position,velocity", sizeof(pps.particle_fields)));

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&position));
    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void *)&velocity));
    for (PetscInt p = 0; p < 5; ++p) {
        position[p][0] = (PetscReal)p;
        position[p][1] = 0.0;
        position[p][2] = 0.0;
        velocity[p][0] = 10.0 + (PetscReal)p;
        velocity[p][1] = 20.0 + (PetscReal)p;
        velocity[p][2] = 30.0 + (PetscReal)p;
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&position));
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void *)&velocity));

    PetscCall(PrepareOutputParticleData(user, &pps, &meta, &n_total));
    if (simCtx->rank == 0) {
        PetscCall(PicurvAssertIntEqual(5, n_total, "PrepareOutputParticleData should report total particles before subsampling"));
        PetscCall(PicurvAssertIntEqual(3, meta.npoints, "Stride-2 subsampling over 5 particles should keep 3 output particles"));
        PetscCall(PicurvAssertIntEqual(VTK_POLYDATA, meta.fileType, "Particle output metadata should be PolyData"));
        PetscCall(PicurvAssertIntEqual(1, meta.num_point_data_fields, "Velocity field should be exported as one point-data field"));
        PetscCall(PicurvAssertBool((PetscBool)(strcmp(meta.point_data_fields[0].name, "velocity") == 0),
                                   "Particle field name should propagate to VTK metadata"));

        PetscCall(PicurvAssertRealNear(0.0, PetscRealPart(meta.coords[0]), 1.0e-12, "First exported particle should be index 0"));
        PetscCall(PicurvAssertRealNear(2.0, PetscRealPart(meta.coords[3]), 1.0e-12, "Second exported particle should be index 2"));
        PetscCall(PicurvAssertRealNear(4.0, PetscRealPart(meta.coords[6]), 1.0e-12, "Third exported particle should be index 4"));

        PetscCall(PetscFree(meta.coords));
        for (PetscInt f = 0; f < meta.num_point_data_fields; ++f) {
            PetscCall(PetscFree(meta.point_data_fields[f].data));
        }
        PetscCall(PetscFree(meta.connectivity));
        PetscCall(PetscFree(meta.offsets));
    }

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests structured-grid VTK file creation from assembled metadata.
 */

static PetscErrorCode TestCreateVTKFileFromMetadataWritesStructuredGrid(void)
{
    char tmpdir[PETSC_MAX_PATH_LEN];
    char vtk_path[PETSC_MAX_PATH_LEN];
    VTKMetaData meta;
    FILE *file = NULL;
    char probe[512];
    size_t nread = 0;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&meta, sizeof(meta)));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(vtk_path, sizeof(vtk_path), "%s/field_00001.vts", tmpdir));

    meta.fileType = VTK_STRUCTURED;
    meta.mx = 2;
    meta.my = 2;
    meta.mz = 2;
    meta.npoints = 8;
    meta.num_point_data_fields = 1;

    PetscCall(PetscMalloc1(3 * meta.npoints, &meta.coords));
    PetscCall(PetscMalloc1(meta.npoints, &meta.point_data_fields[0].data));
    PetscCall(PetscStrncpy(meta.point_data_fields[0].name, "P_nodal", sizeof(meta.point_data_fields[0].name)));
    meta.point_data_fields[0].num_components = 1;

    for (PetscInt p = 0; p < meta.npoints; ++p) {
        meta.coords[3 * p + 0] = (PetscScalar)(p % 2);
        meta.coords[3 * p + 1] = (PetscScalar)((p / 2) % 2);
        meta.coords[3 * p + 2] = (PetscScalar)(p / 4);
        meta.point_data_fields[0].data[p] = (PetscScalar)p;
    }

    PetscCall(CreateVTKFileFromMetadata(vtk_path, &meta, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertFileExists(vtk_path, "CreateVTKFileFromMetadata should emit a .vts file"));

    file = fopen(vtk_path, "rb");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open generated VTK file '%s'.", vtk_path);
    nread = fread(probe, 1, sizeof(probe) - 1, file);
    fclose(file);
    probe[nread] = '\0';
    PetscCall(PicurvAssertBool((PetscBool)(strstr(probe, "StructuredGrid") != NULL),
                               "Generated VTK header should declare StructuredGrid type"));

    PetscCall(PetscFree(meta.coords));
    PetscCall(PetscFree(meta.point_data_fields[0].data));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-post-vtk PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"prepare-output-coordinates-subsamples-interior-grid", TestPrepareOutputCoordinatesSubsamplesInteriorGrid},
        {"prepare-output-eulerian-field-data-subsamples-scalar", TestPrepareOutputEulerianFieldDataSubsamplesScalar},
        {"prepare-output-particle-data-subsampling", TestPrepareOutputParticleDataSubsampling},
        {"create-vtk-file-from-metadata-writes-structured-grid", TestCreateVTKFileFromMetadataWritesStructuredGrid},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv VTK I/O tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-post-vtk", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
