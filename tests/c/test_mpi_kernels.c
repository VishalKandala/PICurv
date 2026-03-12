/**
 * @file test_mpi_kernels.c
 * @brief C unit tests for MPI collective helper routines.
 */

#include "test_support.h"

#include "ParticleSwarm.h"
#include "grid.h"

#include <stdlib.h>
/**
 * @brief Tests collective particle distribution consistency across MPI ranks.
 */

static PetscErrorCode TestDistributeParticlesCollectiveConsistency(void)
{
    PetscMPIInt rank = 0, size = 1;
    PetscInt local_particles = 0;
    PetscInt remainder = 0;
    PetscInt global_particles = 0;
    PetscInt remainder_min = 0, remainder_max = 0;
    const PetscInt total_particles = 137;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(PicurvAssertBool((PetscBool)(size >= 2), "unit-mpi requires at least two MPI ranks"));

    PetscCall(DistributeParticles(total_particles, rank, size, &local_particles, &remainder));

    PetscCall(PicurvAssertIntEqual(
        total_particles / size + (((PetscInt)rank < remainder) ? 1 : 0),
        local_particles,
        "local particle share should match quotient+remainder policy"));

    PetscCallMPI(MPI_Allreduce(&local_particles, &global_particles, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertIntEqual(total_particles, global_particles, "distributed particle count must conserve total particles"));

    PetscCallMPI(MPI_Allreduce(&remainder, &remainder_min, 1, MPIU_INT, MPI_MIN, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(&remainder, &remainder_max, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertIntEqual(remainder_min, remainder_max, "all ranks should report the same remainder"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests multi-rank bounding-box gather and broadcast helpers.
 */

static PetscErrorCode TestBoundingBoxCollectivesMultiRank(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BoundingBox local_bbox;
    BoundingBox *boxes = NULL;
    PetscMPIInt rank = 0, size = 1;
    PetscReal global_min_x = 0.0;
    PetscReal global_max_x = 0.0;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(PicurvAssertBool((PetscBool)(size >= 2), "unit-mpi requires at least two MPI ranks"));

    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 6, 4));
    PetscCall(ComputeLocalBoundingBox(user, &local_bbox));
    PetscCall(GatherAllBoundingBoxes(user, &boxes));
    PetscCall(BroadcastAllBoundingBoxes(user, &boxes));
    PetscCall(PicurvAssertBool((PetscBool)(boxes != NULL), "all ranks should hold the gathered bounding-box table"));

    for (PetscMPIInt r = 0; r < size; ++r) {
        PetscCall(PicurvAssertBool((PetscBool)(boxes[r].min_coords.x <= boxes[r].max_coords.x), "bbox x-range should be valid"));
        PetscCall(PicurvAssertBool((PetscBool)(boxes[r].min_coords.y <= boxes[r].max_coords.y), "bbox y-range should be valid"));
        PetscCall(PicurvAssertBool((PetscBool)(boxes[r].min_coords.z <= boxes[r].max_coords.z), "bbox z-range should be valid"));
    }

    PetscCall(PicurvAssertRealNear(local_bbox.min_coords.x, boxes[rank].min_coords.x, 1.0e-10, "local bbox min x preserved"));
    PetscCall(PicurvAssertRealNear(local_bbox.max_coords.x, boxes[rank].max_coords.x, 1.0e-10, "local bbox max x preserved"));

    PetscCallMPI(MPI_Allreduce(&local_bbox.min_coords.x, &global_min_x, 1, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(&local_bbox.max_coords.x, &global_max_x, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertBool((PetscBool)(global_min_x <= 0.0), "global min x should include domain start"));
    PetscCall(PicurvAssertBool((PetscBool)(global_max_x >= 7.0), "global max x should include domain end"));

    free(boxes);
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-mpi PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"distribute-particles-collective-consistency", TestDistributeParticlesCollectiveConsistency},
        {"bounding-box-collectives-multi-rank", TestBoundingBoxCollectivesMultiRank},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv MPI-focused runtime tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-mpi", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
