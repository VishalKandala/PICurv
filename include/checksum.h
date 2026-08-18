/**
 * @file checksum.h
 * @brief Small dependency-free SHA-256 utility for persistent metadata identity.
 */

#ifndef PICURV_CHECKSUM_H
#define PICURV_CHECKSUM_H

#include <petscsys.h>
#include <stddef.h>
#include <stdint.h>

/** @brief Incremental SHA-256 state. */
typedef struct {
    uint32_t state[8];
    uint64_t bit_count;
    unsigned char buffer[64];
    size_t buffer_size;
} PicurvSHA256Context;

/**
 * @brief Initialize an incremental SHA-256 calculation.
 * @param[out] context Hash state to initialize.
 */
void PicurvSHA256Init(PicurvSHA256Context *context);

/**
 * @brief Add bytes to an incremental SHA-256 calculation.
 * @param[in,out] context Active hash state.
 * @param[in] data Bytes to add.
 * @param[in] length Number of bytes to add.
 */
void PicurvSHA256Update(PicurvSHA256Context *context, const void *data, size_t length);

/**
 * @brief Finish a SHA-256 calculation and return a lowercase hexadecimal digest.
 * @param[in,out] context Active hash state.
 * @param[out] digest_hex Null-terminated 64-character lowercase digest.
 */
void PicurvSHA256FinalHex(PicurvSHA256Context *context, char digest_hex[65]);

/**
 * @brief Compute the lowercase SHA-256 digest of a file.
 * @param[in] path File to hash.
 * @param[out] digest_hex Null-terminated 64-character lowercase digest.
 * @return Zero on success or a PETSc filesystem error.
 */
PetscErrorCode PicurvSHA256File(const char *path, char digest_hex[65]);

#endif /* PICURV_CHECKSUM_H */
