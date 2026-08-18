/**
 * @file checksum.c
 * @brief Dependency-free SHA-256 implementation used for small persistent metadata.
 */

#include "checksum.h"

#include <stdio.h>

#define ROTATE_RIGHT(value, bits) (((value) >> (bits)) | ((value) << (32u - (bits))))

static const uint32_t kRoundConstants[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
    0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
    0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
    0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
    0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
    0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
};

/** @brief Transform one complete 64-byte SHA-256 block. */
static void PicurvSHA256Transform(PicurvSHA256Context *context, const unsigned char block[64])
{
    uint32_t words[64];
    uint32_t a, b, c, d, e, f, g, h;

    for (size_t index = 0; index < 16; ++index) {
        const size_t offset = 4 * index;
        words[index] = ((uint32_t)block[offset] << 24) |
                       ((uint32_t)block[offset + 1] << 16) |
                       ((uint32_t)block[offset + 2] << 8) |
                       (uint32_t)block[offset + 3];
    }
    for (size_t index = 16; index < 64; ++index) {
        const uint32_t s0 = ROTATE_RIGHT(words[index - 15], 7) ^ ROTATE_RIGHT(words[index - 15], 18) ^ (words[index - 15] >> 3);
        const uint32_t s1 = ROTATE_RIGHT(words[index - 2], 17) ^ ROTATE_RIGHT(words[index - 2], 19) ^ (words[index - 2] >> 10);
        words[index] = words[index - 16] + s0 + words[index - 7] + s1;
    }

    a = context->state[0]; b = context->state[1]; c = context->state[2]; d = context->state[3];
    e = context->state[4]; f = context->state[5]; g = context->state[6]; h = context->state[7];

    for (size_t index = 0; index < 64; ++index) {
        const uint32_t sum1 = ROTATE_RIGHT(e, 6) ^ ROTATE_RIGHT(e, 11) ^ ROTATE_RIGHT(e, 25);
        const uint32_t choose = (e & f) ^ ((~e) & g);
        const uint32_t temp1 = h + sum1 + choose + kRoundConstants[index] + words[index];
        const uint32_t sum0 = ROTATE_RIGHT(a, 2) ^ ROTATE_RIGHT(a, 13) ^ ROTATE_RIGHT(a, 22);
        const uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
        const uint32_t temp2 = sum0 + majority;

        h = g; g = f; f = e; e = d + temp1;
        d = c; c = b; b = a; a = temp1 + temp2;
    }

    context->state[0] += a; context->state[1] += b; context->state[2] += c; context->state[3] += d;
    context->state[4] += e; context->state[5] += f; context->state[6] += g; context->state[7] += h;
}

/** @brief Implementation of \ref PicurvSHA256Init(). */
void PicurvSHA256Init(PicurvSHA256Context *context)
{
    if (!context) return;
    context->state[0] = 0x6a09e667u; context->state[1] = 0xbb67ae85u;
    context->state[2] = 0x3c6ef372u; context->state[3] = 0xa54ff53au;
    context->state[4] = 0x510e527fu; context->state[5] = 0x9b05688cu;
    context->state[6] = 0x1f83d9abu; context->state[7] = 0x5be0cd19u;
    context->bit_count = 0;
    context->buffer_size = 0;
}

/** @brief Implementation of \ref PicurvSHA256Update(). */
void PicurvSHA256Update(PicurvSHA256Context *context, const void *data, size_t length)
{
    const unsigned char *bytes = (const unsigned char *)data;

    if (!context || (!bytes && length > 0)) return;
    context->bit_count += (uint64_t)length * 8u;

    while (length > 0) {
        const size_t room = 64 - context->buffer_size;
        const size_t take = length < room ? length : room;
        for (size_t index = 0; index < take; ++index) {
            context->buffer[context->buffer_size + index] = bytes[index];
        }
        context->buffer_size += take;
        bytes += take;
        length -= take;
        if (context->buffer_size == 64) {
            PicurvSHA256Transform(context, context->buffer);
            context->buffer_size = 0;
        }
    }
}

/** @brief Implementation of \ref PicurvSHA256FinalHex(). */
void PicurvSHA256FinalHex(PicurvSHA256Context *context, char digest_hex[65])
{
    static const char hex[] = "0123456789abcdef";
    unsigned char digest[32];
    const uint64_t message_bits = context ? context->bit_count : 0;
    unsigned char padding[128] = {0x80};
    size_t padding_length;

    if (!context || !digest_hex) return;
    padding_length = context->buffer_size < 56 ? 56 - context->buffer_size : 120 - context->buffer_size;
    PicurvSHA256Update(context, padding, padding_length);
    context->bit_count = message_bits;

    for (size_t index = 0; index < 8; ++index) {
        context->buffer[56 + index] = (unsigned char)(message_bits >> (56 - 8 * index));
    }
    PicurvSHA256Transform(context, context->buffer);

    for (size_t index = 0; index < 8; ++index) {
        digest[4 * index] = (unsigned char)(context->state[index] >> 24);
        digest[4 * index + 1] = (unsigned char)(context->state[index] >> 16);
        digest[4 * index + 2] = (unsigned char)(context->state[index] >> 8);
        digest[4 * index + 3] = (unsigned char)context->state[index];
    }
    for (size_t index = 0; index < sizeof(digest); ++index) {
        digest_hex[2 * index] = hex[digest[index] >> 4];
        digest_hex[2 * index + 1] = hex[digest[index] & 0x0f];
    }
    digest_hex[64] = '\0';
}

/** @brief Implementation of \ref PicurvSHA256File(). */
PetscErrorCode PicurvSHA256File(const char *path, char digest_hex[65])
{
    FILE *file = NULL;
    unsigned char buffer[8192];
    PicurvSHA256Context context;

    PetscFunctionBeginUser;
    PetscCheck(path != NULL && digest_hex != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "SHA-256 file path and digest output are required.");
    file = fopen(path, "rb");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
               "Unable to open '%s' while computing SHA-256.", path);

    PicurvSHA256Init(&context);
    while (!feof(file)) {
        const size_t count = fread(buffer, 1, sizeof(buffer), file);
        if (count > 0) PicurvSHA256Update(&context, buffer, count);
        PetscCheck(!ferror(file), PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
                   "Unable to read '%s' while computing SHA-256.", path);
    }
    PetscCheck(fclose(file) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
               "Unable to close '%s' after computing SHA-256.", path);
    PicurvSHA256FinalHex(&context, digest_hex);
    PetscFunctionReturn(0);
}

#undef ROTATE_RIGHT
