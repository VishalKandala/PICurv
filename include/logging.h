/**
 * @file logging.h
 * @brief Logging utilities and macros for PETSc-based applications.
 *
 * This header defines logging levels, scopes, and macros for consistent logging throughout the application.
 * It provides functions to retrieve the current logging level and macros to simplify logging with scope control.
 */

#ifndef LOGGING_H
#define LOGGING_H

// Include necessary headers
#include <petsc.h>   // PETSc library header
#include <stdlib.h>
#include <string.h>
#include <petscsys.h>
#include <ctype.h>
#include "variables.h"
#include "statistics_window.h"
#include "Boundaries.h"
// --------------------- Logging Levels Definition ---------------------

/**
 * @brief Enumeration of logging levels.
 *
 * Defines various severity levels for logging messages.
 */
typedef enum {
    LOG_ERROR = 0,   /**< Critical errors that may halt the program */
    LOG_WARNING,     /**< Non-critical issues that warrant attention */
    LOG_INFO,        /**< Informational messages about program execution */
    LOG_DEBUG,       /**< Detailed debugging information */
    LOG_TRACE,        /**< Very fine-grained tracing information for in-depth debugging */
    LOG_VERBOSE      /**< Extremely detailed logs, typically for development use only */
} LogLevel;

// -------------------- Logging Scope Definitions ------------------

/**
 * @brief Logging scope definitions for controlling message output.
 *
 * - LOCAL: Logs on the current process using MPI_COMM_SELF.
 * - GLOBAL: Logs across all processes using MPI_COMM_WORLD.
 */
#define LOCAL  0  ///< Scope for local logging on the current process.
#define GLOBAL 1  ///< Scope for global logging across all processes.

//----------------------- Custom KSP Monitor Struct ------------

/**
 * @brief Context for a dual-purpose KSP monitor.
 *
 * This struct holds a file viewer for unconditional logging and a boolean
 * flag to enable/disable optional logging to the console.
 */
typedef struct {
    FILE        *file_handle;        // C file handling for logging.
    PetscBool   log_to_console;   // Flag to enable console output.
    PetscReal   bnorm;            // Stores the norm of the initial RHS vector.
    PetscInt    step;             // Timestep
    PetscInt    block_id;         // the ID of the block this monitor is for.   
    Vec         residual_work[2]; // KSPBuildResidual() workspaces, created on first use and kept with the KSP.
} DualMonitorCtx;

// --------------------- Logging Macros ---------------------

/**
 * @brief Logging macro for PETSc-based applications with scope control.
 *
 * This macro provides a convenient way to log messages with different scopes
 * (LOCAL or GLOBAL) and severity levels. It utilizes PETSc's `PetscPrintf` 
 * function for message output.
 *
 * @param scope Specifies the logging scope:
 *              - LOCAL: Logs on the current process using MPI_COMM_SELF.
 *              - GLOBAL: Logs on all processes using MPI_COMM_WORLD.
 * @param level The severity level of the message (e.g., LOG_INFO, LOG_ERROR).
 * @param fmt   The format string for the message (similar to printf).
 * @param ...   Additional arguments for the format string (optional).
 *
 * Example usage:
 *     LOG(LOCAL, LOG_ERROR, "An error occurred at index %ld.\n", idx);
 *     LOG(GLOBAL, LOG_INFO, "Grid size: %ld x %ld x %ld.\n", nx, ny, nz);
 */
#define LOG(scope, level, fmt, ...) \
    do { \
        /* Determine the MPI communicator based on the scope */ \
        MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
        /* Check if the log level is within the allowed range */ \
        if ((int)(level) <= (int)get_log_level()) { \
            /* Print the message to the specified communicator */ \
            PetscPrintf(comm, fmt, ##__VA_ARGS__); \
        } \
    } while (0)

/**
 * @brief Default logging macro for PETSc-based applications.
 *
 * This macro simplifies logging by defaulting the scope to GLOBAL
 * (i.e., `MPI_COMM_WORLD`) and providing a convenient interface for
 * common logging needs.
 *
 * @param level The severity level of the message (e.g., LOG_ERROR, LOG_INFO).
 * @param fmt   The format string for the log message (similar to printf).
 * @param ...   Additional arguments for the format string (optional).
 *
 * Example usage:
 *     LOG_DEFAULT(LOG_ERROR, "Error occurred at index %ld.\n", idx);
 *     LOG_DEFAULT(LOG_INFO, "Grid size: %ld x %ld x %ld.\n", nx, ny, nz);
 *
 * @note
 * - By default, this macro logs across all MPI processes using `MPI_COMM_WORLD`.
 * - If finer control (e.g., local logging) is required, use the more general `LOG` macro.
 * - The log level is filtered based on the value returned by `get_log_level()`.
 */
#define LOG_DEFAULT(level, fmt, ...) \
    do { \
        /* Set the communicator to global (MPI_COMM_WORLD) by default */ \
        MPI_Comm comm = MPI_COMM_WORLD; \
        /* Check if the log level is within the allowed range */ \
        if ((int)(level) <= (int)get_log_level()) { \
            /* Print the message using PetscPrintf with the global communicator */ \
            PetscPrintf(comm, fmt, ##__VA_ARGS__); \
        } \
    } while (0)

/**
 * @brief Logging macro for PETSc-based applications with scope control, 
 *        using synchronized output across processes.
 *
 * This macro uses `PetscSynchronizedPrintf` and `PetscSynchronizedFlush` to 
 * ensure messages from different ranks are printed in a synchronized (rank-by-rank) 
 * manner, preventing interleaved outputs.
 *
 * @param scope Specifies the logging scope:
 *              - LOCAL:  Logs on the current process using MPI_COMM_SELF.
 *              - GLOBAL: Logs on all processes using MPI_COMM_WORLD.
 * @param level The severity level of the message (e.g., LOG_INFO, LOG_ERROR).
 * @param fmt   The format string for the message (similar to printf).
 * @param ...   Additional arguments for the format string (optional).
 *
 * Example usage:
 *     LOG_SYNC(LOCAL, LOG_ERROR, "An error occurred at index %ld.\n", idx);
 *     LOG_SYNC(GLOBAL, LOG_INFO, "Synchronized info: rank = %ld.\n", rank);
 */
#define LOG_SYNC(scope, level, fmt, ...) \
    do { \
        /* Determine the MPI communicator based on the scope */ \
        MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
        /* Check if the log level is within the allowed range */ \
        if ((int)(level) <= (int)get_log_level()) { \
            /* Synchronized print (collective) on the specified communicator */ \
            PetscSynchronizedPrintf(comm, fmt, ##__VA_ARGS__); \
            /* Ensure all ranks have finished printing before continuing */ \
            PetscSynchronizedFlush(comm, PETSC_STDOUT); \
        } \
    } while (0)

/**
 * @brief Default synchronized logging macro for PETSc-based applications.
 *
 * This macro simplifies logging by defaulting the scope to GLOBAL 
 * (i.e., `MPI_COMM_WORLD`) and provides synchronized output across 
 * all processes.
 *
 * @param level The severity level of the message (e.g., LOG_ERROR, LOG_INFO).
 * @param fmt   The format string for the log message (similar to printf).
 * @param ...   Additional arguments for the format string (optional).
 *
 * Example usage:
 *     LOG_SYNC_DEFAULT(LOG_ERROR, "Error at index %ld.\n", idx);
 *     LOG_SYNC_DEFAULT(LOG_INFO,  "Process rank: %ld.\n", rank);
 *
 * @note
 * - By default, this macro logs across all MPI processes using `MPI_COMM_WORLD`.
 * - If local (per-process) logging is required, use the more general `LOG_SYNC` macro.
 * - The log level is filtered based on the value returned by `get_log_level()`.
 */
#define LOG_SYNC_DEFAULT(level, fmt, ...) \
    do { \
        if ((int)(level) <= (int)get_log_level()) { \
            PetscSynchronizedPrintf(MPI_COMM_WORLD, fmt, ##__VA_ARGS__); \
            PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT); \
        } \
    } while (0)



/**
 * @brief True when a code-level diagnostic at `level` prints from the calling function.
 *
 * The level must be admitted by the run's verbosity. A warning or an error then prints from
 * any function, since an operator running at WARNING must see it; INFO and below print only
 * from allow-listed functions, so instrumenting one function does not drag in the others.
 * The value is the same on every rank.
 */
#define LOG_LEVEL_ADMITS(level) \
    ((int)(level) <= (int)get_log_level() && ((int)(level) <= (int)LOG_WARNING || is_function_allowed(__func__)))

/**
 * @brief Logging macro that checks both the log level and whether the calling function
 *        is in the allowed-function list before printing. Useful for selective, per-function logging.
 *        Warnings and errors ignore the allow-list (see LOG_LEVEL_ADMITS()).
 *
 * @param scope Specifies the logging scope (LOCAL or GLOBAL).
 * @param level The severity level of the message (e.g., LOG_INFO, LOG_ERROR).
 * @param fmt   The format string for the message (similar to printf).
 * @param ...   Additional arguments for the format string (optional).
 *
 * Example usage:
 *     LOG_ALLOW(LOCAL, LOG_DEBUG, "Debugging info in function: %s\n", __func__);
 */
#define LOG_ALLOW(scope, level, fmt, ...) \
    do { \
        MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
        if (LOG_LEVEL_ADMITS(level)) { \
            PetscPrintf(comm, "[%s] " fmt, __func__, ##__VA_ARGS__); \
        } \
    } while (0)


/* ------- DEBUG ------------------------------------------
#define LOG_ALLOW(scope, level, fmt, ...) \
    do { \
        MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
        PetscInt __current_level_val = get_log_level(); \
        PetscBool __allowed_func_val = is_function_allowed(__func__); \
        // Print BEFORE the check  \
        if (strcmp(__func__, "LocateAllParticlesInGrid") == 0) { \
             printf("[DEBUG LOG_ALLOW in %s] Checking: level=%d, get_log_level() returned %d, func_allowed=%d\n", \
                    __func__, (int)level, (int)__current_level_val, (int)__allowed_func_val); \
        } \
        if ((int)(level) <= (int)__current_level_val && __allowed_func_val) { \
             // Print AFTER passing the check // \
             if (strcmp(__func__, "LocateAllParticlesInGrid") == 0) { \
                  printf("[DEBUG LOG_ALLOW in %s] Check PASSED. Printing log.\n", __func__); \
             } \
             PetscPrintf(comm, "[%s] " fmt, __func__, ##__VA_ARGS__); \
        } \
    } while (0)
-------------------------------------------------------------------------------
*/

/**
 * @brief Synchronized logging macro that checks both the log level 
 *        and whether the calling function is in the allow-list.
 *
 * This macro uses `PetscSynchronizedPrintf` and `PetscSynchronizedFlush` to
 * ensure messages from different ranks are printed in a rank-ordered fashion
 * (i.e., to avoid interleaving). It also filters out messages if the current
 * function is not in the allow-list (`is_function_allowed(__func__)`) or the
 * requested log level is higher than `get_log_level()`.
 *
 * With GLOBAL scope every rank prints its own line, in rank order, so the message should
 * carry per-rank data; a value every rank shares belongs in `LOG_ALLOW(GLOBAL, ...)`, which
 * prints it once. The flush is collective over the scope's communicator and runs only when
 * the message is enabled. The level and the allow-list are the same on every rank, so either
 * every rank of the communicator flushes or none does: a disabled message costs no
 * communication. Every rank of that communicator must still reach the call.
 *
 * @param scope  Either LOCAL (MPI_COMM_SELF) or GLOBAL (MPI_COMM_WORLD).
 * @param level  One of LOG_ERROR, LOG_WARNING, LOG_INFO, LOG_DEBUG, LOG_TRACE, LOG_VERBOSE.
 * @param fmt    A `printf`-style format string (e.g., "Message: %ld\n").
 * @param ...    Variadic arguments to fill in `fmt`.
 *
 * Example usage:
 *
 * \code{.c}
 * LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "[Rank %d] owns %d cells\n", rank, ncells);
 * \endcode
 */
#define LOG_ALLOW_SYNC(scope, level, fmt, ...)                                     \
do {                                                                               \
    /* ------------------------------------------------------------------ */      \
    /* Validate scope and pick communicator *before* any early exits.     */      \
    /* ------------------------------------------------------------------ */      \
    MPI_Comm _comm;                                                                \
    if      ((scope) == LOCAL)  _comm = MPI_COMM_SELF;                             \
    else if ((scope) == GLOBAL) _comm = MPI_COMM_WORLD;                            \
    else {                                                                        \
        fprintf(stderr, "LOG_ALLOW_SYNC ERROR: invalid scope (%d) at %s:%d\n",     \
                (scope), __FILE__, __LINE__);                                      \
        MPI_Abort(MPI_COMM_WORLD, 1);                                              \
    }                                                                              \
                                                                                   \
    /* ------------------------------------------------------------------ */      \
    /* Decide whether *this* rank should actually print.                   */      \
    /* ------------------------------------------------------------------ */      \
    PetscBool _doPrint =                                                          \
        LOG_LEVEL_ADMITS(level);   \
                                                                                   \
    if (_doPrint) {                                                                \
        PetscSynchronizedPrintf(_comm, "[%s] " fmt, __func__, ##__VA_ARGS__);      \
        /* Collective; _doPrint has the same value on every rank. */               \
        PetscSynchronizedFlush(_comm, PETSC_STDOUT);                               \
    }                                                                              \
} while (0)

/**
 * @brief True when `pid` is one of the run's traced particles (`simCtx->tracePIDs`, sorted).
 * @param[in] simCtx The simulation context holding the traced IDs.
 * @param[in] pid    The particle ID to test.
 * @return PETSC_TRUE when the particle is traced.
 */
static inline PetscBool ParticleIsTraced(const SimCtx *simCtx, PetscInt64 pid)
{
    PetscInt lo = 0, hi = simCtx->nTracePIDs - 1;
    while (lo <= hi) {
        const PetscInt mid = lo + (hi - lo) / 2;
        if (simCtx->tracePIDs[mid] == pid) return PETSC_TRUE;
        if (simCtx->tracePIDs[mid] < pid) lo = mid + 1;
        else hi = mid - 1;
    }
    return PETSC_FALSE;
}

/** @brief True when the calling function's per-particle lines print for every particle (VERBOSE and allow-listed). */
#define PARTICLE_LOG_VERBOSE_ENABLED() \
    (PetscUnlikely((int)get_log_level() >= (int)LOG_VERBOSE) && is_function_allowed(__func__))

/**
 * @brief True when LOG_PARTICLE() would print for this particle in the calling function.
 * @details For a stage whose line is costly to assemble (several fields): test this, build the
 *          text, then emit it with LOG_PARTICLE().
 */
#define PARTICLE_LOG_ACTIVE(simCtx, pid) \
    ((PetscBool)((PetscUnlikely((simCtx)->nTracePIDs > 0) && ParticleIsTraced((simCtx), (pid))) || PARTICLE_LOG_VERBOSE_ENABLED()))

/**
 * @brief Reports one particle at one stage of a step: the consolidated per-particle line.
 *
 * Each stage that processes particles (position update, location, interpolation, model
 * update, removal) emits at most one such line per particle, carrying every value the stage
 * computed for it. The line prints, on the rank holding the particle, when either
 * - the particle is traced (`monitor.yml -> logging.trace_particles`), whatever the log level
 *   and allow-list, prefixed `[trace pid=<id> step=<n> <stage>]`; or
 * - the log level is VERBOSE and the calling function is allow-listed, prefixed with the
 *   function name.
 * With neither, it costs two loads and two compares.
 *
 * @param simCtx The simulation context (`const SimCtx *`).
 * @param pid    The particle ID (`PetscInt64`).
 * @param stage  A short stage label, such as "locate".
 * @param fmt    A `printf`-style format string ending in a newline.
 */
#define LOG_PARTICLE(simCtx, pid, stage, fmt, ...)                                                   \
    do {                                                                                             \
        const PetscBool _picurv_traced =                                                             \
            (PetscBool)(PetscUnlikely((simCtx)->nTracePIDs > 0) && ParticleIsTraced((simCtx), (pid))); \
        if (_picurv_traced || PARTICLE_LOG_VERBOSE_ENABLED()) {                                      \
            PetscPrintf(PETSC_COMM_SELF, "[%s pid=%lld step=%" PetscInt_FMT " %s] " fmt,             \
                        _picurv_traced ? "trace" : __func__, (long long)(pid), (simCtx)->step,       \
                        (stage), ##__VA_ARGS__);                                                     \
        }                                                                                            \
    } while (0)

/**
 * @brief Logs a message inside a loop, but only every `interval` iterations.
 *
 * @param scope     LOCAL or GLOBAL.
 * @param level     LOG_* level.
 * @param iterVar   The loop variable (e.g., i).
 * @param interval  Only log when (iterVar % interval == 0).
 * @param fmt       printf-style format string.
 * @param ...       Variadic arguments to include in the formatted message.
 *
 * Example:
 *    for (int i = 0; i < 100; i++) {
 *        LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, i, 10, "Value of i=%d\n", i);
 *    }
 */
#define LOG_LOOP_ALLOW(scope, level, iterVar, interval, fmt, ...)              \
    do {                                                                       \
        if (LOG_LEVEL_ADMITS(level)) { \
            if ((iterVar) % (interval) == 0) {                                 \
                MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
                PetscPrintf(comm, "[%s] [%s=%" PetscInt_FMT "] " fmt,                       \
                            __func__, #iterVar, (iterVar), ##__VA_ARGS__);               \
            }                                                                  \
        }                                                                      \
    } while (0)
/*
#define LOG_LOOP_ALLOW(scope,level, iterVar, interval, fmt, ...)              \
    do {                                                                       \
        if ((int)(level) <= (int)get_log_level() && is_function_allowed(__func__)) { \
            if ((iterVar) % (interval) == 0) {                                 \
                MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
                PetscPrintf(comm, "[%s] [Iter=%d] " fmt,                       \
                            __func__, (iterVar), ##__VA_ARGS__);               \
            }                                                                  \
        }                                                                      \
    } while (0)
*/

/**
 * @brief Logs a custom message if a variable equals a specific value.
 *
 * This is a variadic macro for logging a single event when a condition is met.
 * It is extremely useful for printing debug information at a specific iteration
 * of a loop or when a state variable reaches a certain value.
 *
 * @param scope   Either LOCAL or GLOBAL.
 * @param level   The logging level.
 * @param var     The variable to check (e.g., a loop counter 'k').
 * @param val     The value that triggers the log (e.g., 6). The log prints if var == val.
 * @param fmt     A printf-style format string.
 * @param ...     A printf-style format string and its corresponding arguments.
 */
#define LOG_LOOP_ALLOW_EXACT(scope, level, var, val, fmt, ...)                  \
    do {                                                                       \
        /* First, perform the cheap, standard gatekeeper checks. */            \
        if (LOG_LEVEL_ADMITS(level)) { \
            /* Only if those pass, check the user's specific condition. */      \
            if ((var) == (val)) {                                              \
                MPI_Comm comm = ((scope) == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
                /* Print the standard prefix, then the user's custom message. */ \
                PetscPrintf(comm, "[%s] [%s=%" PetscInt_FMT "] " fmt,                          \
                    __func__, #var, (var), ##__VA_ARGS__);                      \
            }                                                                  \
        }                                                                      \
    } while (0)

/**
 * @brief Logs a single element of an array, given an index.
 *
 * @param scope   Either LOCAL or GLOBAL.
 * @param level   LOG_ERROR, LOG_WARNING, LOG_INFO, or LOG_DEBUG.
 * @param arr     Pointer to the array to log from.
 * @param length  The length of the array (to prevent out-of-bounds).
 * @param idx     The index of the element to print.
 * @param fmt     The printf-style format specifier (e.g. "%g", "%f", etc.).
 *
 * This macro only logs if:
 *  1) The current function is in the allow-list (`is_function_allowed(__func__)`).
 *  2) The requested logging `level` <= the current global `get_log_level()`.
 *  3) The index `idx` is valid (0 <= idx < length).
 */
#define LOG_ARRAY_ELEMENT_ALLOW(scope,level, arr, length, idx, fmt)          \
    do {                                                                      \
        if (LOG_LEVEL_ADMITS(level)) { \
            if ((idx) >= 0 && (idx) < (length)) {                             \
                MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
                PetscPrintf(comm, "[%s] arr[%d] = " fmt "\n",                 \
                            __func__, (idx), (arr)[idx]);                     \
            }                                                                 \
        }                                                                     \
    } while (0)

/**
 * @brief Logs a consecutive subrange of an array.
 *
 * @param scope   Either LOCAL or GLOBAL.
 * @param level   LOG_ERROR, LOG_WARNING, LOG_INFO, or LOG_DEBUG.
 * @param arr     Pointer to the array to log from.
 * @param length  Total length of the array.
 * @param start   Starting index of the subrange.
 * @param end     Ending index of the subrange (inclusive).
 * @param fmt     The printf-style format specifier (e.g., "%g", "%f").
 *
 * This macro prints each element arr[i] for i in [start, end], bounded by [0, length-1].
 */
#define LOG_ARRAY_SUBRANGE_ALLOW(scope,level, arr, length, start, end, fmt)  \
    do {                                                                      \
        if (LOG_LEVEL_ADMITS(level)) { \
            MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
            PetscInt _start = (start) < 0 ? 0 : (start);                      \
            PetscInt _end   = (end) >= (length) ? (length) - 1 : (end);       \
            for (PetscInt i = _start; i <= _end; i++) {                       \
                PetscPrintf(comm, "[%s] arr[%d] = " fmt "\n", __func__, i, (arr)[i]); \
            }                                                                 \
        }                                                                     \
    } while (0)

// --------------------- Function Declarations ---------------------

/**
 * @brief Read-only view of the cached logging level; -1 until the first get_log_level() call.
 * @details Only logging.c writes the level. The view lets get_log_level() be inlined into every
 *          logging macro, so a disabled message in a per-particle loop costs a load and a
 *          compare rather than a function call.
 */
extern const LogLevel *const _picurv_log_level;

/**
 * @brief Internal hook: read `LOG_LEVEL` from the environment and cache it, once.
 * @details Supported levels are "ERROR", "WARNING", "INFO", "DEBUG", "TRACE", and "VERBOSE".
 *          Unset or unrecognized values default to "ERROR". Call get_log_level() instead.
 * @return LogLevel The cached logging level.
 */
LogLevel _ResolveLogLevel(void);

/**
 * @brief Retrieves the current logging level from the environment variable `LOG_LEVEL`.
 *
 * The first call reads `LOG_LEVEL` (see _ResolveLogLevel()); every later call returns the
 * cached value.
 *
 * @return LogLevel The current logging level.
 */
static inline LogLevel get_log_level(void)
{
    const LogLevel level = *_picurv_log_level;
    return PetscLikely((int)level != -1) ? level : _ResolveLogLevel();
}

/**
 * @brief Prints the current logging level to the console.
 *
 * This function retrieves the log level using `get_log_level()` and prints
 * the corresponding log level name. It helps verify the logging configuration
 * at runtime.
 * The log levels supported are:
 * - `LOG_ERROR`   (0) : Logs only critical errors.
 * - `LOG_WARNING` (1) : Logs warnings and errors.
 * - `LOG_INFO`    (2) : Logs general information, warnings, and errors.
 * - `LOG_DEBUG`   (3) : Logs debugging information, info, warnings, and errors.
 * - `LOG_TRACE`   (4) : Logs fine-grained trace information.
 * - `LOG_VERBOSE` (5) : Logs very detailed developer output.
 * If `LOG_LEVEL` is not set, it defaults to `LOG_ERROR`.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode print_log_level(void);

/**
 * @brief Sets the global list of function names that are allowed to log.
 *
 * You can replace the entire list of allowed function names at runtime.
 *
 * @param functionList Replacement array of permitted function names.
 * @param count Number of entries in `functionList`.
 */
void set_allowed_functions(const char** functionList, int count);

/**
 * @brief Checks if a given function is in the allow-list.
 *
 * This helper is used internally by the LOG_ALLOW macro.
 *
 * @param functionName Function name to query.
 * @return PETSC_TRUE when the name is enabled by the allow-list.
 */
PetscBool is_function_allowed(const char* functionName);

/**
 * @brief Prints the coordinates of a cell's vertices.
 *
 * This function iterates through the eight vertices of a given cell and prints their
 * coordinates. It is primarily used for debugging purposes to verify the correctness
 * of cell vertex assignments.
 *
 * @param[in]  cell     Pointer to a `Cell` structure representing the cell, containing its vertices.
 * @param[in]  rank     MPI rank for identification (useful in parallel environments).
 * @return PetscErrorCode Returns 0 to indicate successful execution. Non-zero on failure.
 *
 * @note
 * - Ensure that the `cell` pointer is not `NULL` before calling this function..
 */
PetscErrorCode LOG_CELL_VERTICES(const Cell *cell, PetscMPIInt rank);

/**
 * @brief Prints the signed distances to each face of the cell.
 *
 * This function iterates through the six signed distances from a point to each face of a given cell
 * and prints their values. It is primarily used for debugging purposes to verify the correctness
 * of distance calculations.
 *
 * @param[in]  d        An array of six `PetscReal` values representing the signed distances.
 *                      The indices correspond to:
 *                      - d[LEFT]: Left Face
 *                      - d[RIGHT]: Right Face
 *                      - d[BOTTOM]: Bottom Face
 *                      - d[TOP]: Top Face
 *                      - d[FRONT]: Front Face
 *                      - d[BACK]: Back Face
 *
 * @return PetscErrorCode Returns 0 to indicate successful execution. Non-zero on failure.
 *
 * @note
 * - Ensure that the `d` array is correctly populated with signed distances before calling this function.
 */
PetscErrorCode LOG_FACE_DISTANCES(PetscReal* d);

/**
 * @brief Prints particle fields in a table that automatically adjusts its column widths.
 *
 * Collective. A row is printed for each particle whose ID is a multiple of `printInterval`,
 * at most PICURV_SNAPSHOT_ROWS_PER_RANK per rank, so successive snapshots show the same
 * particles and the table's size does not grow with the swarm. Column widths are agreed
 * across ranks; one header precedes every rank's rows, and a closing line states how many
 * rows of how many particles were shown.
 *
 * @param[in] user           Pointer to the UserCtx structure.
 * @param[in] printInterval  Sampling interval on particle IDs (values below 1 mean every particle).
 *
 * @return PetscErrorCode Returns 0 on success.
 */
PetscErrorCode LOG_PARTICLE_FIELDS(UserCtx* user, PetscInt printInterval);

/**
 * @brief Returns whether periodic particle console snapshots are enabled.
 *
 * This checks only the reporting contract (particles exist, cadence is enabled,
 * and the global log level is at least INFO).
 *
 * @param simCtx Simulation context controlling the operation.
 * @return PetscBool indicating the result of `IsParticleConsoleSnapshotEnabled()`.
 */
PetscBool IsParticleConsoleSnapshotEnabled(const SimCtx *simCtx);

/**
 * @brief Returns whether a particle console snapshot should be emitted for the
 *
 * completed timestep.
 *
 * @param simCtx Simulation context controlling the operation.
 * @param completed_step Completed step index used by the decision helper.
 * @return PetscBool indicating the result of `ShouldEmitPeriodicParticleConsoleSnapshot()`.
 */
PetscBool ShouldEmitPeriodicParticleConsoleSnapshot(const SimCtx *simCtx, PetscInt completed_step);

/**
 * @brief Emits one particle console snapshot into the main solver log.
 *
 * @param user Primary `UserCtx` input for the operation.
 * @param simCtx Simulation context controlling the operation.
 * @param step Step index associated with the operation.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode EmitParticleConsoleSnapshot(UserCtx *user, SimCtx *simCtx, PetscInt step);

/**
 * @brief Reports whether the periodic statistics console snapshot is enabled.
 *
 * Mirrors the particle console gate: the subsystem must be active, the configured
 * cadence positive, and the effective log level at least `LOG_INFO`.
 *
 * @param[in] simCtx Simulation context to inspect.
 * @return `PETSC_TRUE` when snapshots should be emitted.
 */
PetscBool IsStatisticsConsoleSnapshotEnabled(const struct SimCtx *simCtx);

/**
 * @brief Reports whether a completed step falls on the console snapshot cadence.
 * @param[in] simCtx         Simulation context to inspect.
 * @param[in] completed_step Step that has just completed.
 * @return `PETSC_TRUE` when a snapshot is due for this step.
 */
PetscBool ShouldEmitPeriodicStatisticsConsoleSnapshot(const struct SimCtx *simCtx, PetscInt completed_step);

/**
 * @brief Emits one console snapshot of window progress.
 *
 * Reports window-level scalars only; it never dumps field data.
 *
 * @param[in] user   Finest-level block array supplying accumulator state.
 * @param[in] simCtx Simulation context carrying the window array.
 * @param[in] step   Step the snapshot describes.
 * @return Zero on success.
 */
PetscErrorCode EmitStatisticsConsoleSnapshot(UserCtx *user, const struct SimCtx *simCtx, PetscInt step);

/* ------------------------------------------------------------------------- */
/**
 * @brief Free an array previously returned by LoadAllowedFunctionsFromFile().
 *
 * @param[in,out] funcs Array of strings to release (may be @c NULL).
 * @param[in]     n     Number of entries in @p funcs.  Ignored if @p funcs is
 *                      @c NULL.
 *
 * @return 0 on success or a PETSc error code.
 */
PetscErrorCode FreeAllowedFunctions(char **funcs, PetscInt n);

/**
 * @brief Load function names from a text file.
 *
 * The file is expected to contain **one identifier per line**.  Blank lines and
 * lines whose first non‑blank character is a <tt>#</tt> are silently skipped so
 * the file can include comments.  Example:
 *
 * @code{.txt}
 * # Allowed function list
 * main
 * InitializeSimulation
 * InterpolateAllFieldsToSwarm  # inline comments are OK, too
 * @endcode
 *
 * The routine allocates memory as needed (growing an internal buffer with
 * @c PetscRealloc()) and returns the resulting array and its length to the
 * caller.  Use FreeAllowedFunctions() to clean up when done.
 *
 * @param[in]  filename  Path of the configuration file to read.
 * @param[out] funcsOut  On success, points to a freshly‑allocated array of
 *                       <tt>char*</tt> (size @p nOut).
 * @param[out] nOut      Number of valid entries in @p funcsOut.
 *
 * @return 0 on success, or a PETSc error code on failure (e.g. I/O error, OOM).
 */
PetscErrorCode LoadAllowedFunctionsFromFile(const char   filename[],
                                            char      ***funcsOut,
                                            PetscInt   *nOut);


/**
 * @brief Returns the canonical log token for a boundary-face enum value.
 * @param[in] face The BCFace enum value.
 * @return Pointer to a constant string representing the face.
 */
const char* BCFaceToString(BCFace face);

/**
 * @brief Convert an initial-condition mode to a string representation.
 * @param[in] mode Initial-condition mode value.
 * @return Pointer to a constant string representing the initial-condition mode.
 */
const char* InitialConditionModeToString(InitialConditionMode mode);

/**
 * @brief Convert a FlowDirection enum value to its YAML token string.
 * @param[in] fd FlowDirection value.
 * @return Token string such as "+Zeta", or "from INLET" when FLOW_DIR_UNSET.
 */
const char* FlowDirectionToString(FlowDirection fd);

/**
 * @brief Returns the canonical log token for a particle-initialization mode.
 * @param[in] ParticleInitialization The ParticleInitialization enum value.
 * @return Pointer to a constant string representing the particle initialization type.
 */
const char* ParticleInitializationToString(ParticleInitializationType ParticleInitialization);

/**
 * @brief Returns the canonical log token for an LES model selector.
 * @param[in] LESFlag The LES flag value.
 * @return Pointer to a constant string representing the LES Flag.
 */
const char* LESModelToString(LESModelType LESFlag);

/**
 * @brief Returns the user-facing name of an LES grid-filter-width model.
 * @param[in] model The filter-width model selector.
 * @return Pointer to a constant string naming the model.
 */
const char* LESFilterWidthModelToString(LESFilterWidthModel model);

/**
 * @brief Returns the user-facing name of an LES test-filter kernel.
 * @param[in] kernel The test-filter kernel selector.
 * @return Pointer to a constant string naming the kernel.
 */
const char* LESTestFilterKernelToString(LESTestFilterKernel kernel);

/**
 * @brief Returns the user-facing name of an LES coefficient-averaging mode.
 * @param[in] mode The averaging mode selector.
 * @return Pointer to a constant string naming the mode.
 */
const char* LESAveragingModeToString(LESAveragingMode mode);

/**
 * @brief Returns the user-facing name of an LES coefficient-limiting mode.
 * @param[in] mode The limiting mode selector.
 * @return Pointer to a constant string naming the mode.
 */
const char* LESClipModeToString(LESClipMode mode);

/**
 * @brief Returns the user-facing name of a wall-function model.
 * @param[in] model The wall-model selector.
 * @return Pointer to a constant string naming the model, or "none" when disabled.
 */
const char* WallFunctionModelToString(WallFunctionModel model);

/**
 * @brief Returns the canonical log token for a momentum-solver selector.
 * @param[in] SolverFlag The Momentum Solver flag value.
 * @return Pointer to a constant string representing the MomentumSolverType.
 */
const char* MomentumSolverTypeToString(MomentumSolverType SolverFlag);

/**
 * @brief Returns the canonical log token for a boundary mathematical type.
 * @param[in] type The BCType enum value.
 * @return Pointer to a constant string representing the BC type.
 */
const char* BCTypeToString(BCType type);

/**
 * @brief Returns the configuration token for a particle boundary policy.
 * @param[in] policy The ParticleBoundaryPolicy enum value.
 * @return Pointer to a constant string: "remove", "reflect" or "wrap".
 */
const char* ParticleBoundaryPolicyToString(ParticleBoundaryPolicy policy);

/**
 * @brief Returns the display name of a grid-to-particle interpolation method.
 * @param[in] method The InterpolationMethod enum value.
 * @return Pointer to a constant string naming the method, as the banner and logs print it.
 */
const char* InterpolationMethodToString(InterpolationMethod method);

/**
 * @brief Converts a BCHandlerType enum to its string representation.
 *
 * Provides a descriptive string for a specific boundary condition implementation strategy.
 * This is crucial for logging the exact behavior configured for a face.
 *
 * @param handler_type The BCHandlerType enum value (e.g., BC_HANDLER_WALL_NOSLIP).
 * @return A constant character string corresponding to the enum. Returns
 *         "UNKNOWN_HANDLER" if the enum value is not recognized.
 */
const char* BCHandlerTypeToString(BCHandlerType handler_type);

/**
 * @brief A custom KSP monitor that logs to a file and optionally to the console.
 *
 * This function unconditionally calls the standard true residual monitor to log to a
 * file viewer provided in the context. It also checks a flag in the context
* and, if true, calls the monitor again to log to standard output.
 *
 * @param ksp   The Krylov subspace context.
 * @param it    The current iteration number.
 * @param rnorm The preconditioned residual norm.
 * @param ctx   A pointer to the DualMonitorCtx structure.
 * @return      PetscErrorCode 0 on success.
 */
PetscErrorCode DualKSPMonitor(KSP ksp, PetscInt it, PetscReal rnorm, void *ctx);

/**
 * @brief Destroys the DualMonitorCtx.
 *
 * This function is passed to KSPMonitorSet to ensure the viewer is
 * properly destroyed and the context memory is freed when the KSP is destroyed.
 * @param ctx a pointer to the context pointer to be destroyed
 * @return PetscErrorCode
 */
PetscErrorCode DualMonitorDestroy(void **ctx);

/**
 * @brief Logs continuity metrics for a single block to a file.
 *
 * This function should be called for each block, once per timestep. It opens a
 * central log file in append mode. To ensure the header is written only once,
 * it checks if it is processing block 0 on the simulation's start step.
 *
 * @param user A pointer to the UserCtx for the specific block whose metrics
 *             are to be logged. The function accesses both global (SimCtx)
 *             and local (user->...) data.
 * @return     PetscErrorCode 0 on success.
 */
PetscErrorCode LOG_CONTINUITY_METRICS(UserCtx *user);

/**
 * @brief Logs physical solution-convergence metrics once per completed timestep.
 *
 * This logger is intentionally separate from the detailed inner solver-health
 * logs. Depending on the configured mode, it records deterministic step drift,
 * periodic phase-aligned drift, or statistical window drift into
 * `logs/solution_convergence.log`.
 *
 * The logger is solver-mode only and is intended to run after a completed
 * flow-solve / projection step but before history vectors are shifted forward.
 * Warmup rows are handled internally:
 * - steady/transient mode: the first logged step has `has_reference = 0`
 * - periodic mode: the first cycle fills the phase reference ring
 * - statistical mode: window-drift fields remain zero until enough samples
 *   exist to form the requested windows
 *
 * Existing detailed inner logs remain the authoritative source for momentum,
 * Poisson, and continuity health. This logger only summarizes physical
 * solution-drift metrics.
 *
 * @param[in,out] simCtx Master simulation context controlling the logging mode
 *                       and owning any supporting runtime buffers.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode LOG_SOLUTION_CONVERGENCE(SimCtx *simCtx);

/**
 * @brief A function that outputs the name of the current level in the ParticleLocation enum.
 * @param level The ParticleLocation enum value.
 * @return A constant character string corresponding to the enum. Returns
 *        "UNKNOWN_LEVEL" if the enum value is not recognized.
 */
const char* ParticleLocationStatusToString(ParticleLocationStatus level); 

/*================================================================================*
 *                          PROGRESS BAR UTILITY                                  *
 *================================================================================*/

/**
 * @brief Prints a progress bar to the console.
 *
 * This function should only be called by the root process (rank 0). It uses
 * a carriage return `\r` to overwrite the same line in the terminal, creating
 * a dynamic progress bar.
 *
 * @param step             The current step index from the loop (e.g., from 0 to N-1).
 * @param startStep        The global starting step number of the simulation.
 * @param totalSteps       The total number of steps to be run in this simulation instance.
 * @param currentTime      The current simulation time to display.
 * @param seconds_per_step Mean wall seconds per completed step, with the time left at that rate;
 *                         omitted when not positive.
 * @param lost_cumulative  Particles lost so far in the run; omitted when negative.
 */
void PrintProgressBar(PetscInt step, PetscInt startStep, PetscInt totalSteps, PetscReal currentTime,
                      PetscReal seconds_per_step, PetscInt64 lost_cumulative);

/**
 * @brief Applies the profiling configuration from SimCtx to PICurv's PETSc log events.
 *
 * Records the per-step `selected` function list and the enabled helper timers
 * (`profiling.helpers`), and switches each helper event registered so far on or off to
 * match. Starts PETSc's default log handler if nothing has yet. Called once from
 * CreateSimulationContext(), before the main time loop.
 *
 * @param simCtx The master simulation context holding the profiling configuration; its
 *               name lists must outlive profiling.
 * @return PetscErrorCode
 */
PetscErrorCode ProfilingInitialize(SimCtx *simCtx);

/**
 * @brief Starts a new per-step profiling interval.
 *
 * Records every PETSc event's running totals; the next ProfilingLogTimestepSummary()
 * reports the difference. Cumulative totals are untouched.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ProfilingResetTimestepCounters(void);

/**
 * @brief Logs the performance summary for the current timestep and resets timers.
 *
 * Depending on the configured profiling timestep mode, this function writes
 * per-step profiling rows to the configured profiling log file:
 * - `off`: writes nothing
 * - `selected`: writes only functions marked for per-step reporting
 * - `all`: writes every instrumented function that ran in the timestep
 *
 * It must be called once per timestep, typically at the end of the main loop. A step's
 * values are the change in each PICurv event's PETSc running totals since the previous
 * call (or ProfilingResetTimestepCounters()); an event still open, such as the enclosing
 * time loop, adds its time only when it ends. Only rank 0 writes.
 *
 * @param simCtx The simulation context holding profiling output settings.
 * @param step The current simulation step number, for logging context.
 * @return PetscErrorCode
 */
PetscErrorCode ProfilingLogTimestepSummary(SimCtx *simCtx, PetscInt step);

/**
 * @brief Append a reduced runtime memory sample to the configured memory log.
 *
 * The log is intentionally compact: every rank samples process/PETSc allocator
 * memory, one max reduction summarizes the job, and rank 0 appends one
 * terminal-readable row.
 *
 * @param simCtx The simulation context holding log directory and memory log settings.
 * @param step The solver/post step associated with the sample.
 * @param event Human-readable event label such as "Step", "Post", "Shutdown", or "Final".
 * @param reason Human-readable reason, or NULL/"-" when not applicable.
 * @return PetscErrorCode
 */
PetscErrorCode RuntimeMemoryLogSample(SimCtx *simCtx, PetscInt step, const char *event, const char *reason);


/**
 * @brief Writes the end-of-run profiling summary across all ranks; collective on PETSC_COMM_WORLD.
 *
 * Two sections: PICurv's events (aggregates, and helpers that were enabled), then PETSc's own
 * events, whose inclusive times nest inside PICurv's. In each, the names any rank recorded are
 * united by name, each name's total time is reduced to its maximum and minimum (with the ranks
 * holding them) and its mean over ranks, and its calls are summed and reduced to their minimum
 * and maximum per rank (with ranks), so a time spread can be told apart from a call-count spread.
 * Rank 0 writes `ProfilingSummary_<mode>.log`, slowest rank first. This is the profiler's only
 * communication; the timed path stays local to each rank.
 *
 * @param simCtx  The Simulation Context Structure that can contains all the data regarding the simulation.
 * 
 * @return        PetscErrorCode 0 on success.
 */
PetscErrorCode ProfilingFinalize(SimCtx *simCtx);

// --- Internal functions, do not call directly ---
// Called by the macros below.
/**
 * @brief Internal profiling hook: register a call site's function as a PETSc log event.
 * @details Called once per call site, on its first execution. Registers the PICurv event
 *          classes and starts PETSc's default log handler on the first call anywhere. A helper
 *          event is switched off as it registers unless its name is enabled
 *          (`profiling.helpers`). Before PetscInitialize, after PetscFinalize, or on failure, it
 *          registers nothing and the call site stays untimed, trying again on its next call.
 * @param func_name Function name (`__func__`), the event's name.
 * @param helper    PETSC_TRUE for a per-particle or per-cell helper, PETSC_FALSE for an aggregate.
 * @param event     Call site's cached event; set when registration succeeds, else left at -1.
 */
void _ProfilingRegister(const char *func_name, PetscBool helper, PetscLogEvent *event);

/**
 * @brief Read-only view: whether `profiling.helpers` enables any helper timer.
 * @details False until ProfilingInitialize() applies a non-empty list. While it is false every
 *          helper event is switched off, so PROFILE_HELPER_BEGIN skips PETSc's per-event check.
 */
extern const PetscBool *const _picurv_profile_helpers_on;

/**
 * @brief Marks the beginning of a profiled function: an always-on aggregate timer.
 *
 * Place this macro at the very beginning of a function called once or a few times per step
 * or solver pass. The function is a PETSc log event of class `PICurv`, named by `__func__`,
 * registered on its first call and cached in a function-local static; it appears in the
 * final profiling summary, the per-step output and `log_view`.
 */
#define PROFILE_FUNCTION_BEGIN \
    static PetscLogEvent _picurv_profile_event = -1; \
    if (PetscUnlikely(_picurv_profile_event < 0)) _ProfilingRegister(__func__, PETSC_FALSE, &_picurv_profile_event); \
    const PetscBool _picurv_profile_timed = (PetscBool)(_picurv_profile_event >= 0); \
    if (_picurv_profile_timed) (void)PetscLogEventBegin(_picurv_profile_event, 0, 0, 0, 0)

/**
 * @brief Marks the beginning of a profiled helper: timed only when enabled by name.
 *
 * For functions called per particle or per cell. The event (class `PICurvHelper`) is off
 * unless `profiling.helpers` names it (or is `all`). With no helper enabled, a call costs one
 * flag test; with some enabled, PETSc's per-event check decides.
 * End it with PROFILE_FUNCTION_END (or its alias PROFILE_HELPER_END).
 */
#define PROFILE_HELPER_BEGIN \
    static PetscLogEvent _picurv_profile_event = -1; \
    if (PetscUnlikely(_picurv_profile_event < 0)) _ProfilingRegister(__func__, PETSC_TRUE, &_picurv_profile_event); \
    const PetscBool _picurv_profile_timed = (PetscBool)(*_picurv_profile_helpers_on && _picurv_profile_event >= 0); \
    if (_picurv_profile_timed) (void)PetscLogEventBegin(_picurv_profile_event, 0, 0, 0, 0)

/**
 * @brief Marks the end of a profiled code block.
 *
 * Place this macro just before every return point in a function that starts with
 * PROFILE_FUNCTION_BEGIN or PROFILE_HELPER_BEGIN. A missed one leaves the PETSc event open,
 * and the function's later calls are then not recorded. It ends the event only if the BEGIN
 * macro began it, and reads the state that macro declares, so it must sit in that macro's
 * scope; a function that ends a timer it never began does not compile.
 */
#define PROFILE_FUNCTION_END \
    do { if (_picurv_profile_timed) (void)PetscLogEventEnd(_picurv_profile_event, 0, 0, 0, 0); } while (0)

/** @brief Alias of PROFILE_FUNCTION_END for a function that began with PROFILE_HELPER_BEGIN. */
#define PROFILE_HELPER_END PROFILE_FUNCTION_END


/**
 * @brief Computes and logs the local and global min/max values of a 3-component vector field.
 *
 * This utility function inspects a PETSc Vec associated with a DMDA and calculates the
 * minimum and maximum values for each of its three components (e.g., x, y, z) both for the
 * local data on the current MPI rank and for the entire global domain.
 *
 * It uses the same "smart" logic as the flow solver, ignoring the padding nodes at the
 * IM, JM, and KM boundaries of the grid. The results are printed to the standard output
 * in a formatted, easy-to-read table.
 *
 * @param[in] user      Pointer to the user-defined context. Used for grid information (IM, JM, KM)
 *                      and MPI rank.
 * @param[in] field_id Typed identity of the field being analyzed.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode LOG_FIELD_MIN_MAX(UserCtx *user, FieldId field_id);

/**
 * @brief Logs the anatomy of a specified field at key boundary locations,
 *        respecting the solver's specific grid and variable architecture.
 *
 * This intelligent diagnostic function inspects a PETSc Vec and prints its values
 * at critical boundary locations (-Xi/+Xi, -Eta/+Eta, -Zeta/+Zeta). It is "architecture-aware":
 *
 * - **Cell-Centered Fields ("Ucat", "P"):** It correctly applies the "Shifted Index Architecture,"
 *   where the value for geometric `Cell i` is stored at array index `i+1`. It labels
 *   the output to clearly distinguish between true physical values and ghost values.
 * - **Single-Face-Family Fields:** `Csi/ICsi/IEta/IZet/Centx` belong to the
 *   I-face family, with corresponding J- and K-face families.
 * - **Component-Staggered Fields ("Ucont" and histories):** x/y/z components
 *   live on I/J/K faces respectively.
 * - **Node-Centered Fields ("Coordinates"):** It uses a direct index mapping, where the value for
 *   `Node i` is stored at index `i`.
 *
 * The output is synchronized across MPI ranks to ensure readability and focuses on a
 * slice through the center of the domain to be concise.
 *
 * @param user       A pointer to the UserCtx structure containing the DMs and Vecs.
 * @param field_id   Typed identity of the persistent Eulerian field to log.
 * @param stage_name A string identifier for the current simulation stage (e.g., "After Advection").
 * @return           PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode LOG_FIELD_ANATOMY(UserCtx *user, FieldId field_id, const char *stage_name);

/**
 * @brief Logs the node-layout anatomy of the transient center-to-corner interpolation field.
 * @param user       Context owning the corner-staging workspace.
 * @param corner_field_id Which corner workspace the caller populated.
 * @param stage_name Printable simulation stage.
 * @return Zero on success.
 */
PetscErrorCode LOG_CORNER_FIELD_ANATOMY(UserCtx *user, FieldId corner_field_id, const char *stage_name);

/**
 * @brief Logs the interpolation error between the analytical and computed solutions.
 *
 * @param user Primary `UserCtx` input for the operation.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode LOG_INTERPOLATION_ERROR(UserCtx *user);

/**
 * @brief Logs particle-to-grid scatter verification metrics for the prescribed scalar truth path.
 *
 * @param user Primary `UserCtx` input for the operation.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode LOG_SCATTER_METRICS(UserCtx *user);

/**
 * @brief Resets the aggregate per-timestep search instrumentation counters.
 *
 * @param simCtx Simulation context whose search metrics should be zeroed.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ResetSearchMetrics(SimCtx *simCtx);

/**
 * @brief Computes advanced particle statistics and stores them in SimCtx.
 *
 * This function calculates:
 * - Particle load imbalance across MPI ranks.
 * - The total number of grid cells occupied by at least one particle.
 *
 * It requires that CalculateParticleCountPerCell() has been called prior to its
 * execution. It uses collective MPI operations and must be called by all ranks.
 *
 * @param user Pointer to the UserCtx.
 * @return     PetscErrorCode 0 on success.
 */
PetscErrorCode CalculateAdvancedParticleMetrics(UserCtx *user);

/**
 * @brief Writes compact runtime search metrics to CSV and optionally to console.
 *
 * The CSV artifact is always written for particle-enabled runs. Console output
 * remains gated by normal logging level and function allow-listing.
 *
 * @param user Pointer to the UserCtx.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode LOG_SEARCH_METRICS(UserCtx *user);

/**
 * @brief Logs particle swarm metrics, adapting its behavior based on a boolean flag in SimCtx.
 *
 * This function serves a dual purpose:
 * 1. If simCtx->isInitializationPhase is PETSC_TRUE, it logs settlement
 *    diagnostics to "Initialization_Metrics.log", using the provided stageName.
 * 2. If simCtx->isInitializationPhase is PETSC_FALSE, it logs regular
 *    timestep metrics to "Particle_Metrics.log".
 *
 * @param user      A pointer to the UserCtx.
 * @param stageName A descriptive label recorded in the metrics log (for example,
 *                  initialization stage name or "Timestep Metrics").
 * @return          PetscErrorCode 0 on success.
 */
PetscErrorCode LOG_PARTICLE_METRICS(UserCtx *user, const char *stageName);

/**
 * @brief Opens a per-run diagnostics CSV in the run's analysis directory for appending.
 *
 * The runtime diagnostics files share a lifecycle rather than a schema: they live in
 * `simCtx->analysis_dir`, are appended to across a run, carry their column header only
 * once, and mark the seam where a continuation resumed so a reader can tell a restart
 * from a discontinuity in the data. This owns that lifecycle; the caller owns the
 * columns.
 *
 * Call on rank 0 only - the file is a serial artifact, and a collective write to it
 * would interleave. The caller closes the handle.
 *
 * @param[in]  simCtx   Simulation context supplying the analysis directory, step, and
 *                      continuation state.
 * @param[in]  filename Basename of the CSV, including the extension.
 * @param[in]  header   Column header written when the file is created, without a
 *                      trailing newline.
 * @param[out] file     Open append handle, never NULL on success.
 * @return PetscErrorCode 0 on success; `PETSC_ERR_FILE_OPEN` if the file cannot be opened.
 */
PetscErrorCode PicurvOpenDiagnosticsCsv(const SimCtx *simCtx, const char *filename,
                                        const char *header, FILE **file);
#endif // LOGGING_H
