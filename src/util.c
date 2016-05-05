#include <signal.h>
#include "util.h"

/*
 * Do not clean at exit when debugging.
 */
#ifdef NDEBUG
#define EXIT(CODE)\
    do {exit((CODE));} while (0)
#else
#define EXIT(CODE)\
    do {(void)(CODE); abort();} while (0)
#endif

void
die(int exit_code, const char *msg, ...)
{
    va_list va;
    va_start(va, msg);
    vfprintf(stderr, msg, va);
    va_end(va);
    fputc('\n', stderr);
    EXIT(exit_code);
}

void *
safe_malloc(size_t size)
{
    void *result = malloc(size);
    if (result == NULL){
        perror("malloc failed.");
        die(1, "Unable to allocate memory");
    }
    return result;
}

void *
safe_calloc(size_t nitems, size_t size)
{
    void *result = calloc(nitems, size);
    if (result == NULL){
        perror("calloc failed.");
        die(1, "Unable to allocate memory");
    }
    return result;
}

void *
safe_realloc(void *ptr, size_t size)
{
    void *result = realloc(ptr, size);
    if (result == NULL){
        perror("realloc failed.");
        die(1, "Unable to reallocate memory");
    }
    return result;
}
