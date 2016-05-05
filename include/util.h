#ifndef UTIL_H
#define UTIL_H
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

void die(int exit_code, const char *msg, ...);
void *safe_malloc(size_t size);
void *safe_calloc(size_t nitems, size_t size);
void *safe_realloc(void *ptr, size_t size);

#endif /* defined UTIL_H */
