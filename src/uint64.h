// Declare 64-bit unsigned integer type
// uint64 is the type while UINT64MAX is the largest value that uint64 
// can hold.

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif
#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif

#if HAVE_UINT64_T

typedef uint64_t uint64;
#define UINT64MAX UINT64_MAX

#elif HAVE_UINT_LEAST64_T

typedef uint_least64_t uint64;
#define UINT64MAX UINT_LEAST64_MAX

#elif HAVE_UNSIGNED_LONG_LONG

typedef unsigned long long uint64;
#ifdef ULLONG_MAX
#define UINT64MAX ULLONG_MAX
#else
#define UINT64MAX 0xffffffffffffffffull
#endif // #ifdef ULLONG_MAX

#else

typedef unsigned long uint64;
#ifdef ULONG_MAX
#define UINT64MAX ULONG_MAX
#else
#if SIZEOF_UNSIGNED_LONG >= 8
#define UINT64MAX 0xfffffffffffffffful
#else
#define UINT64MAX 0xfffffffful
#endif // #if SIZEOF_ULONG >= 8
#endif // #ifdef ULONG_MAX

#endif
