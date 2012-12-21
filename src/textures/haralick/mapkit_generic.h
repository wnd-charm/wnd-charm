/* MapKit, Generic functions */

/*
 * Copyright (c) 2002-2004, Jean-Sebastien Roy (js@jeannot.org)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/* @(#) $Jeannot: mapkit_generic.h,v 1.31 2004/04/10 14:17:03 js Exp $ */

/*
Important defines for MapKit:
  MAPKIT_EXITONERROR : if defined, will exit upon failure.
    The return code is the MAPKIT error number.
  MAPKIT_DEBUG : if defined, will print information about errors
    and various events. Messages get printed on stderr.
  MAPKIT_COLLISIONS : if defined, statistics about hash collisions will be
    gathered.
*/

#ifndef _MAPKIT_GENERIC_
#define _MAPKIT_GENERIC_

#define MAPKIT_VERSION "1.4"

#include <stddef.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>

/* MS Visual C does not recognize 'inline' */
#if defined(_MSC_VER) && ! defined(__cplusplus) && ! defined(inline)
#define inline __inline
#endif

#if defined(__cplusplus)
extern "C" {
#endif

typedef enum
{
  MAPKIT_OK = 0,
  MAPKIT_ETOOBIG = 1,
  MAPKIT_ENOMEM = 2,
  MAPKIT_EBADKEY = 3,
  MAPKIT_EKEYNOTFOUND = 4
}
mapkit_error;

/* Constants */
extern const char *mapkit_error_msg[5];

/* Macros */

/* print an appropriate error message if _err != MAPKIT_OK */
#define CHKMAPKIT(_err) { mapkit_error __err = (_err); if (__err) fprintf(stderr, "MAPKIT: ERROR: %s (%s, line %d)\n", mapkit_error_msg[__err], __FILE__, __LINE__); }

/* if MAPKIT_DEBUG is set, most mapkit functions will output various
   informational & error messages */
#ifdef MAPKIT_DEBUG
#define MAPKIT_MSG(_err) fprintf(stderr, "MAPKIT: ERROR: %s (%s, line %d)\n", mapkit_error_msg[_err], __FILE__, __LINE__)
#else
#define MAPKIT_MSG(_err)
#endif

#define MAPKIT_FATAL_ERROR(_err) { MAPKIT_MSG(_err); exit(_err); }

/* if MAPKIT_EXITONERROR is set, the program will quit instead of returning
   an error message */
#ifdef MAPKIT_EXITONERROR
#define MAPKIT_ERROR(_err) MAPKIT_FATAL_ERROR(_err)
#define MAPKIT_ERROR_NORET(_err) MAPKIT_FATAL_ERROR(_err)
#else
#define MAPKIT_ERROR(_err) { MAPKIT_MSG(_err); return (_err); }
#define MAPKIT_ERROR_NORET(_err) MAPKIT_MSG(_err)
#endif

/* return index value if key not found */
#define MAPKIT_KEYNOTFOUND (-1)

/* special key values */
#define MAPKIT_FULLSLOT (0)
#define MAPKIT_FREESLOT (-1)
#define MAPKIT_DELETEDSLOT (-2)

/* default initial size */
#define MAPKIT_DEFAULT_EXPECTEDUSED (100)

/* types */
/* sizes, signed to accomodate MAPKIT_KEYNOTFOUND */
#define mapkit_size_t long
/* hash values, must be unsigned (for index, decrement) */
#define mapkit_hash_t unsigned long

/* find the lowest number in mapkit_primes greater than n */
extern mapkit_size_t mapkit_nextprime(const mapkit_size_t n);

/* string hash function */
extern mapkit_hash_t mapkit_strhash(const char *string);

/* memory hash function */
extern mapkit_hash_t mapkit_memhash(const void *data, const size_t len);

/* integer hash function */
static mapkit_hash_t mapkit_hash(mapkit_hash_t key);

/* Thomas Wang's integer hash function */
#if LONG_MAX <= 0x7fffffff
/* 32 bits version */
mapkit_hash_t mapkit_hash(mapkit_hash_t key)
{
  key += ~(key << 15);
  key ^=  (key >> 10);
  key +=  (key << 3);
  key ^=  (key >> 6);
  key += ~(key << 11);
  key ^=  (key >> 16);
  return key;
}
#else
/* 64 bits version */
mapkit_hash_t mapkit_hash(mapkit_hash_t key)
{
  key += ~(key << 32);
  key ^=  (key >> 22);
  key += ~(key << 13);
  key ^=  (key >> 8);
  key +=  (key << 3);
  key ^=  (key >> 15);
  key += ~(key << 27);
  key ^=  (key >> 31);
  return key;
}
#endif

#if defined(__cplusplus)
}
#endif

#endif /* _MAPKIT_GENERIC_ */
