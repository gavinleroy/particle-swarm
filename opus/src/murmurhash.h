/**
 * `murmurhash.h' - murmurhash
 *
 * copyright (c) 2014-2022 joseph werle <joseph.werle@gmail.com>
 * 
 * MIT License
 * https://github.com/jwerle/murmurhash.c
 */

#ifndef MURMURHASH_H
#define MURMURHASH_H

#include <stdint.h>

#define MURMURHASH_VERSION "0.1.0"



/**
 * Returns a murmur hash of `key' based on `seed'
 * using the MurmurHash3 algorithm
 */

uint32_t
murmurhash (const char *, uint32_t, uint32_t);

#endif
