#ifndef MATH_TABLES_H
#define MATH_TABLES_H

#include <stdint.h>

#define AUBANDS 25

extern const int32_t lfe_fir[512];
extern const int32_t band_interpolation[2][512];
extern const int32_t band_spectrum[2][8];

extern const int32_t cos_table[2048];
extern const int bitrev[256];
extern const int32_t auf[9][AUBANDS][256];
extern const int32_t cb_to_level[2048];
extern const int32_t cb_to_add[256];

extern const int quant_levels_cb[27];

static inline int32_t cos_t(int x)
{
	return cos_table[x & 2047];
}

static inline int32_t sin_t(int x)
{
	return cos_t(x - 512);
}

#endif
