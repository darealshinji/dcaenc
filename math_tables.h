/*
 * This file is part of dcaenc.
 *
 * Copyright (c) 2008-2012 Alexander E. Patrakov <patrakov@gmail.com>
 *
 * dcaenc is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * dcaenc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with dcaenc; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

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
