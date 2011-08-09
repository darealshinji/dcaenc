/* 
 * This file is part of dcaenc.
 *
 * Copyright (c) 2008 Alexander E. Patrakov <patrakov@gmail.com>
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

#ifndef WAVFILE_H
#define WAVFILE_H

#include <stdio.h>

typedef struct {
	FILE * file;
	unsigned int channels;
	unsigned int bits_per_sample;
	unsigned int sample_rate;
	unsigned int samples_left;
} wavfile;

wavfile * wavfile_open(const char * filename);
int wavfile_read_s32(wavfile * f, int32_t * samples);
void wavfile_close(wavfile * f);

#endif
