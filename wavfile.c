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

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "wavfile.h"

static uint32_t find_chunk(FILE * file, const uint8_t chunk_id[4])
{
    uint8_t buffer[8];
    while (1) {
	size_t chunksize;
        size_t s = fread(buffer, 1, 8, file);
	if (s < 8)
	    return 0;
	chunksize = (uint32_t)buffer[4] | ((uint32_t)buffer[5] << 8) |
	    ((uint32_t)buffer[6] << 16) | ((uint32_t)buffer[7] << 24);
	if (!memcmp(buffer, chunk_id, 4))
	    return chunksize;
	fseek(file, chunksize, SEEK_CUR);
    }
}

wavfile * wavfile_open(const char * filename)
{
    wavfile *result;
    size_t s;
    uint8_t buffer[8];
    uint8_t *fmt;
    uint32_t v;
    uint32_t avg_bps;
    uint32_t block_align;
    static const uint8_t riff[4] = {'R', 'I', 'F', 'F'};
    static const uint8_t wave[4] = { 'W', 'A', 'V', 'E'};
    static const uint8_t fmt_[4] = {'f', 'm', 't', ' '};
    static const uint8_t data[4] = {'d', 'a', 't', 'a'};
    
    result = (wavfile *)calloc(1, sizeof(wavfile));
    if (!result)
	goto err0;
    
    result->file = fopen(filename, "rb");
    if (!result->file)
	goto err1;

    s = fread(buffer, 1, 8, result->file);
    if (s < 8)
	goto err2;

    if (memcmp(buffer, riff, 4))
	goto err2;

    /* TODO: check size (in buffer[4..8]) */
    s = fread(buffer, 1, 4, result->file);
    if (s < 4)
	goto err2;

    if (memcmp(buffer, wave, 4))
	goto err2;

    s = find_chunk(result->file, fmt_);
    if (s != 16 && s != 40)
	goto err2;

    fmt = (uint8_t*)malloc(s);
    if (!fmt)
	goto err2;

    if (fread(fmt, 1, s, result->file) != s)
	goto err3;

    /* wFormatTag */
    v = (uint32_t)fmt[0] | ((uint32_t)fmt[1] << 8);
    if (v != 1 && v != 0xfffe)
	goto err3;

    /* wChannels */
    v = (uint32_t)fmt[2] | ((uint32_t)fmt[3] << 8);
    if (v != 1 && v != 2 && v != 4 && v != 5 && v !=6)
	goto err3;
    result->channels = v;
    /* dwSamplesPerSec */
    result->sample_rate = (uint32_t)fmt[4] | ((uint32_t)fmt[5] << 8) |
	((uint32_t)fmt[6] << 16) | ((uint32_t)fmt[7] << 24);

    /* dwAvgBytesPerSec */
    avg_bps = (uint32_t)fmt[8] | ((uint32_t)fmt[9] << 8) |
	((uint32_t)fmt[10] << 16) | ((uint32_t)fmt[11] << 24);

    /* wBlockAlign */
    block_align = (uint32_t)fmt[12] | ((uint32_t)fmt[13] << 8);

    /* wBitsPerSample */
    result->bits_per_sample = (uint32_t)fmt[14] | ((uint32_t)fmt[15] << 8);
    if (result->bits_per_sample != 16 && result->bits_per_sample != 32)
	goto err3;

    if (block_align != result->channels * (result->bits_per_sample / 8))
	goto err3;

    if (avg_bps != block_align * result->sample_rate)
	goto err3;

    v = find_chunk(result->file, data);
    if (v == 0 || v % block_align != 0)
	goto err3;

    result->samples_left = v / block_align;
    free(fmt);
    return result;

    err3:
	free(fmt);
    err2:
	fclose(result->file);
    err1:
	free(result);
    err0:
	return NULL;
}

void wavfile_close(wavfile * f)
{
    fclose(f->file);
    free(f);
}

static int32_t get_s32_sample(const wavfile * f, const uint8_t *buffer, int sample, int channel)
{
    int offset = (f->bits_per_sample / 8) * (f->channels * sample + channel);
    uint32_t v;
    switch (f->bits_per_sample) {
	case 16:
	    v = (uint32_t)buffer[offset + 0] | ((uint32_t)buffer[offset + 1] << 8);
	    return v << 16;
	    break;
	case 32:
	    v = (uint32_t)buffer[offset + 0] | ((uint32_t)buffer[offset + 1] << 8) |
		((uint32_t)buffer[offset + 2] << 16) | ((uint32_t)buffer[offset + 3] << 24);
	    return v;
	    break;
	default:
	    return 0;
    }
}

int wavfile_read_s32(wavfile * f, int32_t *samples)
{
    uint8_t buffer[512 * 6 * 4];
    int32_t smpte_sample[6];
    int samples_to_read;
    int bytes_to_read;
    unsigned int i, ch;
    
    memset(buffer, 0, 512 * 6 * 4);
    samples_to_read = f->samples_left < 512 ? f->samples_left : 512;
    bytes_to_read = samples_to_read * f->channels * (f->bits_per_sample / 8);
    f->samples_left -= samples_to_read;
    if (fread(buffer, 1, bytes_to_read, f->file) != bytes_to_read) {
	f->samples_left = 0;
    }
    
    for (i = 0; i < 512; i++) {
	for (ch = 0; ch < f->channels; ch++)
	    smpte_sample[ch] = get_s32_sample(f, buffer, i, ch);
	switch(f->channels) {
	case 1:
	case 2:
	case 4:
	    for (ch = 0; ch < f->channels; ch++)
		    *(samples++) = smpte_sample[ch];
	    break;
	case 5:
	    *(samples++) = smpte_sample[2];
	    *(samples++) = smpte_sample[0];
	    *(samples++) = smpte_sample[1];
	    *(samples++) = smpte_sample[3];
	    *(samples++) = smpte_sample[4];
	    break;
	case 6:
	    *(samples++) = smpte_sample[2];
	    *(samples++) = smpte_sample[0];
	    *(samples++) = smpte_sample[1];
	    *(samples++) = smpte_sample[4];
	    *(samples++) = smpte_sample[5];
	    *(samples++) = smpte_sample[3];
	    break;
	}
    }
    
    return f->samples_left;
}
