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

#ifndef DCAENC_PRIVATE_H
#define DCAENC_PRIVATE_H

#include <stdint.h>
#include "softfloat.h"
#define DCAENC_MAX_CHANNELS 6

/* The standard allows up to 4 subsubframes per subframe,
 * but more than 2 don't fit in MPEG program stream at 1536 kbps
 */
#define DCAENC_SUBSUBFRAMES 2
#define DCAENC_LFE_SAMPLES 8
#define DCAENC_SUBBAND_SAMPLES 16

struct dcaenc_context_s {
	int samplerate_index;
	int channel_config;
	int channels;
	int fullband_channels;
	int flags;
	int bitrate_index;
	int frame_bits;
	const int32_t *band_interpolation;
	const int32_t *band_spectrum;
	int32_t pcm_history[512][DCAENC_MAX_CHANNELS];
/*	int32_t subband_history[???][???][DCAENC_MAX_CHANNELS]; */

	/* usage: subband_samples[sample][band][channel], peak[band][channel] */
	int32_t subband_samples[DCAENC_SUBBAND_SAMPLES][32][DCAENC_MAX_CHANNELS];
	int32_t quantized_samples[DCAENC_SUBBAND_SAMPLES][32][DCAENC_MAX_CHANNELS];
	int32_t peak_cb[32][DCAENC_MAX_CHANNELS];

	int32_t downsampled_lfe[DCAENC_LFE_SAMPLES];
	int32_t lfe_peak_cb;

	int32_t masking_curve_cb[DCAENC_SUBSUBFRAMES][256];
	int abits[32][DCAENC_MAX_CHANNELS];
	int nscale[32][DCAENC_MAX_CHANNELS];
	softfloat quant[32][DCAENC_MAX_CHANNELS];
	int lfe_nscale;
	softfloat lfe_quant;

	int32_t eff_masking_curve_cb[256];
	int32_t band_masking_cb[32];
	int32_t worst_quantization_noise;
	int32_t worst_noise_ever;
	int consumed_bits;

	uint32_t word;
	int wbits;
	uint8_t *output;
	int wrote;
};

#endif
