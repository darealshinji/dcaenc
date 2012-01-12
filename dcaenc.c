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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "config.h"
#include "dcaenc.h"
#include "dcaenc_private.h"
#include "int_data.h"
#include "math_tables.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define div_round_up(a, b) (((a) + (b) - 1) / (b))
#define round_up(a, b) ((((a) + (b) - 1) / (b)) * (b))

int dcaenc_channel_config_to_count(int channel_config)
{
	if (channel_config > DCAENC_CHANNELS_3FRONT_2REAR)
		return -1;

	if (channel_config < 0)
		return -1;

	return channels_table[channel_config];
}

dcaenc_context dcaenc_create(int sample_rate, int channel_config,
                             int approx_bitrate, int flags)
{
	int i, frame_bits, bit_step, fir, useful_bitrate, min_frame_bits;
	dcaenc_context result;

	i = 0;
	while (i < 9 && sample_rate != sample_rates[i])
		i++;
	if (i == 9)
		return NULL;

	if (approx_bitrate < 32000)
		return NULL;
	if (approx_bitrate > 6144000)
		return NULL;

	if (channel_config < 0)
		return NULL;
	if (channel_config > DCAENC_CHANNELS_3FRONT_2REAR)
		return NULL;


	if (flags & DCAENC_FLAG_28BIT)
		useful_bitrate = div_round_up(approx_bitrate * 7, 8);
	else
		useful_bitrate = approx_bitrate;

	frame_bits = div_round_up(useful_bitrate * 512, sample_rate);

	bit_step = (flags & DCAENC_FLAG_28BIT) ? (7 * 32) : 32;
	fir = (flags & DCAENC_FLAG_PERFECT_QMF) ? 0 : 1;
	/* Round frame_bits up to the next permitted value */
	frame_bits = round_up(frame_bits, bit_step);
	
	min_frame_bits = 132 + (493 + 28 * 32) * channels_table[channel_config];
	if (flags & DCAENC_FLAG_LFE)
		min_frame_bits += 72;

	if (frame_bits < min_frame_bits)
		return NULL;

	if (frame_bits > 131072)
		return NULL;

	if ((flags & DCAENC_FLAG_IEC_WRAP) && frame_bits > 16320)
		return NULL;

	result = (dcaenc_context)calloc(1, sizeof(struct dcaenc_context_s));
	if (!result)
		return NULL;

	result->samplerate_index = i;
	result->frame_bits = frame_bits;
	result->channel_config = channel_config;
	result->flags = flags;
	result->channels = result->fullband_channels = channels_table[channel_config];
	if (flags & DCAENC_FLAG_LFE)
		result->channels++;

	
	for (i = 0; target_bitrate_table[i] < approx_bitrate; i++)
		;
	result->bitrate_index = i;
	result->band_interpolation = band_interpolation[fir];
	result->band_spectrum = band_spectrum[fir];
	result->worst_quantization_noise = -2047;
	result->worst_noise_ever = -2047;
	return result;
}

int dcaenc_input_size(dcaenc_context c)
{
	return 512;
}

int dcaenc_output_size(dcaenc_context c)
{
	if (c->flags & DCAENC_FLAG_IEC_WRAP)
		return 2048;
	return c->frame_bits / ((c->flags & DCAENC_FLAG_28BIT) ? 7 : 8);
}

static inline const int32_t *pcm_sample(dcaenc_context c,
                                        const int32_t *container,
                                        int sample, int channel)
{
	return &container[sample * c->channels + channel];
}

static inline int32_t half32(int32_t a)
{
	return (a + 1) >> 1;
}

static inline int32_t mul32(int32_t a, int32_t b)
{
	int64_t r = (int64_t)a * b + 0x80000000ULL;
	return r >> 32;
}

static void dcaenc_subband_transform(dcaenc_context c, const int32_t *input)
{
	int ch, subs, i, k, j;

	for (ch = 0; ch < c->fullband_channels; ch++) {
		/* History is copied because it is also needed for PSY */
		int32_t hist[512];
		int hist_start = 0;

		for (i = 0; i < 512; i++)
			hist[i] = c->pcm_history[i][ch];

		for (subs = 0; subs < DCAENC_SUBBAND_SAMPLES; subs++) {
			int32_t accum[64];
			int32_t resp;
			int band;

			/* Calculate the convolutions at once */
			for (i = 0; i < 64; i++)
				accum[i] = 0;

			for (k = 0, i = hist_start, j = 0;
			     i < 512; k = (k + 1) & 63, i++, j++)
				accum[k] += mul32(hist[i], c->band_interpolation[j]);
			for (i = 0; i < hist_start; k = (k + 1) & 63, i++, j++)
				accum[k] += mul32(hist[i], c->band_interpolation[j]);

			for (k = 16; k < 32; k++)
				accum[k] = accum[k] - accum[31 - k];
			for (k = 32; k < 48; k++)
				accum[k] = accum[k] + accum[95 - k];

			for (band = 0; band < 32; band++) {
				resp = 0;
				for (i = 16; i < 48; i++) {
					int s = (2 * band + 1) * (2 * (i + 16) + 1);
					resp += mul32(accum[i], cos_t(s << 3)) >> 3;
				}

				c->subband_samples[subs][band][ch] = ((band + 1) & 2) ?
					(-resp) : resp;
			}

			/* Copy in 32 new samples from input */
			for (i = 0; i < 32; i++)
				hist[i + hist_start] = *pcm_sample(c, input, subs * 32 + i, ch);

			hist_start = (hist_start + 32) & 511;
		}
	}
}

static void dcaenc_lfe_downsample(dcaenc_context c, const int32_t *input)
{
	/* FIXME: make 128x LFE downsampling possible */
	int i, j, lfes;

	int32_t hist[512];
	int32_t accum;
	int hist_start = 0;

	for (i = 0; i < 512; i++)
		hist[i] = c->pcm_history[i][c->channels - 1];

	for (lfes = 0; lfes < DCAENC_LFE_SAMPLES; lfes++) {
		/* Calculate the convolution */
		accum = 0;

		for (i = hist_start, j = 0; i < 512; i++, j++)
			accum += mul32(hist[i], lfe_fir[j]);
		for (i = 0; i < hist_start; i++, j++)
			accum += mul32(hist[i], lfe_fir[j]);

		c->downsampled_lfe[lfes] = accum;

		/* Copy in 64 new samples from input */
		for (i = 0; i < 64; i++)
			hist[i + hist_start] = *pcm_sample(c, input,
			                                   lfes * 64 + i, c->channels - 1);

		hist_start = (hist_start + 64) & 511;
	}
}

typedef struct {
	int32_t re;
	int32_t im;
} cplx32;

static void fft(const int32_t in[2 * 256], cplx32 out[256])
{
	cplx32 buf[256], rin[256], rout[256];
	int i, j, k, l;

	/* do two transforms in parallel */
	for (i = 0; i < 256; i++) {
		/* Apply the Hann window */
		rin[i].re = mul32(in[2 * i], 0x3fffffff - (cos_t(8 * i + 2) >> 1));
		rin[i].im = mul32(in[2 * i + 1], 0x3fffffff - (cos_t(8 * i + 6) >> 1));
	}
	/* pre-rotation */
	for (i = 0; i < 256; i++) {
		buf[i].re = mul32(cos_t(4 * i + 2), rin[i].re)
			- mul32(sin_t(4 * i + 2), rin[i].im);
		buf[i].im = mul32(cos_t(4 * i + 2), rin[i].im)
			+ mul32(sin_t(4 * i + 2), rin[i].re);
	}

	for (j = 256, l = 1; j != 1; j >>= 1, l <<= 1)
		for (k = 0; k < 256; k += j)
			for (i = k; i < k + j / 2; i++) {
				cplx32 sum, diff;

				int t = 8 * l * i;

				sum.re = buf[i].re + buf[i + j / 2].re;
				sum.im = buf[i].im + buf[i + j / 2].im;

				diff.re = buf[i].re - buf[i + j / 2].re;
				diff.im = buf[i].im - buf[i + j / 2].im;

				buf[i].re = half32(sum.re);
				buf[i].im = half32(sum.im);

				buf[i + j / 2].re = mul32(diff.re, cos_t(t))
					- mul32(diff.im, sin_t(t));
				buf[i + j / 2].im = mul32(diff.im, cos_t(t))
					+ mul32(diff.re, sin_t(t));
			}
	/* post-rotation */
	for (i = 0; i < 256; i++) {
		int b = bitrev[i];
		rout[i].re = mul32(buf[b].re, cos_t(4 * i))
			- mul32(buf[b].im, sin_t(4 * i));
		rout[i].im = mul32(buf[b].im, cos_t(4 * i))
			+ mul32(buf[b].re, sin_t(4 * i));
	}

	for (i = 0; i < 256; i++) {
		/* separate the results of the two transforms */
		cplx32 o1, o2;

		o1.re = rout[i].re - rout[255 - i].re;
		o1.im = rout[i].im + rout[255 - i].im;

		o2.re = rout[i].im - rout[255 - i].im;
		o2.im = -rout[i].re - rout[255 - i].re;

		/* combine them into one long transform */
		out[i].re = mul32(o1.re + o2.re, cos_t(2 * i + 1))
			+ mul32(o1.im - o2.im, sin_t(2 * i + 1));
		out[i].im = mul32(o1.im + o2.im, cos_t(2 * i + 1))
			+ mul32(-o1.re + o2.re, sin_t(2 * i + 1));
	}
}

/* in: anything, result: between -2047 and 0 */
static int32_t get_cb(int32_t in)
{
	int i, res;
	res = 0;
	if (in < 0)
		in = -in;
	for (i = 1024; i > 0; i >>= 1) {
		if (cb_to_level[i + res] >= in)
			res += i;
	}
	return -res;
}

/* accepts any input, gives no guarantees about the range of the result,
   except that result >= a && result >= b */
static int32_t add_cb(int32_t a, int32_t b)
{
	if (a < b) {
		int32_t tmp = a;
		a = b;
		b = tmp;
	}
	if (a - b >= 256)
		return a;

	return a + cb_to_add[a - b];
}

/* accepts anything, out_cb[i] can only grow */
static void adjust_jnd(int samplerate_index,
                       const int32_t in[512], int32_t out_cb[256])
{
	int i, j;
	int32_t power[256];
	cplx32 out[256];
	int32_t out_cb_unnorm[256];
	int32_t denom;
	const int32_t ca_cb = -1114;
	const int32_t cs_cb = 928;

	fft(in, out);

	for (j = 0; j < 256; j++) {
		power[j] = add_cb(get_cb(out[j].re), get_cb(out[j].im));
		out_cb_unnorm[j] = -2047;	/* and can only grow */
	}

	for (i = 0; i < AUBANDS; i++) {
		denom = ca_cb;	/* and can only grow */
		for (j = 0; j < 256; j++)
			denom = add_cb(denom, power[j] + auf[samplerate_index][i][j]);
		for (j = 0; j < 256; j++)
			out_cb_unnorm[j] = add_cb(out_cb_unnorm[j],
			                          -denom + auf[samplerate_index][i][j]);
	}

	for (j = 0; j < 256; j++)
		out_cb[j] = add_cb(out_cb[j], -out_cb_unnorm[j] - ca_cb - cs_cb);
}


typedef void (*walk_band_t)(dcaenc_context c, int band1, int band2, int f,
                            int32_t spectrum1, int32_t spectrum2, int channel,
                            int32_t * arg);

static void walk_band_low(dcaenc_context c, int band, int channel,
                          walk_band_t walk, int32_t * arg)
{
	int f;
	if (band == 0) {
		for (f = 0; f < 4; f++)
			walk(c, 0, 0, f, 0, -2047, channel, arg);
	} else {
		for (f = 0; f < 8; f++)
			walk(c, band, band - 1, 8 * band - 4 + f,
			     c->band_spectrum[7 - f], c->band_spectrum[f], channel, arg);
	}
}

static void walk_band_high(dcaenc_context c, int band, int channel,
                           walk_band_t walk, int32_t * arg)
{
	int f;
	if (band == 31) {
		for (f = 0; f < 4; f++)
			walk(c, 31, 31, 256 - 4 + f, 0, -2047, channel, arg);
	} else {
		for (f = 0; f < 8; f++)
			walk(c, band, band + 1, 8 * band + 4 + f,
			     c->band_spectrum[f], c->band_spectrum[7 - f], channel, arg);
	}
}

static void walk_whole_spectrum(dcaenc_context c, int channel,
                                walk_band_t walk, int32_t * arg)
{
	int band;
	for (band = 0; band < 32; band++)
		walk_band_low(c, band, channel, walk, arg);
	walk_band_high(c, 31, channel, walk, arg);
}

static void update_band_masking(dcaenc_context c, int band1, int band2,
							int f, int32_t spectrum1, int32_t spectrum2,
							int channel, int32_t * arg)
{
	int32_t value = c->eff_masking_curve_cb[f] - spectrum1;
	if (value < c->band_masking_cb[band1])
		c->band_masking_cb[band1] = value;
}

static void dcaenc_calc_masking(dcaenc_context c, const int32_t *input)
{
	int i, k, band, ch, ssf;
	int32_t data[512];

	for (i = 0; i < 256; i++)
		for (ssf = 0; ssf < DCAENC_SUBSUBFRAMES; ssf++)
			c->masking_curve_cb[ssf][i] = -2047;

	for (ssf = 0; ssf < DCAENC_SUBSUBFRAMES; ssf++)
		for (ch = 0; ch < c->fullband_channels; ch++) {
			for (i = 0, k = 128 + 256 * ssf; k < 512; i++, k++)
				data[i] = c->pcm_history[k][ch];
			for (k -= 512; i < 512; i++, k++)
				data[i] = *pcm_sample(c, input, k, ch);
			adjust_jnd(c->samplerate_index, data, c->masking_curve_cb[ssf]);
		}
	for (i = 0; i < 256; i++) {
		int32_t m = 2048;
		for (ssf = 0; ssf < DCAENC_SUBSUBFRAMES; ssf++)
			if (c->masking_curve_cb[ssf][i] < m)
				m = c->masking_curve_cb[ssf][i];
		c->eff_masking_curve_cb[i] = m;
	}
	
	for (band = 0; band < 32; band++) {
		c->band_masking_cb[band] = 2048;
		walk_band_low(c, band, 0, update_band_masking, NULL);
		walk_band_high(c, band, 0, update_band_masking, NULL);
	}
}

static void dcaenc_find_peaks(dcaenc_context c)
{
	int band, ch;

	for (band = 0; band < 32; band++)
		for (ch = 0; ch < c->fullband_channels; ch++) {
			int sample;
			int32_t m = 0;
			for (sample = 0; sample < DCAENC_SUBBAND_SAMPLES; sample++) {
				int32_t s = abs(c->subband_samples[sample][band][ch]);
				if (m < s)
					m = s;
			}
			c->peak_cb[band][ch] = get_cb(m);
		}

	if (c->flags & DCAENC_FLAG_LFE) {
		int sample;
		int32_t m = 0;
		for (sample = 0; sample < DCAENC_LFE_SAMPLES; sample++)
			if (m < abs(c->downsampled_lfe[sample]))
				m = abs(c->downsampled_lfe[sample]);
		c->lfe_peak_cb = get_cb(m);
	}
}

static const int snr_fudge = 128;
static const int USED_1ABITS = 1;
static const int USED_NABITS = 2;
static const int USED_26ABITS = 4;

static int init_quantization_noise(dcaenc_context c, int noise)
{
	int ch, band;
	int ret = 0;
	
	c->consumed_bits = 132 + 493 * c->fullband_channels;
	if (c->flags & DCAENC_FLAG_LFE)
		c->consumed_bits += 72;

	if (c->flags & DCAENC_FLAG_IEC_WRAP)
		c->consumed_bits += (c->flags & DCAENC_FLAG_28BIT) ? 56 : 64;

	/* attempt to guess the bit distribution based on the prevoius frame */
	for (ch = 0; ch < c->fullband_channels; ch++) {
		for (band = 0; band < 32; band++) {
			int snr_cb = c->peak_cb[band][ch]
				- c->band_masking_cb[band]
				- noise;

			if (snr_cb >= 1312) {
				c->abits[band][ch] = 26;
				ret |= USED_26ABITS;
			} else if (snr_cb >= 222) {
				c->abits[band][ch] = 8 + mul32(snr_cb - 222, 69000000);
				ret |= USED_NABITS;
			} else if (snr_cb >= 0) {
				c->abits[band][ch] = 2 + mul32(snr_cb, 106000000);
				ret |= USED_NABITS;
			} else {
				c->abits[band][ch] = 1;
				ret |= USED_1ABITS;
			}
		}
	}

	for (band = 0; band < 32; band++)
		for (ch = 0; ch < c->fullband_channels; ch++) {
			c->consumed_bits += bit_consumption[c->abits[band][ch]];
	}

	return ret;
}

static void dcaenc_assign_bits(dcaenc_context c)
{
	/* Find the bounds where the binary search should work */
	int low, high, down;
	int used_abits = 0;
	init_quantization_noise(c, c->worst_quantization_noise);
	low = high = c->worst_quantization_noise;
	if (c->consumed_bits > c->frame_bits) {
		while (c->consumed_bits > c->frame_bits) {
			assert(("Too low bitrate should have been rejected in dcaenc_create", used_abits != USED_1ABITS));
			low = high;
			high += snr_fudge;
			used_abits = init_quantization_noise(c, high);
		}
	} else {
		while (c->consumed_bits <= c->frame_bits) {
			high = low;
			if (used_abits == USED_26ABITS)
				goto out; /* The requested bitrate is too high, pad with zeros */
			low -= snr_fudge;
			used_abits = init_quantization_noise(c, low);
		}
	}

	/* Now do a binary search between low and high to see what fits */
	for (down = snr_fudge >> 1; down; down >>= 1) {
		init_quantization_noise(c, high - down);
		if (c->consumed_bits <= c->frame_bits)
			high -= down;
	}
	init_quantization_noise(c, high);
out:
	c->worst_quantization_noise = high;
	if (high > c->worst_noise_ever)
		c->worst_noise_ever = high;
}

static void dcaenc_shift_history(dcaenc_context c, const int32_t *input)
{
	int k, ch;

	for (k = 0; k < 512; k++)
		for (ch = 0; ch < c->channels; ch++)
			c->pcm_history[k][ch] = *pcm_sample(c, input, k, ch);
}

static void bitstream_init(dcaenc_context c, uint8_t *output)
{
	c->word = 0;
	c->wbits = 0;
	c->output = output;
	c->wrote = 0;
}

static void bitstream_put(dcaenc_context c, uint32_t bits, int nbits)
{
	int max_bits = (c->flags & DCAENC_FLAG_28BIT) ? 28 : 32;
	assert(bits < (1 << nbits));
	c->wrote += nbits;
	bits &= ~(0xffffffff << nbits);
	if (nbits + c->wbits >= max_bits) {
		uint8_t b1, b2, b3, b4;
		c->word |= bits >> (nbits + c->wbits - max_bits);
		if (c->flags & DCAENC_FLAG_28BIT) {
			b1 = (c->word >> 22) & 0x3f;
			if (b1 & 0x20)
				b1 |= 0xc0;
			b2 = (c->word >> 14) & 0xff;
			b3 = (c->word >> 8) & 0x3f;
			if (b3 & 0x20)
				b3 |= 0xc0;
			b4 = (c->word) & 0xff;
		} else {
			b1 = (c->word >> 24) & 0xff;
			b2 = (c->word >> 16) & 0xff;
			b3 = (c->word >> 8) & 0xff;
			b4 = (c->word) & 0xff;
		}
		if (c->flags & DCAENC_FLAG_BIGENDIAN) {
			*(c->output++) = b1;
			*(c->output++) = b2;
			*(c->output++) = b3;
			*(c->output++) = b4;
		} else {
			*(c->output++) = b2;
			*(c->output++) = b1;
			*(c->output++) = b4;
			*(c->output++) = b3;
		}
		c->wbits = nbits + c->wbits - max_bits;
		if (c->wbits)
			c->word = (bits << (32 - c->wbits)) >> (32 - max_bits);
		else
			c->word = 0;
	} else {
		c->word |= bits << (max_bits - nbits - c->wbits);
		c->wbits += nbits;
	}
}

static void bitstream_flush(dcaenc_context c)
{
	int max_bits = (c->flags & DCAENC_FLAG_28BIT) ? 28 : 32;
	bitstream_put(c, 0, max_bits - 1);
}

static int32_t dcaenc_quantize_value(int32_t value, softfloat quant)
{
	int32_t offset = 1 << (quant.e - 1);
	value = mul32(value, quant.m) + offset;
	value = value >> quant.e;
	return value;
}

static int32_t dcaenc_quantize(dcaenc_context c, int sample, int band, int ch)
{
	int32_t result = dcaenc_quantize_value(c->subband_samples[sample][band][ch],
					       c->quant[band][ch]);

	assert(result <= (quant_levels[c->abits[band][ch]] - 1) / 2);
	assert(result >= -(quant_levels[c->abits[band][ch]] / 2));
	return result;
}

static int dcaenc_calc_one_scale(int32_t peak_cb, int abits, softfloat *quant)
{
	int32_t peak;
	int our_nscale, try_remove;
	softfloat our_quant;

	assert(peak_cb <= 0);
	assert(peak_cb >= -2047);

	our_nscale = 127;
	peak = cb_to_level[-peak_cb];

	for (try_remove = 64; try_remove > 0; try_remove >>= 1) {
		if (scalefactor_inv[our_nscale - try_remove].e + stepsize_inv[abits].e <= 17)
			continue;
		our_quant.m = mul32(scalefactor_inv[our_nscale - try_remove].m, stepsize_inv[abits].m);
		our_quant.e = scalefactor_inv[our_nscale - try_remove].e + stepsize_inv[abits].e - 17;
		if ((quant_levels[abits] - 1) / 2 < dcaenc_quantize_value(peak, our_quant))
			continue;
		our_nscale -= try_remove;
	}
	if (our_nscale >= 125)
		our_nscale = 124;
	quant->m = mul32(scalefactor_inv[our_nscale].m, stepsize_inv[abits].m);
	quant->e = scalefactor_inv[our_nscale].e + stepsize_inv[abits].e - 17;
	assert((quant_levels[abits] - 1) / 2 >= dcaenc_quantize_value(peak, *quant));
	return our_nscale;
}

static void dcaenc_calc_scales(dcaenc_context c)
{
	int band, ch;

	for (band = 0; band < 32; band++)
		for (ch = 0; ch < c->fullband_channels; ch++)
			c->nscale[band][ch] = dcaenc_calc_one_scale(c->peak_cb[band][ch], c->abits[band][ch], &c->quant[band][ch]);

	if (c->flags & DCAENC_FLAG_LFE)
		c->lfe_nscale = dcaenc_calc_one_scale(c->lfe_peak_cb, 11, &c->lfe_quant);
}

static void dcaenc_quantize_all(dcaenc_context c)
{
	int sample, band, ch;
	for (sample = 0; sample < DCAENC_SUBBAND_SAMPLES; sample++)
		for (band = 0; band < 32; band++)
			for (ch = 0; ch < c->fullband_channels; ch++)
				c->quantized_samples[sample][band][ch] = dcaenc_quantize(c, sample, band, ch);
}

static void put_frame_header(dcaenc_context c, int normal)
{
	/* SYNC */
	bitstream_put(c, 0x7ffe, 16);	//16
	bitstream_put(c, 0x8001, 16);	//32

	/* Frame type */
	bitstream_put(c, normal, 1);	//33

	/* Deficit sample count: none */
	bitstream_put(c, 31, 5);	//38

	/* CRC is not present */
	bitstream_put(c, 0, 1);	//39

	/* Number of PCM sample blocks: 16 */
	bitstream_put(c, 15, 7);	//46

	/* Primary frame byte size */
	bitstream_put(c, c->frame_bits / 8 - 1, 14);	//60

	/* Audio channel arrangement */
	bitstream_put(c, c->channel_config, 6);	//66

	/* Core audio sampling frequency */
	bitstream_put(c, bitstream_sfreq[c->samplerate_index], 4);	//70

	/* Transmission bit rate */
	bitstream_put(c, c->bitrate_index, 5);	//75

	/* Embedded down mix: disabled */
	bitstream_put(c, 0, 1);	//76

	/* Embedded dynamic range flag: not present */
	bitstream_put(c, 0, 1);	//77

	/* Embedded time stamp flag: not present */
	bitstream_put(c, 0, 1);	//78

	/* Auxiliary data flag: not present */
	bitstream_put(c, 0, 1);	//79

	/* HDCD source: no */
	bitstream_put(c, 0, 1);	//80

	/* Extension audio ID: N/A */
	bitstream_put(c, 0, 3);	//83

	/* Extended audio data: not present */
	bitstream_put(c, 0, 1);	//84

	/* Audio sync word insertion flag: after each sub-frame */
	bitstream_put(c, 0, 1);	//85

	/* Low frequency effects flag: not present or 64x subsampling */
	bitstream_put(c, (c->flags & DCAENC_FLAG_LFE) ? 2 : 0, 2);	//87

	/* Predictor history switch flag: on */
	bitstream_put(c, 1, 1);	//88

	/* No CRC */
	/* Multirate interpolator switch */
	bitstream_put(c, (c->flags & DCAENC_FLAG_PERFECT_QMF) ? 1 : 0, 1);	//89

	/* Encoder software revision: 7 */
	bitstream_put(c, 7, 4);	//93

	/* Copy history: 0 */
	bitstream_put(c, 0, 2);	//95

	/* Source PCM resolution: 16 bits (even though the library accepts only 32-bit samples), not DTS ES */
	bitstream_put(c, 0, 3);	//98

	/* Front sum/difference coding: no */
	bitstream_put(c, 0, 1);	//99

	/* Surrounds sum/difference coding: no */
	bitstream_put(c, 0, 1);	//100

	/* Dialog normalization: 0 dB */
	bitstream_put(c, 0, 4);	//104
}

static void put_primary_audio_header(dcaenc_context c)
{
	int ch;
	/* Number of subframes: 1 */
	bitstream_put(c, 0, 4);	//108

	/* Number of primary audio channels */
	bitstream_put(c, c->fullband_channels - 1, 3);	//111

	/* Subband activity count: 32 in each channel */
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 30, 5);	//111 + 5 * c

	/* High frequency VQ start subband: 32 (i.e., no VQ) */
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 31, 5);	// 111 + 10 * c

	/* Joint intensity coding index: 0 */
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 0, 3);	// 111 + 13 * c

	/* Transient mode codebook: A4 (arbitrary) */
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 0, 2);	// 111 + 15 * c

	/* Scale factor code book: 7 bit linear, 7-bit sqrt table (for each channel) */
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 6, 3);	// 111 + 18 * c

	/* Bit allocation quantizer select: linear 5-bit */
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 6, 3);	// 111 + 21 * c

	/* Quantization index codebook select: dummy data
	   to avoid transmission of scale factor adjustment */
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 1, 1);	// 111 + 22 * c
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 3, 2);	// 111 + 24 * c
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 3, 2);	// 111 + 26 * c
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 3, 2);	// 111 + 28 * c
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 3, 2);	// 111 + 30 * c
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 7, 3);	// 111 + 33 * c
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 7, 3);	// 111 + 36 * c
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 7, 3);	// 111 + 39 * c
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 7, 3);	// 111 + 42 * c
	for (ch = 0; ch < c->fullband_channels; ch++)
		bitstream_put(c, 7, 3);	// 111 + 45 * c

	/* Scale factor adjustment index: not transmitted */
	/* Audio header CRC check word: not transmitted */
}

static void put_subframe_samples(dcaenc_context c, int ss, int band, int ch)
{
	if (c->abits[band][ch] <= 7) {
		int sum, i, j;
		for (i = 0; i < 8; i += 4) {
			sum = 0;
			for (j = 3; j >= 0; j--) {
				sum *= quant_levels[c->abits[band][ch]];
				sum += c->quantized_samples[ss * 8 + i + j][band][ch];
				sum += (quant_levels[c->abits[band][ch]] - 1) / 2;
			}
			bitstream_put(c, sum, bit_consumption[c->abits[band][ch]] / 4);
		}
	} else {
		int i;
		for (i = 0; i < 8; i++) {
			int bits = bit_consumption[c->abits[band][ch]] / 16;
			int32_t mask = (1 << bits) - 1;
			bitstream_put(c,
				c->quantized_samples[ss * 8 + i][band][ch] & mask,
				bits);
		}
	}
}

static void put_subframe(dcaenc_context c)
{
	int ch, band, ss;

	/* Subsubframe count: 2 */
	bitstream_put(c, 1, 2);	// 113 + 45 * c

	/* Partial subsubframe sample count: dummy */
	bitstream_put(c, 0, 3);	// 116 + 45 * c

	/* Prediction mode: no ADPCM, in each channel and subband */
	for (ch = 0; ch < c->fullband_channels; ch++)
		for (band = 0; band < 32; band++)
			bitstream_put(c, 0, 1);	// 116 + 77 * c

	/* Prediction VQ addres: not transmitted */
	/* Bit allocation index for each channel and subband */
	for (ch = 0; ch < c->fullband_channels; ch++)
		for (band = 0; band < 32; band++)
			bitstream_put(c, c->abits[band][ch], 5);	// 116 + 237 * c

	/* Transition mode: none for each channel and subband */
	for (ch = 0; ch < c->fullband_channels; ch++)
		for (band = 0; band < 32; band++)
			bitstream_put(c, 0, 1);	/* according to Huffman codebook A4 */// 116 + 269 * c

	/* Scale factors */
	for (ch = 0; ch < c->fullband_channels; ch++)
		for (band = 0; band < 32; band++)
			bitstream_put(c, c->nscale[band][ch], 7);	// 116 + 493 * c

	/* Joint subband scale factor codebook select: not transmitted */
	/* Scale factors for joint subband coding: not transmitted */
	/* Stereo down-mix coefficients: not transmitted */
	/* Dynamic range coefficient: not transmitted */
	/* Side information CRC check word: not transmitted */
	/* VQ encoded high frequency subbands: not transmitted */
	/* LFE data: 8 samples and scalefactor */
	if (c->flags & DCAENC_FLAG_LFE) {
		for (ss = 0; ss < DCAENC_LFE_SAMPLES; ss++)
			bitstream_put(c, dcaenc_quantize_value(c->downsampled_lfe[ss], c->lfe_quant) & 0xff, 8);
		bitstream_put(c, c->lfe_nscale, 8);
	}
	/* Audio data: 2 subsubframes */
	for (ss = 0; ss < 2; ss++)
		for (ch = 0; ch < c->fullband_channels; ch++)
			for (band = 0; band < 32; band++)
				put_subframe_samples(c, ss, band, ch);
	/* DSYNC */
	bitstream_put(c, 0xffff, 16);
}

static void dump_bits(dcaenc_context c)
{
	int band, ch;
	printf("bits:");
	for (band = 0; band < 32; band++) {
		printf("(");
		for (ch = 0; ch < c->fullband_channels; ch++) {
			printf("%02d", c->abits[band][ch]);
			if (ch < c->fullband_channels - 1)
				printf(" ");
		}
		printf(")");
	}
	printf("\n");
}

static void dump_masking_curve(dcaenc_context c)
{
	int f;
	
	printf("\nWorst noise: %d (%d)\n", c->worst_quantization_noise, c->worst_noise_ever);
	printf("masking curve: ");
	for (f = 0; f < 256; f++)
		printf(" %d", c->eff_masking_curve_cb[f]);
	printf("\n");
}

static int dcaenc_convert_frame(dcaenc_context c, const int32_t *input, uint8_t *output, int normal)
{
	int i;

	int l = dcaenc_output_size(c);
	for (i = 0; i < l; i++)
		output[i] = 0;

	if (c->flags & DCAENC_FLAG_IEC_WRAP) {
		if (c->flags & DCAENC_FLAG_BIGENDIAN) {
			*(output++) = 0xf8; *(output++) = 0x72;
			*(output++) = 0x4e; *(output++) = 0x1f;
			*(output++) = 0x00; *(output++) = 0x0b;
			*(output++) = c->frame_bits >> 8;
			*(output++) = c->frame_bits & 0xff;
		} else {
			*(output++) = 0x72; *(output++) = 0xf8;
			*(output++) = 0x1f; *(output++) = 0x4e;
			*(output++) = 0x0b; *(output++) = 0x00;
			*(output++) = c->frame_bits & 0xff;
			*(output++) = c->frame_bits >> 8;
		}
	}

	bitstream_init(c, output);
	dcaenc_subband_transform(c, input);
	if (c->flags & DCAENC_FLAG_LFE)
		dcaenc_lfe_downsample(c, input);
	dcaenc_calc_masking(c, input);
	dcaenc_find_peaks(c);

	dcaenc_assign_bits(c);
	dcaenc_calc_scales(c);
	dcaenc_quantize_all(c);


	dcaenc_shift_history(c, input);
	put_frame_header(c, normal);
	put_primary_audio_header(c);
	put_subframe(c);
	bitstream_flush(c);

	/* FIXME: this will be invalid for VBR. We should return the actual
	 * number of bytes that we wrote.
	 */
	// dump_bits(c);
	return dcaenc_output_size(c);
}

int dcaenc_convert_s32(dcaenc_context c, const int32_t *input, uint8_t *output)
{
	return dcaenc_convert_frame(c, input, output, 1);
}

static int32_t zero[1024 * DCAENC_MAX_CHANNELS];

int dcaenc_destroy(dcaenc_context c, uint8_t *output)
{
	int retval = output ? dcaenc_convert_frame(c, zero, output, 0) : 0;
	free(c);
	return retval;
}
