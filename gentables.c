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

/* Visual C++ brokenness */
#define _USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include "config.h"
#include "int_data.h"
#include "float_data.h"

#define AUBANDS 25

/* Dump the table of cosines, scaled to 0x7fffffff
 * cos_table[i] == (int32_t)(0x7fffffff * cos(M_PI * i / 1024))
 */
void print_cos_table(void)
{
	int i;
	printf("const int32_t cos_table[2048] = {");

	for (i = 0; i < 2048; i++) {
		if (i % 4 == 0)
			printf("\n   ");
		printf(" %d,", (int32_t)(0x7fffffff * cos(M_PI * i / 1024)));
	}
	printf("\n};\n\n");
}

/* Dumps the table that is used to reverse bits of small integers */
void print_bitrev_table(void)
{
	int i, j, k;
	int br[256];

	for (k = 0; k < 256; k++)
		br[k] = 0;

	for (i = 128, j = 1; i != 0; i >>= 1, j <<= 1)
		for (k = 0; k < 256; k++)
			if (k & j)
				br[k] += i;

	printf("const int bitrev[256] = {");
	for (k = 0; k < 256; k++) {
		if (k % 8 == 0)
			printf("\n   ");
		printf(" %d,", br[k]);
	}
	printf("\n};\n\n");
}

/* Transfer function of outer and middle ear, Hz -> dB */
static double hom(double f)
{
	double f1 = f / 1000;

	return -3.64 * pow(f1, -0.8)
		   + 6.8 * exp(-0.6 * (f1 - 3.4) * (f1 - 3.4))
		   - 6.0 * exp(-0.15 * (f1 - 8.7) * (f1 - 8.7))
		   - 0.0006 * (f1 * f1) * (f1 * f1);
}

/* Auditory filter center frequencies and bandwidths, in Hz.
 * The last two are made up, because there is no scientific data.
 */
double fc[AUBANDS] = {
	50, 150, 250, 350, 450, 570, 700, 840,
	1000, 1170, 1370, 1600, 1850, 2150, 2500, 2900,
	3400, 4000, 4800, 5800, 7000, 8500, 10500, 13500,
	17000
};

double erb[AUBANDS] = {
	80, 100, 100, 100, 110, 120, 140, 150,
	160, 190, 210, 240, 280, 320, 380, 450,
	550, 700, 900, 1100, 1300, 1800, 2500, 3500,
	4500
};

/* transfer function of i-th auditory filter, in dB */
static double gammafilter(int i, double f)
{
	double h = (f - fc[i]) / erb[i];
	h = 1 + h * h;
	h = 1 / (h * h);
	return 20 * log10(h);
}

/* Dumps the auditory filter */
void print_auf(void)
{
	int i, j, k;

	printf("const int32_t auf[9][AUBANDS][256] = {\n   ");
	for (i = 0; i < 9; i++) {
		printf(" { /* Sample rate: %d Hz */\n       ", sample_rates[i]);
		for (j = 0; j < AUBANDS; j++) {
			printf(" {");
			for (k = 0; k < 256; k++) {
				double freq = sample_rates[i] * (k + 0.5) / 512;
				if (k % 8 == 0)
					printf("\n           ");
				printf(" %d,", (int)(10 * (hom(freq) + gammafilter(j, freq))));
			}
			printf("\n        },");
		}
		printf("\n    },");
	}
	printf("\n};\n\n");
}

/* Dumps the exponent table,
 * used to convert centibels to 32-bit integer samples
 */
void print_cb_table(void)
{
	int i;
	printf("const int32_t cb_to_level[2048] = {");
	for (i = 0; i < 2048; i++) {
		if (i % 8 == 0)
			printf("\n   ");
		printf(" %d,", (int)(0x7fffffff * pow(10, -0.005 * i)));
	}
	printf("\n};\n\n");
}

/* Dumps the table that is used for addition of centibels */
void print_cb_add_table(void)
{
	int i;
	printf("const int32_t cb_to_add[256] = {");
	for (i = 0; i < 256; i++) {
		double add = 1 + pow(10, -0.01 * i);
		add = 100 * log10(add);
		if (i % 8 == 0)
			printf("\n   ");
		printf(" %d,", (int)add);
	}
	printf("\n};\n\n");
}

/* Dumps the integer version of the LFE FIR table */
void print_lfe_fir(void)
{
	int i;
	printf("const int32_t lfe_fir[512] = {");
	
	for (i = 0; i < 512; i++) {
		if (i % 8 == 0)
			printf("\n   ");
		printf(" %d,", (int32_t)(0x01ffffff * lfe_fir_64[i]));
	}
	
	printf("\n};\n\n");

}

/* Dumps the integer version of the subband FIR table and its spectrum */
void print_subband_fir(void)
{
	int i, j, k;
	double spectrum[2][8];
	
	for (j = 0; j < 2; j++) {

		for (k = 0; k < 8; k++) {
			double accum = 0;
			for (i = 0; i < 512; i++) {
				double reconst = reconstruction_fir[j][i] * ((i & 64) ? (-1) : 1);
				accum += reconst * cos(2 * M_PI * (i + 0.5 - 256) * (k + 0.5) / 512);
			}
			spectrum[j][k] = accum;
		}
	}

	printf("const int32_t band_interpolation[2][512] = {\n");
	for (j = 0; j < 2; j++) {
		printf("  {");
		for (i = 0; i < 512; i++) {
			if (i % 8 == 0)
				printf("\n   ");
			printf(" %d,", (int32_t) (0x1000000000ULL * reconstruction_fir[j][i]));
		}
		printf("\n  },\n");
	}
	printf("};\n\n");

	printf("const int32_t band_spectrum[2][8] = {\n");
	for (j = 0; j < 2; j++) {
		printf("\t{ ");
		for (i = 0; i < 8; i++)
			printf("%d, ", (int32_t) (200 * log10(spectrum[j][i])));
		printf("},\n");
	}
	printf("};\n\n");

}

int main(int argc, char *argv[])
{
	printf("/* GENERATED FILE, DO NOT EDIT */\n\n");
	printf("#include \"config.h\"\n");
	printf("#include \"math_tables.h\"\n\n");
	print_lfe_fir();
	print_subband_fir();
	print_cos_table();
	print_bitrev_table();
	print_cb_table();
	print_cb_add_table();
	print_auf();
	
	return 0;
}
