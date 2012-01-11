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
#include <string.h>
#include "config.h"
#include "dcaenc.h"
#include "wavfile.h"

int main(int argc, char *argv[])
{
	dcaenc_context c;
	int32_t data[512 * 6];
	uint8_t output[16384];
	wavfile * f;
	FILE * outfile;
	int bitrate;
	int wrote;
	int samples_total;
	
	static const int channel_map[6] = {DCAENC_CHANNELS_MONO, DCAENC_CHANNELS_STEREO, 0,
	DCAENC_CHANNELS_2FRONT_2REAR, DCAENC_CHANNELS_3FRONT_2REAR, DCAENC_CHANNELS_3FRONT_2REAR };

	if (argc != 4) {
	    if (argc == 2 && !strcmp(argv[1], "--version")) {
	        fprintf(stderr, PACKAGE_NAME "-" PACKAGE_VERSION "\n");
		fprintf(stderr, PACKAGE_URL "\n\n");
		fprintf(stderr, "Copyrignt (c) 2008-2012 Alexander E. Patrakov <patrakov@gmail.com>\n");
		fprintf(stderr, "License: GNU LGPL version 2.1 or later <http://gnu.org/licenses/lgpl.html>\n");
		fprintf(stderr, "This is free software: you are free to change and redistribute it.\n");
		fprintf(stderr, "There is NO WARRANTY, to the extent permitted by law.\n");
		return 0;
	    } else {
	        fprintf(stderr, "Usage: dcaenc input.wav output.dts bits_per_second\n");
	        return 1;
	    }
	}
	f = wavfile_open(argv[1]);
	if (!f) {
	    printf("Could not open or parse %s\n", argv[1]);
	    return 1;
	}
	bitrate = atoi(argv[3]);
	
	samples_total = f->samples_left;
	c = dcaenc_create(f->sample_rate, channel_map[f->channels - 1], bitrate, f->channels == 6 ? DCAENC_FLAG_LFE : 0);
	
	if (!c) {
	    printf("Wrong bitrate or sample rate\n");
	    return 1;
	}
	outfile = fopen(argv[2], "wb");
	if (!outfile) {
	    printf("Could not open %s\n", argv[2]);
	    return 1;
	}
	while (wavfile_read_s32(f, data)) {
		wrote = dcaenc_convert_s32(c, data, output);
		fwrite(output, 1, wrote, outfile);
	}
	wrote = dcaenc_destroy(c, output);
	fwrite(output, 1, wrote, outfile);
	fclose(outfile);
	wavfile_close(f);
	return 0;
}
