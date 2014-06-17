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

#include "dcaenc.h"
#include <stdint.h>

#include <alsa/asoundlib.h>
#include <alsa/pcm_external.h>

struct dcaplug_info {
        snd_pcm_extplug_t ext;
        int iec61937;
        dcaenc_context enc;
        int32_t pcm_buffer[6 * 512];
        int16_t dts_buffer[2 * 512];
        snd_pcm_uframes_t bufpos;
};

static inline void *area_addr(const snd_pcm_channel_area_t *area,
                              snd_pcm_uframes_t offset)
{
        unsigned int bitofs = area->first + area->step * offset;
        return (char *) area->addr + bitofs / 8;
}


static inline int32_t get_s32(const snd_pcm_channel_area_t *area,
                              snd_pcm_uframes_t offset, snd_pcm_format_t format)
{
        int32_t* sample32;
        int16_t* sample16;
        switch (format) {
        case SND_PCM_FORMAT_S32:
                sample32 = (int32_t*)area_addr(area, offset);
                return *sample32;
        case SND_PCM_FORMAT_S16:
                sample16 = (int16_t*)area_addr(area, offset);
                return *sample16 * 0x10000;
        }
        return 0;
}

static inline void put_s16(const snd_pcm_channel_area_t *area,
                              snd_pcm_uframes_t offset, int16_t value)
{
       int16_t* sample16 = (int16_t*)area_addr(area, offset);
       *sample16 = value;
}

static snd_pcm_sframes_t dcaplug_transfer(snd_pcm_extplug_t *ext,
        const snd_pcm_channel_area_t *dst_areas, snd_pcm_uframes_t dst_offset,
        const snd_pcm_channel_area_t *src_areas, snd_pcm_uframes_t src_offset,
        snd_pcm_uframes_t size)
{
        struct dcaplug_info *dcaplug = (struct dcaplug_info*)ext;

        /* Force samples to the buffer */

        snd_pcm_uframes_t remaining = 512 - dcaplug->bufpos;
        if (size > remaining)
                size = remaining;

        snd_pcm_uframes_t i;
        int channel;
        int srcbufidx = ext->channels * dcaplug->bufpos;
        int dstbufidx = 2 * dcaplug->bufpos;
        for (i = 0; i < size; i++) {
                if (ext->channels == 4) {
                        dcaplug->pcm_buffer[srcbufidx++] = get_s32(&src_areas[0], i + src_offset, ext->format);
                        dcaplug->pcm_buffer[srcbufidx++] = get_s32(&src_areas[1], i + src_offset, ext->format);
                        dcaplug->pcm_buffer[srcbufidx++] = get_s32(&src_areas[2], i + src_offset, ext->format);
                        dcaplug->pcm_buffer[srcbufidx++] = get_s32(&src_areas[3], i + src_offset, ext->format);
                } else {
                        dcaplug->pcm_buffer[srcbufidx++] = get_s32(&src_areas[4], i + src_offset, ext->format);
                        dcaplug->pcm_buffer[srcbufidx++] = get_s32(&src_areas[0], i + src_offset, ext->format);
                        dcaplug->pcm_buffer[srcbufidx++] = get_s32(&src_areas[1], i + src_offset, ext->format);
                        dcaplug->pcm_buffer[srcbufidx++] = get_s32(&src_areas[2], i + src_offset, ext->format);
                        dcaplug->pcm_buffer[srcbufidx++] = get_s32(&src_areas[3], i + src_offset, ext->format);
                        dcaplug->pcm_buffer[srcbufidx++] = get_s32(&src_areas[5], i + src_offset, ext->format);
                }

                put_s16(&dst_areas[0], i + dst_offset, dcaplug->dts_buffer[dstbufidx++]);
                put_s16(&dst_areas[1], i + dst_offset, dcaplug->dts_buffer[dstbufidx++]);
        }
        dcaplug->bufpos += size;

        if (dcaplug->bufpos == 512) {
                dcaenc_convert_s32(dcaplug->enc, dcaplug->pcm_buffer, (uint8_t*)dcaplug->dts_buffer);
                dcaplug->bufpos = 0;
        }

        return size;
}

static const int32_t zero[512 * 6];
static int dcaplug_init(snd_pcm_extplug_t *ext)
{
        struct dcaplug_info *dcaplug = (struct dcaplug_info*)ext;

        if (ext->rate != 44100 && ext->rate != 48000) {
                SNDERR("Wrong sample rate, must be 44100 or 48000 Hz");
                return -EINVAL;
        }

        if (ext->channels == 2) {
                /* TODO: passthrough? */
                SNDERR("Conversion from stereo to DTS is pointless");
                return -EINVAL;
        }

        if (ext->channels != 4 && ext->channels != 6) {
                SNDERR("Wrong number of channels");
                return -EINVAL;
        }

        if (dcaplug->iec61937) {
                dcaplug->enc = dcaenc_create(ext->rate,
                        (ext->channels == 4) ? DCAENC_CHANNELS_2FRONT_2REAR : DCAENC_CHANNELS_3FRONT_2REAR,
                        ext->rate * 503 / 16, /* same as DVD */
                        ((ext->channels == 4) ? 0 : DCAENC_FLAG_LFE) | DCAENC_FLAG_IEC_WRAP);
        } else {
                dcaplug->enc = dcaenc_create(ext->rate,
                        (ext->channels == 4) ? DCAENC_CHANNELS_2FRONT_2REAR : DCAENC_CHANNELS_3FRONT_2REAR,
                        ext->rate * 32,       /* same as S16 stereo on a CD */
                        (ext->channels == 4) ? 0 : DCAENC_FLAG_LFE);
        }

        if (!dcaplug->enc) {
                SNDERR("Failed to create DCA encoder");
                return -ENOMEM;
        }

        if (dcaenc_output_size(dcaplug->enc) != 2048) {
                SNDERR("The dcaenc library is incompatible");
                return -EINVAL;
        }

        /* Create a dummy frame of silence */
        dcaenc_convert_s32(dcaplug->enc, zero, (uint8_t*)dcaplug->dts_buffer);

        return 0;
}

static int dcaplug_close(snd_pcm_extplug_t *ext)
{
        struct dcaplug_info *dcaplug = (struct dcaplug_info*)ext;

        dcaenc_destroy(dcaplug->enc, NULL);
        dcaplug->enc = NULL;
        return 0;
}

static const snd_pcm_extplug_callback_t dcaplug_callback = {
        .transfer = dcaplug_transfer,
        .init = dcaplug_init,
        .close = dcaplug_close,
};


SND_PCM_PLUGIN_DEFINE_FUNC(dca)
{
        snd_config_iterator_t i, next;
        snd_config_t *slave = NULL;
        int iec61937 = 0;
        struct dcaplug_info *dcaplug;
        int err;

        if (stream != SND_PCM_STREAM_PLAYBACK) {
                SNDERR("dca is only for playback");
                return -EINVAL;
        }
        snd_config_for_each(i, next, conf) {
                snd_config_t *n = snd_config_iterator_entry(i);
                const char *id;
                if (snd_config_get_id(n, &id) < 0)
                        continue;
                if (strcmp(id, "comment") == 0 || strcmp(id, "type") == 0 || strcmp(id, "hint") == 0)
                        continue;
                if (strcmp(id, "slave") == 0) {
                        slave = n;
                        continue;
                }
                if (strcmp(id, "iec61937") == 0) {
                        if ((err = snd_config_get_bool(n)) < 0) {
                                SNDERR("Invalid value for %s", id);
                                return -EINVAL;
                        }
                        iec61937 = err;
                        continue;
                }
                SNDERR("Unknown field %s", id);
                return -EINVAL;
        }

        if (!slave) {
                SNDERR("No slave defined for dca");
                return -EINVAL;
        }

        dcaplug = calloc(1, sizeof(*dcaplug));
        if (dcaplug == NULL)
                return -ENOMEM;

        dcaplug->ext.version = SND_PCM_EXTPLUG_VERSION;
        dcaplug->ext.name = "DTS Coherent Acoustics encoder";
        dcaplug->ext.callback = &dcaplug_callback;
        dcaplug->ext.private_data = dcaplug;
        dcaplug->iec61937 = iec61937;

        err = snd_pcm_extplug_create(&dcaplug->ext, name, root, slave, stream, mode);
        if (err < 0) {
                dcaenc_destroy(dcaplug->enc, NULL);
                free(dcaplug);
                return err;
        }

        static const int channels[2] = {4, 6};
        static const int formats[2] = {SND_PCM_FORMAT_S32, SND_PCM_FORMAT_S16};

        snd_pcm_extplug_set_param_list(&dcaplug->ext, SND_PCM_EXTPLUG_HW_CHANNELS,
                2, channels);
        snd_pcm_extplug_set_slave_param(&dcaplug->ext, SND_PCM_EXTPLUG_HW_CHANNELS, 2);

        snd_pcm_extplug_set_param_list(&dcaplug->ext, SND_PCM_EXTPLUG_HW_FORMAT,
                2, formats);
        snd_pcm_extplug_set_slave_param(&dcaplug->ext, SND_PCM_EXTPLUG_HW_FORMAT, SND_PCM_FORMAT_S16);

        *pcmp = dcaplug->ext.pcm;
        return 0;
}

SND_PCM_PLUGIN_SYMBOL(dca);
