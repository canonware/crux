/* -*- mode: c ; c-file-style: "canonware-c-style" -*-
 ******************************************************************************
 *
 * <Copyright = jasone>
 * <License>
 *
 ******************************************************************************
 *
 * Version: Crux <Version = crux>
 *
 ******************************************************************************/

#define MODCRUX_CPU_INIT() modcrux_ia32_cpu_init()

extern bool modcrux_ia32_use_sse2;

void
modcrux_ia32_cpu_init(void);
