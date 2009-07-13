#ifndef CxIa32_h
#  define CxIa32_h

#  include "Cx.h"

#  ifdef CxmCpuIa32
#    define CxmCpuInit() CxIa32CpuInit()

extern bool CxgIa32UseSse2;

void
CxIa32CpuInit(void);

#  endif // CxmCpuIa32
#endif // CxIa32_h
