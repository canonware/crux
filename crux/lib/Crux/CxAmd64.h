#ifndef CxAmd64_h
#  define CxAmd64_h

#  include "Cx.h"

#  ifdef CxmCpuAmd64
#    define CxmCpuInit() CxAmd64CpuInit()

extern bool CxgAmd64UseSse2;

void
CxAmd64CpuInit(void);

#  endif // CxmCpuAmd64
#endif // CxAmd64_h
