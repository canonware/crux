#ifndef CxPpc_h
#  define CxPpc_h

#  include "Cx.h"

#  ifdef CxmCpuPpc
#    define CxmCpuInit() CxPpcCpuInit()

extern bool CxgPpcUseAltivec;

void
CxPpcCpuInit(void);

#  endif // CxmCpuPpc
#endif // CxPpc_h
