#include "CxPpc.h"

#ifdef CxmCpuPpc
#include <sys/sysctl.h>

bool CxgPpcUseAltivec;

void
CxPpcCpuInit(void)
{
    int mib[2];
    int result, error;
    size_t len;

    mib[0] = CTL_HW;
    mib[1] = HW_VECTORUNIT;
    len = sizeof(result);
    error = sysctl(mib, 2, &result, &len, NULL, 0);
    if (error != 0)
    {
	CxgPpcUseAltivec = false;
	goto RETURN;
    }

    // The sysctl interface may return values greater than 1 in the future, so
    // don't check specifically for 1.
    if (result > 0)
    {
	CxgPpcUseAltivec = true;
    }
    else
    {
	CxgPpcUseAltivec = false;
    }

    RETURN:
    ;
}

#endif // CxmCpuPpc
